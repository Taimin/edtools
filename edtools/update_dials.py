import collections
import numpy as np
import scipy.ndimage as ndimage
from skimage.registration import phase_cross_correlation
from pathlib import Path
import shutil

from instamatic import config
from instamatic.formats import read_image
from instamatic.formats import write_adsc
from instamatic.formats import write_tiff
from instamatic.processing import apply_stretch_correction
from instamatic.processing.stretch_correction import affine_transform_ellipse_to_circle
from instamatic.tools import find_beam_center

from .utils import parse_args_for_fns


def update_dials(fn, wavelength, physical_pixelsize, pixelsize, exposure, phi, osc_angle, name,
             integrate=False, center=False, stretch=False):
    shutil.copyfile(fn, fn.with_name("dials_process.bat~"))

    mrc_folder = fn.parent.parent / 'RED'
    smv_folder = fn.parent
    img_name_list = mrc_folder.glob('*.mrc')
    img_list = []
    distance = (1 / wavelength) * (physical_pixelsize / pixelsize)

    for img_name in img_name_list:
        img, _ = read_mrc(img_name)
        img_list.append((img_name.name.split('.')[0], img))

    # TODO: The interpolation will blur out the diffraction pattern or make some artifacts...
    # Consider change to fourier shift
    for i, (img_name, img) in enumerate(img_list, 1): 
        if center:
            if i == 1:
                center_x, center_y = find_beam_center(img)
                template = img[int(round(center_x-16)):int(round(center_x+16)), 
                        int(round(center_y-16)):int(round(center_y+16))].copy()
                center_x_first, center_y_first = find_beam_center(template, sigma=5)
            else:
                center  = find_beam_center(img[int(round(center_x-16)):int(round(center_x+16)), 
                                            int(round(center_y-16)):int(round(center_y+16))], sigma=5)
                shift = (center_x_first-center[0], center_y_first-center[1])
                #shift, error, phasediff = phase_cross_correlation(template, center_area, upsample_factor=10)
                print(shift)
                img = ndimage.shift(img, shift, order=1, output=np.uint16, mode='nearest')

        if stretch:
            if stretch_cent_x is None or stretch_cent_y is None:
                if i == 1:
                    center_x, center_y = find_beam_center(img)
                template = img[int(round(center_x-16)):int(round(center_x+16)), 
                        int(round(center_y-16)):int(round(center_y+16))].copy()
                center_x_new, center_y_new = find_beam_center(template, sigma=5)
                update_cent_x, update_cent_y = center_x-16+center_x_new, center_y-16+center_y_new
                print(update_cent_x, update_cent_y)
            img = apply_stretch_correction(img, center=[update_cent_x, update_cent_y], azimuth=stretch_azimuth, amplitude=stretch_amplitude)

        header = collections.OrderedDict()
        header['HEADER_BYTES'] = 512
        header['DIM'] = 2
        header['BYTE_ORDER'] = 'little_endian'
        header['TYPE'] = 'unsigned_short'
        header['SIZE1'] = img.shape[1]
        header['SIZE2'] = img.shape[0]
        header['PIXEL_SIZE'] = physical_pixelsize
        header['BIN'] = '1x1'
        header['BIN_TYPE'] = 'HW'
        header['ADC'] = 'fast'
        header['CREV'] = 1
        header['BEAMLINE'] = name      # special ID for DIALS
        header['DETECTOR_SN'] = 901         # special ID for DIALS
        header['DATE'] = 0
        header['TIME'] = exposure
        header['DISTANCE'] = f'{distance:.4f}'
        header['TWOTHETA'] = 0.00
        header['PHI'] = f'{phi:.4f}'
        header['OSC_START'] = f'{phi:.4f}'
        header['OSC_RANGE'] = f'{osc_angle:.4f}'
        header['WAVELENGTH'] = f'{wavelength:.4f}'
        # reverse XY coordinates for XDS
        header['BEAM_CENTER_X'] = f'{center_x*physical_pixelsize:.4f}'
        header['BEAM_CENTER_Y'] = f'{center_y*physical_pixelsize:.4f}'
        header['DENZO_X_BEAM'] = f'{center_x*physical_pixelsize:.4f}'
        header['DENZO_Y_BEAM'] = f'{center_y*physical_pixelsize:.4f}'
        write_adsc(smv_folder/(img_name+'.img'), img, header)

def main():
    import argparse

    description = "Program to convert SMV file and update dials processing commands."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of dials_process.bat files or list of directories. If a list of directories is given "
                        "the program will find all dials_process.bat files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("-cent", "--center",
                        action="store", type=bool, dest="center",
                        help="Center ED patterns")

    parser.add_argument("-stre", "--stretch",
                        action="store", type=bool, dest="stretch",
                        help="Stretch ED patterns")

    parser.add_argument("-int", "--integrate",
                        action="store", type=bool, dest="integrate",
                        help="Integrate the data")

    parser.add_argument("-ref", "--refine",
                        action="store", type=bool, dest="refine",
                        help="Refine the data")

    parser.add_argument("-n", "--name",
                        action="store", type=str, dest="name",
                        help="Add camera name")

    parser.set_defaults(integrate=False,
                        center=False,
                        stretch=False,
                        refine=False,
                        name='Medipix3',)


    fns = parse_args_for_fns(fns, name="summary.txt", match=match)

    for fn in fns:
        lines = open(fn, "r", encoding = 'cp1252').readlines()
        for line in lines:
            if 'PixelSize' in line:
                pixelsize = float(line.split(' ')[1])
            elif 'Step size' in line:
                osc_angle = float(line.split(' ')[1])
            elif 'Physical' in line:
                physical_pixelsize = float(line.split(' ')[2])
            elif 'Wavelength' in line:
                wavelength = float(line.split(' ')[1])
            elif 'Exposure' in line:
                exposure = float(line.split(' ')[1])
        lines = open(fn.parent/'RED'/'1.ed3d', "r", encoding = 'cp1252').readlines()
        for line in lines:
            if '00001.mrc' in line:
                phi = float(line.split()[-1])

        print("\033[K", fn, end='\r')  # "\033[K" clears line
        update_dials(fn, wavelength, physical_pixelsize, pixelsize, exposure, phi, osc_angle, name,
                     integrate=integrate,
                     center=center,
                     stretch=stretch,)

    print(f"\033[KUpdated {len(fns)} files")


