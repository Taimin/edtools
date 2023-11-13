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
from instamatic.processing.ImgConversion import rotation_axis_to_xyz
from instamatic.tools import find_beam_center, find_subranges
from instamatic.processing.PETS2_template import PETS2_template

from .utils import parse_args_for_fns
import traceback


def update_dials(fn, wavelength, physical_pixelsize, pixelsize, exposure, phi, osc_angle, name, axis,
             write_smv=False, write_tif=False, integrate=False, center=False, stretch=False, refine=False, 
             symmetry=False, scale=False, merge=False, report=False, export=False, gain=1):
    if fn.exists():
        shutil.copyfile(fn, fn.with_name("dials_process.bat~"))

    mrc_folder = fn.parent.parent / 'RED'
    smv_folder = fn.parent
    tiff_folder = fn.parent.parent / 'tiff'
    img_name_list = mrc_folder.glob('*.mrc')
    img_list = []
    distance = (1 / wavelength) * (physical_pixelsize / pixelsize)

    for img_name in img_name_list:
        img, _ = read_image(img_name)
        img_list.append((img_name.name.split('.')[0], img))

    center_x_first = None
    pixel_num = 16

    # TODO: The interpolation will blur out the diffraction pattern or make some artifacts...
    # Consider change to fourier shift
    for i, (img_name, img) in enumerate(img_list, 1): 
        try:
            if i == 1:
                center_x, center_y = find_beam_center(img)
            if center:
                if i == 1 or center_x_first is None:
                    center_x, center_y = find_beam_center(img)
                    template = img[int(round(center_x-pixel_num)):int(round(center_x+pixel_num)), 
                            int(round(center_y-pixel_num)):int(round(center_y+pixel_num))].copy()
                    center_x_first, center_y_first = find_beam_center(template, sigma=5)
                    update_cent_x, update_cent_y = center_x-pixel_num+center_x_first, center_y-pixel_num+center_y_first
                else:
                    center_pos  = find_beam_center(img[int(round(center_x-pixel_num)):int(round(center_x+pixel_num)), 
                                                int(round(center_y-pixel_num)):int(round(center_y+pixel_num))], sigma=5)
                    shift = (center_x_first-center_pos[0], center_y_first-center_pos[1])
                    #shift, error, phasediff = phase_cross_correlation(template, center_area, upsample_factor=10)
                    print(shift)
                    img = ndimage.shift(img, shift, order=1, output=np.uint16, mode='nearest')
        except:
            traceback.print_exc()
            continue

        if stretch:
            if stretch_cent_x is None or stretch_cent_y is None:
                if i == 1:
                    center_x, center_y = find_beam_center(img)
                template = img[int(round(center_x-pixel_num)):int(round(center_x+pixel_num)), 
                        int(round(center_y-pixel_num)):int(round(center_y+pixel_num))].copy()
                center_x_new, center_y_new = find_beam_center(template, sigma=5)
                update_cent_x, update_cent_y = center_x-pixel_num+center_x_new, center_y-pixel_num+center_y_new
                print(update_cent_x, update_cent_y)
            img = apply_stretch_correction(img, center=[update_cent_x, update_cent_y], azimuth=stretch_azimuth, amplitude=stretch_amplitude)
        if write_smv:
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
            header['BEAM_CENTER_X'] = f'{update_cent_x*physical_pixelsize:.4f}'
            header['BEAM_CENTER_Y'] = f'{update_cent_y*physical_pixelsize:.4f}'
            header['DENZO_X_BEAM'] = f'{update_cent_x*physical_pixelsize:.4f}'
            header['DENZO_Y_BEAM'] = f'{update_cent_y*physical_pixelsize:.4f}'

            write_adsc(smv_folder/'data'/(img_name+'.img'), img, header)

        if write_tif:
            header = {}
            write_tiff(tiff_folder/(img_name+'.tiff'), img, header)

    img_name_list = smv_folder.rglob('*.img')
    img_list = []
    for img_name in img_name_list:
        img_list.append((img_name.name.split('.')[0]))
    scanranges = [int(num) for num in img_list]
    observed_range = set(scanranges)
    complete_range = set(range(min(observed_range), max(observed_range) + 1))
    missing_range = observed_range ^ complete_range
    scanranges = find_subranges(scanranges)
    scanrange = ' '.join(f'scan_range={i},{j}' for i, j in scanranges)
    excludeimages = ','.join(str(n) for n in missing_range)
    rot_x, rot_y, rot_z = rotation_axis_to_xyz(axis, setting='dials')
    with open(fn, 'w', encoding = 'cp1252') as f:
        print('@echo off', file=f)
        print('', file=f)
        print(f'set scan_range={scanrange}', file=f)
        print(f'set exclude_images=exclude_images={excludeimages}', file=f)
        print(f'set rotation_axis=geometry.goniometer.axes={rot_x:.4f},{rot_y:.4f},{rot_z:.4f}', file=f)
        print(f'call dials.import template=./data/#####.img %rotation_axis% lookup.dx=dx.pickle lookup.dy=dy.pickle panel.gain={gain}', file=f)
        print(f'call dials.find_spots imported.expt %scan_range% nproc=4', file=f)
        print(f'call dials.index imported.expt strong.refl max_lattices=3 refinement_protocol.n_macro_cycles=2 restrain.phil', file=f)
        if refine:
            print(f'call dials.refine_bravais_settings indexed.expt indexed.refl', file=f)
            print(f'call dials.refine indexed.expt indexed.refl', file=f)
        if integrate:
            print(f'call dials.integrate refined.expt refined.refl nproc=8', file=f)
        if symmetry:
            print(f'call dials.symmetry integrated.expt integrated.refl', file=f)
        if scale:
            print(f'call dials.scale symmetrized.expt symmetrized.refl', file=f)
        if merge:
            print(f'call dials.merge scaled.expt scaled.refl', file=f)
        if report:
            print(f'call dials.report scaled.expt scaled.refl', file=f)
        if export:
            print(f'call dials.export scaled.expt scaled.refl', file=f)
    
    if write_smv:
        empty = np.zeros_like(img)
        for n in missing_range:
            write_adsc(smv_folder/'data'/(f'{n:05d}'+'.img'), empty, header)

    if write_tif:
        omega = np.degrees(axis)
        # for pets, 0 <= omega <= 360
        if omega < 0:
            omega += 360
        elif omega > 360:
            omega -= 360
        s = PETS2_template.format(
            date="",
            wavelength=wavelength,
            pixelsize=pixelsize,
            precession_angle=0,
            omega=omega,
            indexing_method='diffspace',
            geometry='static',
            ellip_dist_amp=config.defaults.serial['ellip_dist_amp'],
            ellip_dist_phi1=config.defaults.serial['ellip_dist_phi'],
            ellip_dist_phi2=float(config.defaults.serial['ellip_dist_phi'])-45,
            G_gamma=config.defaults.serial['G_gamma'],
            refl_diameter=config.defaults.serial['refl_diameter'],
            psi=config.defaults.serial['psi'],
            I_sigma_search=config.defaults.serial['I_sigma_search'],
            I_sigma_camel=config.defaults.serial['I_sigma_camel'],
        )

        with open(tiff_folder.parent / 'pets.pts', 'w') as f:
            print(s, file=f)
            for name in img_list:
                fn = f'{name}.tiff'
                angle = phi + osc_angle * (int(name) - 1)
                print(f'{tiff_folder}/{fn} {angle:10.4f} 0.00', file=f)
            print('endimagelist', file=f)

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

    parser.add_argument("-m", "--match",
                        action="store", type=str, dest="match",
                        help="Include the dials_process.bat files only if they are in the given directories (i.e. --match DIALS_reprocessed)")

    parser.add_argument("-w", "--write_smv",
                        action="store", type=bool, dest="write_smv",
                        help="Write SMVs")

    parser.add_argument("-w_t", "--write_tif",
                        action="store", type=bool, dest="write_tif",
                        help="Write Tiffs")

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

    parser.add_argument("-sym", "--symmetry",
                        action="store", type=bool, dest="symmetry",
                        help="Check the symmetry of the data")

    parser.add_argument("-meg", "--merge",
                        action="store", type=bool, dest="merge",
                        help="Merge the data")

    parser.add_argument("-sca", "--scale",
                        action="store", type=bool, dest="scale",
                        help="Scale the data")

    parser.add_argument("-rep", "--report",
                        action="store", type=bool, dest="report",
                        help="Report the data")

    parser.add_argument("-exp", "--export",
                        action="store", type=bool, dest="export",
                        help="Export the data")

    parser.add_argument("-sk", "--skip",
                        action="store", type=int, dest="skip",
                        help="skip n frames")

    parser.add_argument("-inc", "--include_frames",
                        action="store", type=int, nargs=2, dest="include_frames",
                        help="Specify frame number to include frames to process")

    parser.add_argument("-g", "--gain",
                        action="store", type=int, nargs=, dest="gain",
                        help="Specify the gain value for the detector")

    parser.set_defaults(write_smv=False, write_tif=False, integrate=False, center=False, stretch=False, refine=False,
                        symmetry=False, scale=False, merge=False, report=False, export=False,
                        name='ADSC', skip=None, include_frames=None, gain=1)

    options = parser.parse_args()
    fns = options.args
    write_smv = options.write_smv
    write_tif = options.write_tif
    integrate = options.integrate
    center = options.center
    stretch = options.stretch
    refine = options.refine
    symmetry = options.symmetry
    scale = options.scale
    merge = options.merge
    report = options.report
    export = options.export
    name = options.name
    match = options.match
    skip = options.skip
    include_frames = options.include_frames
    gain = options.gain

    fns = parse_args_for_fns(fns, name="summary.txt", match=match)
    
    for fn in fns:
        lines = open(fn, "r", encoding = 'cp1252').readlines()
        for line in lines:
            if 'Pixelsize' in line:
                pixelsize = float(line.split()[1])
            elif 'Step size' in line:
                osc_angle = float(line.split()[2])
            elif 'Physical' in line:
                physical_pixelsize = float(line.split()[2])
            elif 'Wavelength' in line:
                wavelength = float(line.split()[1])
            elif 'Exposure' in line:
                exposure = float(line.split()[2])
            elif 'Rotation axis' in line:
                axis = float(line.split()[2])
        lines = open(fn.parent/'RED'/'1.ed3d', "r", encoding = 'cp1252').readlines()
        for line in lines:
            if '00001.mrc' in line:
                phi = float(line.split()[-1])
        (fn.parent/'SMV'/'data').mkdir(parents=True, exist_ok=True)
        (fn.parent/'tiff').mkdir(parents=True, exist_ok=True)
        print("\033[K", fn, end='\r')  # "\033[K" clears line
        update_dials(fn.parent/'SMV'/'dials_process.bat', wavelength, physical_pixelsize, 
                     pixelsize, exposure, phi, osc_angle, name, axis,
                     write_smv=write_smv, write_tif=write_tif, integrate=integrate, center=center, stretch=stretch, refine=refine,
                     symmetry=symmetry, scale=scale, merge=merge, report=report, export=export, gain=gain)

    print(f"\033[KUpdated {len(fns)} files")


