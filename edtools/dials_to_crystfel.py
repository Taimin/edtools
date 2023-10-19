import numpy as np
from pathlib import Path
import shutil
import os
import threading
import h5py
import subprocess as sp
import json

from instamatic import config
from instamatic.formats import read_image
from instamatic.formats import write_hdf5

from .utils import parse_args_for_fns
from .utils import space_group_lib

DEVNULL = open(os.devnull, 'w')
CWD = Path(os.getcwd())
rlock = threading.RLock()
spglib = space_group_lib()

def lattice_type_sym(lattice, unique_axis='c'):
    if lattice[0] == 'a':
        return lattice
    elif lattice[0] == 'm':
        return lattice + unique_axis
    elif lattice[0] == 'o':
        return lattice
    elif lattice[0] == 't':
        return lattice + unique_axis
    elif lattice[0] == 'cubic':
        return lattice
    elif lattice[0] == 'h':
        return lattice + unique_axis
    elif lattice[0] == 'r':
        return lattice
    else:
        warn('Invalid lattice type {}'.format(lattice))
        return 'invalid'

def main():
    import argparse

    description = "Program to convert results from DIALS to files required by Crystfel."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of dials_process.bat files or list of directories. If a list of directories is given "
                        "the program will find all dials_process.bat files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("-m", "--match",
                        action="store", type=str, dest="match",
                        help="Include files only if they are in the given directories (i.e. --match DIALS_reprocessed)")

    parser.add_argument("-sp", "--split",
                        action="store", type=bool, dest="split",
                        help="Split sequences into stills")

    parser.set_defaults(split=False)

    options = parser.parse_args()
    fns = options.args
    match = options.match
    split = options.split

    fns = parse_args_for_fns(fns, name="summary.txt", match=match)

    (CWD/'h5').mkdir(parents=True, exist_ok=True)

    with open(CWD / "files.lst", "w") as f:
        pass

    #with open(CWD / "indexed.sol", "w") as f:
    #    pass

    for fn in fns:
        drc = fn.parent/'SMV'
        cwd_smv = str(drc)
        if split:
            cmd = 'dials.sequence_to_stills.bat indexed.expt indexed.refl'
            try:
                p = sp.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
                p.wait()
            except Exception as e:
                print("ERROR in subprocess call:", e)

        # convert every folder to an h5 file, file name is the relative path, save it to the h5 folder in CWD
        # /entry/data/raw_counts   /entry/data/corrected
        img_list = list(drc.rglob('*.img'))
        center_list = []
        imgs = []
        frame_list = []
        event_list = []
        for index, img_file in enumerate(img_list):
            img, h = read_image(img_file)
            imgs.append(img)
            center = [float(h['BEAM_CENTER_Y']), float(h['BEAM_CENTER_X'])]
            center_list.append(center)
            frame_list.append(index)
            event_list.append(f'entry//{index}')

        imgs = np.array(imgs)
        event_list = [s.encode('utf-8') for s in event_list]
        center_X_mm, center_Y_mm = zip(*center_list)
        center_X_mm = np.array(list(center_X_mm))
        center_Y_mm = np.array(list(center_Y_mm))
        center_X = center_X_mm / float(h['PIXEL_SIZE'])
        center_Y = center_Y_mm / float(h['PIXEL_SIZE'])
        h5_filename = '_'.join(os.path.relpath(img_file.parent, CWD).split(os.sep))
        h5_filename = h5_filename + '.h5'
        with h5py.File(CWD/'h5'/h5_filename, 'w') as f:
            dt = h5py.special_dtype(vlen=bytes)
            f.create_dataset('/entry/data/raw_counts', data=imgs)
            f.create_dataset('/entry/data/index', data=frame_list)
            f.create_dataset('/entry/shots/det_shift_x_mm', data=center_X_mm)
            f.create_dataset('/entry/shots/det_shift_y_mm', data=center_Y_mm)
            f.create_dataset('/entry/shots/center_x', data=center_X)
            f.create_dataset('/entry/shots/center_y', data=center_Y)
            f.create_dataset('/entry/shots/com_x', data=center_X)
            f.create_dataset('/entry/shots/com_y', data=center_Y)
            f.create_dataset('/entry/shots/Event', data=event_list, dtype=dt)
            f.create_dataset('/entry/shots/frame', data=frame_list)

        # generate a .lst file in the directory, append relative path 
        with open(CWD / "files.lst", "a") as f:
            print('./h5/'+h5_filename, file=f)

        # make .sol file, append reciprocal vector: inverse and transpose, in nm-1
        # append center: shift from image center in mm
        with open(drc / 'stills.expt', 'r') as f:
            d = json.load(f)
            crystals = d['crystal']

        with open(CWD / "indexed.sol", "a") as f:
            for index, crystal in enumerate(crystals):
                real_sp_matrix = np.array([crystal['real_space_a'], crystal['real_space_b'], crystal['real_space_c']])
                real_sp_matrix = real_sp_matrix / 10
                reciprocal_sp_matrix = np.linalg.inv(real_sp_matrix).T
                reciprocal_sp_matrix = reciprocal_sp_matrix.flatten().tolist()
                reciprocal_sp_matrix = [str(num) for num in reciprocal_sp_matrix]
                sp = crystal['space_group_hall_symbol'].replace(' ', '')
                for num, tmp in spglib.items():
                    if tmp['name'] == sp:
                        sp_num = num
                lattice = spglib[sp_num]['lattice']
                sym = lattice_type_sym(lattice)
                print(f"h5/{h5_filename} entry//{index} {' '.join(reciprocal_sp_matrix)} 0 0 {sym}", file=f)

    print(f"\033[KUpdated {len(fns)} files")