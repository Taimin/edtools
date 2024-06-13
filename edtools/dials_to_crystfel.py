import numpy as np
from pathlib import Path
import shutil
import os
import threading
import h5py
import subprocess
import json
import traceback
import yaml
import concurrent.futures

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

def process_data(index, fn, split, write_h5, lock, files, d_min, reindex=False, refine=False, integrate=False, \
                integrate_sweep=False, scale_sweep=False, file_exists=False, space_group=None, write_sol=False, merge=False,
                single_crystal=True, gain=1):
    try:
        print(f'Start processing crystal number {index}.')
        drc = fn.parent/'SMV'
        cwd_smv = str(drc)
        if not (drc / 'indexed.expt').is_file() and not (drc / 'indexed.refl').is_file():
            print(f'indexed.expt or indexed.refl file does not exist for crystal number {index}.')
            return -1

        #if reindex:
        #    print(f'Start reindex {fn}')
        #    if space_group is None:
        #        pass
        #    else:
        #        cmd = f'dials.reindex.bat indexed.expt indexed.refl space_group={space_group}'
        #    try:
        #        print(cmd)
        #        p = subprocess.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
        #        p.communicate()
        #    except Exception as e:
        #        print("ERROR in subprocess call:", e)
        if refine:
            print(f'Start refine {fn}')
            if reindex:
                if (drc / 'reindexed.expt').is_file():
                    cmd = 'dials.refine.bat reindexed.expt reindexed.refl detector.fix=distance unit_cell.force_static=True nproc=2'
                else:
                    cmd = 'dials.refine.bat indexed.expt indexed.refl detector.fix=distance unit_cell.force_static=True nproc=2'
            else:
                cmd = 'dials.refine.bat indexed.expt indexed.refl detector.fix=distance unit_cell.force_static=True nproc=2'
            try:
                print(cmd)
                p = subprocess.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
                p.communicate()
            except Exception as e:
                print("ERROR in subprocess call:", e)
        if split:
            print(f'Start split {fn}')
            if refine:
                if (drc / 'refined.expt').is_file():
                    cmd = f'dials.sequence_to_stills.bat refined.expt refined.refl detector.gain={gain}'
                else:
                    cmd = f'dials.sequence_to_stills.bat indexed.expt indexed.refl detector.gain={gain}'
            elif reindex:
                if (drc / 'reindexed.expt').is_file():
                    cmd = f'dials.sequence_to_stills.bat reindexed.expt reindexed.refl detector.gain={gain}'
                else:
                    cmd = f'dials.sequence_to_stills.bat indexed.expt indexed.refl detector.gain={gain}'
            else:
                cmd = f'dials.sequence_to_stills.bat indexed.expt indexed.refl detector.gain={gain}'
            try:
                print(cmd)
                p = subprocess.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
                p.communicate()
            except Exception as e:
                print("ERROR in subprocess call:", e)
        with open(drc / 'stills.expt', 'r') as f:
            d = json.load(f)
            crystals = d['crystal']
        # convert every folder to an h5 file, file name is the relative path, save it to the h5 folder in CWD
        # /entry/data/raw_counts   /entry/data/corrected
        with open(drc/'dials.import.log', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('  template ='):
                    folder = line.split('/')[1]

        img_list = list((drc/folder).glob('*.img'))
        center_list = []
        imgs = []
        frame_list = []
        event_list = []
        if single_crystal:
            if len(crystals) > len(img_list):
                print(f'There are more than one crystal under the beam for crystal number {index}')
                return -1

        if integrate_sweep:
            if refine:
                cmd = 'dials.integrate.bat refined.expt refined.refl nproc=2'
            elif reindex:
                cmd = 'dials.integrate.bat reindexed.expt reindexed.refl nproc=2'
            else:
                cmd = 'dials.integrate.bat indexed.expt indexed.refl nproc=2'
            try:
                p = subprocess.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
                p.communicate()
            except Exception as e:
                print("ERROR in subprocess call:", e)

        if integrate:
            cmd = f'dials.ssx_integrate.bat stills.expt stills.refl prediction.d_min={d_min} mosaicity_max_limit=0.2 ellipsoid.unit_cell.fixed=True min_n_reflections=5 nproc=2'
            try:
                p = subprocess.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
                p.communicate()
            except Exception as e:
                print("ERROR in subprocess call:", e)
        if merge:
            if not (drc/'integrated_1.refl').is_file():
                print(f"{drc/'integrated_1.refl'} file does not exist.")
                return -1
            if space_group is not None:
                cmd = f'dials.reindex.bat integrated_1.refl integrated_1.expt space_group={space_group} output.experiments=re_integrated_1.expt output.reflections=re_integrated_1.refl'
                try:
                    print(f'Start changing space group for {fn}')
                    p = subprocess.Popen(cmd, cwd=cwd_smv, stdout=DEVNULL)
                    p.communicate()
                    target = list(drc.glob('re_integrated_1.expt'))[0]
                except Exception as e:
                    print("ERROR in subprocess call:", e)
            elif space_group is None:
                target = list(drc.glob('integrated_1.expt'))[0]

            if scale_sweep:
                target_1 = list(drc.glob('integrated.expt'))[0]
            if scale_sweep:
                files.append([cwd_smv+'/'+target.name, cwd_smv+'/'+target.name.split('.')[0]+'.refl', cwd_smv+'/'+target_1.name, cwd_smv+'/'+target_1.name.split('.')[0]+'.refl'])
            else:
                files.append([cwd_smv+'/'+target.name, cwd_smv+'/'+target.name.split('.')[0]+'.refl'])

        h5_filename = '_'.join(os.path.relpath(img_list[0].parent, CWD).split(os.sep))
        h5_filename = h5_filename + '.h5'
        if write_h5:
            for index, img_file in enumerate(img_list):
                img, h = read_image(img_file)
                imgs.append(img)
                center = [float(h['BEAM_CENTER_Y']), float(h['BEAM_CENTER_X'])]
                center_list.append(center)
                frame_list.append(index)
                event_list.append(f'entry//{index}')

            file_list = ['./h5/'+h5_filename] * len(img_list)
            file_list = [s.encode('utf-8') for s in file_list]
            imgs = np.array(imgs)
            event_list = [s.encode('utf-8') for s in event_list]
            center_X_mm, center_Y_mm = zip(*center_list)
            center_X_mm = np.array(list(center_X_mm))
            center_Y_mm = np.array(list(center_Y_mm))
            center_X = center_X_mm / float(h['PIXEL_SIZE'])
            center_Y = center_Y_mm / float(h['PIXEL_SIZE'])
            
            with h5py.File(CWD/'h5'/h5_filename, 'w') as f:
                dt = h5py.special_dtype(vlen=bytes)
                f.create_dataset('/entry/data/raw_counts', data=imgs)
                f.create_dataset('/entry/data/index', data=frame_list)
                f.create_dataset('/entry/shots/det_shift_x_mm', data=center_X_mm, dtype=np.float32)
                f.create_dataset('/entry/shots/det_shift_y_mm', data=center_Y_mm, dtype=np.float32)
                f.create_dataset('/entry/shots/center_x', data=center_X, dtype=np.float32)
                f.create_dataset('/entry/shots/center_y', data=center_Y, dtype=np.float32)
                f.create_dataset('/entry/shots/com_x', data=center_X, dtype=np.float32)
                f.create_dataset('/entry/shots/com_y', data=center_Y, dtype=np.float32)
                f.create_dataset('/entry/shots/Event', data=event_list, dtype=dt)
                f.create_dataset('/entry/shots/frame', data=frame_list)
                f.create_dataset('/entry/shots/file', data=file_list, dtype=dt)
        # generate a .lst file in the directory, append relative path 

        with lock:
            if not file_exists:
                with open(CWD / "files.lst", "a") as f:
                    print('./h5/'+h5_filename, file=f)
        # make .sol file, append reciprocal vector: inverse and transpose, in nm-1
        # append center: shift from image center in mm
        if write_sol:
            with lock:
                with open(CWD / "indexed.sol", "a") as f:
                    h5_filename = h5_filename.split('.')[0] + '_hit.h5'
                    for index, crystal in enumerate(crystals):
                        real_sp_matrix = np.array([crystal['real_space_a'], crystal['real_space_b'], crystal['real_space_c']])
                        real_sp_matrix = real_sp_matrix / 10
                        reciprocal_sp_matrix = np.linalg.inv(real_sp_matrix).T
                        reciprocal_sp_matrix[:, 0] = -reciprocal_sp_matrix[:, 0]
                        reciprocal_sp_matrix = reciprocal_sp_matrix.flatten().tolist()
                        reciprocal_sp_matrix = [f'{num:.7f}' for num in reciprocal_sp_matrix]
                        sp = crystal['space_group_hall_symbol'].replace(' ', '')
                        for num, tmp in spglib.items():
                            if tmp['name'] == sp:
                                sp_num = num
                        lattice = spglib[sp_num]['lattice']
                        sym = lattice_type_sym(lattice)
                        number = index % len(img_list)
                        print(f"./h5/{h5_filename} entry//{number} {' '.join(reciprocal_sp_matrix)} 0 0 {sym}", file=f)
    except:
        traceback.print_exc()

def run_parallel(fns, split, write_h5, lock, d_min, thresh, reindex=False, refine=False, merge=False, integrate=False, \
                integrate_sweep=False, scale_sweep=False, file_exists=False, space_group=None, write_sol=False, single_crystal=True, gain=1,
                reference=None):
    futures = []
    FILES = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        for index, fn in enumerate(fns):
            futures.append(executor.submit(process_data, index, fn, split, write_h5, lock, FILES, d_min, reindex, refine, integrate,\
                                         integrate_sweep, scale_sweep, file_exists, space_group, write_sol, merge, single_crystal, gain))
    concurrent.futures.wait(futures, return_when=concurrent.futures.ALL_COMPLETED)

    if merge:
        FILES = list(zip(*FILES))
        experiments = ' '.join(FILES[0])
        reflections = ' '.join(FILES[1])
        if scale_sweep:
            experiments_sweep = ' '.join(FILES[2])
            reflections_sweep = ' '.join(FILES[3])
            cmd = f'dials.scale.bat {experiments_sweep} {reflections_sweep} d_min={d_min} output.unmerged_mtz=scaled.mtz'
            try:
                p = subprocess.Popen(cmd, cwd=CWD)
                p.communicate()
            except Exception as e:
                print("ERROR in subprocess call:", e)

        with open(CWD/'merge.phil', 'w') as f:
            print("input {", file=f)
            for experiment, reflection in zip(FILES[0], FILES[1]):
                print(f'  experiments={experiment}', file=f)
                print(f'  reflections={reflection}', file=f)
            print("}", file=f)
            print('multiprocessing.nproc = 8', file=f)
            print('clustering {', file=f)
            print(f'  absolute_angle_tolerance=None', file=f)
            print(f'  absolute_length_tolerance=None', file=f)
            print('}', file=f)
            print(f'partiality_threshold = {thresh}', file=f)
            print(f'd_min = {d_min}', file=f)
            print('symmetry {', file=f)
            print(f'  lattice_symmetry_max_delta=0', file=f)
            if space_group is not None: 
                print(f'  space_group = {space_group}', file=f)
            print('}', file=f)
            if scale_sweep:
                reference = 'scaled.mtz'
                print(f'reference = {reference}', file=f)

        cmd = f'xia2.ssx_reduce.bat merge.phil'
        try:
            p = subprocess.Popen(cmd, cwd=CWD)
            p.communicate()
        except Exception as e:
            print("ERROR in subprocess call:", e)

        if reference is None:
            cmd = f'dials.export.bat DataFiles/scaled.expt DataFiles/scaled.refl format=shelx partiality_threshold={thresh}'
        else:
            cmd = f'dials.export.bat DataFiles/scaled_batch1.expt DataFiles/scaled_batch1.refl format=shelx partiality_threshold={thresh}'
        try:
            p = subprocess.Popen(cmd, cwd=CWD)
            p.communicate()
        except Exception as e:
            print("ERROR in subprocess call:", e)

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

    parser.add_argument("-w", "--write_h5",
                        action="store", type=bool, dest="write_h5",
                        help="Split sequences into stills")

    parser.add_argument("-re", "--reindex",
                        action="store", type=bool, dest="reindex",
                        help="Convert reindex results instead of index results")

    parser.add_argument("-f", "--file_exists",
                        action="store", type=bool, dest="file_exists",
                        help="A input file that list all the directories")

    parser.add_argument("-int", "--integrate",
                        action="store", type=bool, dest="integrate",
                        help="Integrate the splitted results using DIALS")

    parser.add_argument("-int_s", "--integrate_sweep",
                        action="store", type=bool, dest="integrate_sweep",
                        help="Integrate the dataset as a sweep")

    parser.add_argument("-sca_s", "--scale_sweep",
                        action="store", type=bool, dest="scale_sweep",
                        help="Scale a sweep dataset")

    parser.add_argument("-mg", "--merge",
                        action="store", type=bool, dest="merge",
                        help="Merge the integrated results from DIALS")

    parser.add_argument("-r", "--refine",
                        action="store", type=bool, dest="refine",
                        help="Whether use DIALS to refine the indexed file before reindex")

    parser.add_argument("-d", "--d_min",
                        action="store", type=float, dest="d_min",
                        help="Minimum distance for ssx_reduce.")

    parser.add_argument("-t", "--thresh",
                        action="store", type=float, dest="thresh",
                        help="Partiality threshold for dials.export.")

    parser.add_argument("-sg", "--space_group",
                        action="store", type=str, dest="space_group",
                        help="Set the space group for data reduction.")

    parser.add_argument("-w_sol", "--write_sol",
                        action="store", type=bool, dest="write_sol",
                        help="Write the sol file.")

    parser.add_argument("-sc", "--single_crystal",
                        action="store", type=bool, dest="single_crystal",
                        help="Only use single crystals.")

    parser.add_argument("-g", "--gain",
                        action="store", type=float, dest="gain",
                        help="Gain used in integration.")

    parser.add_argument("-ref", "--reference",
                        action="store", type=str, dest="reference",
                        help="Reference during scaling.")

    parser.set_defaults(split=False, write_h5=False, reindex=False, refine=False, d_min=0.8, thresh=0.6, 
                        file_exists=False, space_group=None, write_sol=False, single_crystal=False, gain=1, reference=None)

    options = parser.parse_args()
    args = options.args
    match = options.match
    split = options.split
    write_h5 = options.write_h5
    reindex = options.reindex
    file_exists = options.file_exists
    refine = options.refine
    d_min = options.d_min
    thresh = options.thresh
    merge = options.merge
    integrate = options.integrate
    integrate_sweep = options.integrate_sweep
    space_group = options.space_group
    write_sol = options.write_sol
    scale_sweep = options.scale_sweep
    single_crystal = options.single_crystal
    gain = options.gain
    reference = options.reference

    if args:
        fns = []
        for arg in args:
            if arg.split('.')[-1] == "yaml":
                ds = yaml.load(open(arg, "r"), Loader=yaml.Loader)
                for d in ds:
                    fns.append(Path(d['directory']))
            elif arg.split('.')[-1] == "lst":
                with open(arg, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.split('/')[-1].split('_')
                        folder = ['_'.join(line[0:3]), '_'.join(line[3:5])]
                        folder = '/'.join(folder)
                        file = Path('./' + folder) / "summary.txt"
                        fns.append(file)
        fns = [fn.resolve() for fn in fns]
    else:
        fns = parse_args_for_fns(args, name="summary.txt", match=match)

    (CWD/'h5').mkdir(parents=True, exist_ok=True)

    if not file_exists:
        with open(CWD / "files.lst", "w") as f:
            pass

    with open(CWD / "indexed.sol", "w") as f:
        pass
    lock = threading.Lock()
    run_parallel(fns, split, write_h5, lock, d_min, thresh, reindex, refine, merge, integrate, integrate_sweep, \
                scale_sweep, file_exists, space_group, write_sol, single_crystal, gain, reference)
    print(f"\033[KUpdated {len(fns)} files")
