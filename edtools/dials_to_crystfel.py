import numpy as np
from pathlib import Path
import shutil

from instamatic import config
from instamatic.formats import read_image
from instamatic.formats import write_hdf5

from .utils import parse_args_for_fns

DEVNULL = open(os.devnull, 'w')
CWD = Path(os.getcwd())
rlock = threading.RLock()

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

    options = parser.parse_args()
    fns = options.args

    fns = parse_args_for_fns(fns, name="summary.txt", match=match)

    (CWD/'h5').mkdir(parents=True, exist_ok=True)

    with open(CWD / "files.lst", "w") as f:
        pass

    #with open(CWD / "indexed.sol", "w") as f:
    #    pass

    for fn in fns:
        drc = fn.parent/'SMV'
        cwd = str(drc)
        cmd = 'dials.sequence_to_stills'
        print("\033[K", fn, end='\r')  # "\033[K" clears line
        try:
            p = sp.Popen(cmd, cwd=cwd, stdout=DEVNULL)
            p.wait()
        except Exception as e:
            print("ERROR in subprocess call:", e)

        # convert every folder to an h5 file, file name is the relative path, save it to the h5 folder in CWD
        # /entry/data/raw_counts   /entry/data/corrected
        img_list = list(drc.rglob('*.img'))
        center_list = []
        imgs = []
        header = {}
        for img_file in img_list:
            img, h = read_image(img_file)
            imgs.append(img)
            center = [h['BEAM_CENTER_Y'], h['BEAM_CENTER_X']]
            center_list.append(center)

        imgs = np.array(imgs)
        center_X, center_Y = list(zip(center))
        header['det_shift_x_mm'] = center_X
        header['det_shift_y_mm'] = center_Y
        h5_filename = '_'.join(os.path.relpath(img_file.parent, CWD).split(os.sep))
        h5_filename = h5_filename + '.h5'
        with h5py.File(CWD/'h5'/h5_filename, 'w') as f:
            h5data = f.create_dataset('/entry/data', data=imgs)
            h5data.attrs.update(header)

        # generate a .lst file in the directory, append relative path 
        with open(CWD / "files.lst", "a") as f:
            print(CWD/'h5'/h5_filename, file=f)

        # make .sol file, append reciprocal vector: inverse and transpose, in nm-1
        # append center: shift from image center in mm
        #with open(CWD / "indexed.sol", "a") as f:
        #    pass



    print(f"\033[KUpdated {len(fns)} files")