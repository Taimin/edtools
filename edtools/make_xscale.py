from pathlib import Path, PurePosixPath
from sys import argv
import os, sys, time
from math import radians, cos
import numpy as np
import yaml
from collections import Counter
from .utils import space_group_lib

platform = sys.platform

spglib = space_group_lib()

threshold = 2.0


def parse_xds_ascii(fn):
    d = {"xds_ascii": fn.absolute()}
    with open(fn, "r") as f:
        for line in f:
            if not line.startswith("!"):
                break
            if "UNIT_CELL_CONSTANTS" in line:
                inp = line.split()
                cell = [float(val) for val in inp[-6:]]
            if "SPACE_GROUP_NUMBER" in line:
                inp = line.strip().split("=")
                spgr = int(inp[-1])

        d["space_group"] = spgr
        d["unit_cell"] = cell
    
    return d


def get_xds_ascii_names(lst):
    ret = []
    for d in lst:
        if "xds_ascii" in d:
            ret.append(d["xds_ascii"])
        else:
            ret.append(Path(d["directory"]) / "XDS_ASCII.HKL")
    return ret


def write_xscale_inp(fns, unit_cell, space_group, resolution, sigma, cc_half, reindex):
    cwd = Path(".").resolve()

    cell_str = " ".join((f"{val:.3f}" for val in unit_cell))

    print("\nUsing:")
    print(f"  SPACE_GROUP_NUMBER= {space_group}")
    print(f"  UNIT_CELL_CONSTANTS= {cell_str}")
    print()

    with open("XSCALE.INP", "w") as f:

        print("SNRC= 2", file=f)
        print("SAVE_CORRECTION_IMAGES= FALSE", file=f)  # prevent local directory being littered with .cbf files
        print(f"SPACE_GROUP_NUMBER= {space_group}", file=f)
        print(f"UNIT_CELL_CONSTANTS= {cell_str}", file=f)
        
        print(file=f)
        print("OUTPUT_FILE= MERGED.HKL", file=f)
        print(file=f)
        
        for i, fn in enumerate(fns):
            corr_fn = fn.parent/'CORRECT.LP'
            with open(corr_fn, 'r') as f_corr:
                lines = f_corr.readlines()
            if sigma is not None or cc_half is not None and len(resolution) == 0:
                resolution = [0, 0]
                
                for index, line in enumerate(lines):
                    if 'WILSON STATISTICS OF DATA SET  "XDS_ASCII.HKL"' in line:
                        line_num = index
                table = []
                for i in range(40):
                    try:
                        float(lines[line_num - i].split()[0])
                        table.append([float(a.strip('%').strip('*')) for a in lines[line_num - i].split()])
                    except:
                        pass
                resolution[0] = 80
                resolution[1] = table[0][0]
                for row in table:
                    if sigma is not None:
                        if row[8] < sigma:
                            resolution[1] = row[0]
                    if cc_half is not None:
                        if row[10] < cc_half:
                            resolution[1] = row[0] 

            reindex_matrix = None
            if reindex is not None:
                for index, line in enumerate(lines):
                    if '  LATTICE-  BRAVAIS-   QUALITY  UNIT CELL CONSTANTS (ANGSTROEM & DEGREES)    REINDEXING TRANSFORMATION' in line:
                        line_num = index
                for i in range(44):
                    split = lines[line_num+3+i].split()
                    if len(split) == 22:
                        if int(split[1]) == reindex:
                            reindex_matrix = '  '.join(split[10:])

            fn = fn.absolute()
            try:
                fn = fn.relative_to(cwd).as_posix()
            except ValueError:
                if platform == "win32":
                    drive = fn.drive
                    drive_letter = drive.lower()[0]
                    fn = fn.as_posix().replace(f"{drive}", f"/mnt/{drive_letter}")

            if reindex is None:
                print(f"    INPUT_FILE= {fn}", file=f)
                print(f"    INCLUDE_RESOLUTION_RANGE= {resolution[0]} {resolution[1]}", file=f)
            elif reindex is not None and reindex_matrix is not None:
                print(f"    INPUT_FILE= {fn}", file=f)
                print(f"    INCLUDE_RESOLUTION_RANGE= {resolution[0]} {resolution[1]}", file=f)
                print(f"    REIDX_ISET=  {reindex_matrix} ", file=f)
            print(file=f)

    print(f"Wrote file {f.name}")


def write_xdsconv_inp():
    with open("XDSCONV.INP", "w") as f:
        print(f"""
INPUT_FILE= MERGED.HKL
OUTPUT_FILE= shelx.hkl  SHELX    ! Warning: do _not_ name this file "temp.mtz" !
FRIEDEL'S_LAW= FALSE             ! default is FRIEDEL'S_LAW=TRUE""", file=f)

    print(f"Wrote file {f.name}")


def main():
    import argparse

    description = "Program to make an input file for XSCALE."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="Path to a cells.yaml / XDS_ASCII.HKL files")

    parser.add_argument("-s","--spgr",
                        action="store", type=int, dest="spgr",
                        help="Space group number (default: most common one)")

    parser.add_argument("-c","--cell",
                        action="store", type=float, nargs=6, dest="cell",
                        help="Override the unit cell parameters (default: mean unit cell)")

    parser.add_argument("-r","--resolution",
                        action="store", type=float, nargs=2, dest="resolution",
                        help="Override the resolution (default: mean unit cell)")

    parser.add_argument("-sig", "--sigma",
                        action="store", type=float, dest="sigma",
                        help="Use sigma to determine resolution")

    parser.add_argument("-cc", "--cc_half",
                        action="store", type=float, dest="cc_half",
                        help="Use CC_HALF to determine resolution")

    parser.add_argument("-rdx", "--reindex",
                        action="store", type=int, dest="reindex",
                        help="Reindex according to the lattice character")

    parser.set_defaults(cell=None,
                        spgr=None,
                        resolution=None,
                        sigma=None,
                        cc_half=None,
                        reindex=None)
    
    options = parser.parse_args()
    spgr = options.spgr
    cell = options.cell
    resolution = options.resolution
    args = options.args
    sigma = options.sigma
    cc_half = options.cc_half
    reindex = options.reindex
    
    if not args:  # attempt to populate args
        if os.path.exists("cells.yaml"):
            args = ["cells.yaml"]
        else:
            args = list(Path(".").glob("*XDS_ASCII.HKL"))
        
    if not args:
        exit()
    else:
        lst = []
        for arg in args:
            fn = Path(arg)
            extension = fn.suffix.lower()
            if extension == ".yaml":
                d = yaml.load(open(fn, "r"), Loader=yaml.Loader)
                lst.extend(d)
            if extension == ".hkl":
                lst.append(parse_xds_ascii(fn))

    print(f"Loaded {len(lst)} cells")

    fns = get_xds_ascii_names(lst)

    cells = np.array([d["unit_cell"] for d in lst])

    if resolution is None:
        resolution = []

    if not cell:
        cell = np.mean(cells, axis=0)

    if not spgr:
        c = Counter([spglib[d["space_group"]]["laue_symmetry"] for d in lst])
        for key, count in c.most_common(10):
            d = spglib[key]
            lattice = d["lattice"]
            laue_symm = d["laue_symmetry"]
            print(f"Lattice type `{lattice}` was found {count:3} times (space group number: {laue_symm})")
        spgr =  spglib[c.most_common()[0][0]]["laue_symmetry"]

    else:
        d = spglib[spgr]
        lattice = d["lattice"]
        laue_symm = d["laue_symmetry"]
        print(f"Lowest possible symmetry for {spgr} ({lattice}): {laue_symm}")

    write_xscale_inp(fns, unit_cell=cell, space_group=spgr, resolution=resolution, sigma=sigma, cc_half=cc_half, reindex=reindex)
    write_xdsconv_inp()


if __name__ == '__main__':
    main()
