from pathlib import Path
import os
import time
import shutil
import re
from .utils import volume, parse_args_for_fns
from .utils import space_group_lib
import copy

spglib = space_group_lib()


class dials_parser:
    """docstring for dials_parser"""
    def __init__(self, filename, type='index'):
        self.ios_threshold = 0.8

        self.filename = Path(filename).resolve()
        self.d = self.parse(type)

    def parse(self, type):
        fn = self.filename
        with open(fn, "r") as f:
            lines = f.readlines()

        crystals = []
        line_start = []
        line_end = []

        d = {}
        d_list = []

        for index, line in enumerate(lines):
            if line.startswith("Indexed crystal models:"):
                line_start.append(index)
            elif line.startswith("Starting refinement (macro-cycle 1)"):
                line_end.append(index)

        for start, end in zip(line_start, line_end):
            if end > start:
                crystals.append(lines[start:end])
            else:
                print('Check the index.log file')
                return -1
                
        for crystal in crystals:
            cell, spgr = None, None
            for index, line in enumerate(crystal):
                if line.startswith("model"):
                    model_num = int(line.split()[1])
                elif line.startswith("    Unit cell"):
                    line = re.sub(r'\([^)]*\)', '', line)
                    line = line.strip("\n").split()[2:8]
                    line = [l.strip(',') for l in line]
                    cell = list(map(float, line))
                elif line.startswith("    Space group"):
                    spgr = ''.join(line.strip("\n").split()[2:])
                elif line.startswith("|   Imageset"):
                    infos = [info.strip() for info in crystal[index+2].split('|')]
                    d['indexed'] = int(infos[2])
                    d['unindexed'] = int(infos[3])
            d["volume"] = volume(cell)
            d["cell"] = cell
            d["spgr"] = spgr
            d["model_num"] = model_num
            d["fn"] = fn
            d_list.append(copy.deepcopy(d))
        
        return d_list

    @staticmethod
    def info_header(hline=True):
        s  = "   #   dmax  dmin    ntot   nuniq   compl   i/sig   rmeas CC(1/2)     ISa   B(ov)\n"
        if hline:
            s += "---------------------------------------------------------------------------------\n"
        return s

    def print_filename(self):
        print("#", self.filename)

    def cell_info(self, sequence=0):
        d = self.d
        i = sequence
        s_list = []
        for crystal in d:
            fn = self.filename
            model_num = crystal["model_num"]
            s = f"{i: 4d}: {fn.parents[0]} # {model_num} # {time.ctime(os.path.getmtime(fn))}\n"
            s += "Spgr {} - Cell {:10.2f}{:10.2f}{:10.2f}{:10.2f}{:10.2f}{:10.2f} - Vol {:10.2f} Indexed: {} Unindexed: {}\n"\
                    .format(crystal["spgr"], *crystal["cell"], crystal["volume"], crystal["indexed"], crystal["unindexed"])
            s_list.append(s)
        return s_list

    def integration_info(self, sequence=0, outer_shell=True, filename=False):
        d = self.d
        k = sequence

        s = ""

        if k == 0:
            s += self.info_header()

        dmax, dmin = d["res_range"]

        s += "{k: 4d} {dmax: 6.2f}{dmin: 6.2f}{ntot: 8d}{nuniq: 8d}{completeness: 8.1f}{ios: 8.2f}{rmeas: 8.1f}{cchalf: 8.1f}{ISa: 8.2f}{Boverall: 8.2f}".format(
        k=k, dmax=dmax, dmin=dmin, **d["total"], **d)

        if filename:
            s += f"  # {d['fn']}\n"
        else:
            s += "\n"

        if outer_shell:
            outer = d["outer"]
            dmax_sh, dmin_sh = d["outer_shell"]
            s +="   - {dmax: 6.2f}{dmin: 6.2f}{ntot: 8d}{nuniq: 8d}{completeness: 8.1f}{ios: 8.2f}{rmeas: 8.1f}{cchalf: 8.1f}\n".format(
                k=k, dmax=dmax_sh, dmin=dmin_sh, **d[outer])

        return s


def cells_to_excel(ps, fn="cells.xlsx"):
    """Takes a list of `dials_parser` instances and writes the cell
    parameters to an excel file `cells.xlsx`.
    """
    d = {}
    for i, p in enumerate(ps):
        i += 1
        d[i] = p.cell_as_dict()

    import pandas as pd
    df = pd.DataFrame(d).T
    df = df["spgr a b c al be ga volume rotation_angle total_completeness file".split()]
    if not os.path.exists(fn):
        df.to_excel(fn)
    else:
        """To address that cells.xlsx does not overwrite"""
        os.remove(fn)
        df.to_excel(fn)

    print(f"Wrote {i} cells to file {fn}")


def cells_to_cellparm(ps):
    """Takes a list of `dials_parser` instances and writes the cell
    parameters to an instruction file `CELLPARM.INP` for the program
    `cellparm`.
    """
    fn = "CELLPARM.INP"
    # write cellparm input file
    with open(fn, "w") as f:
        for i, p in enumerate(ps):
            i += 1
            fn = p.filename
            # cell = p.unit_cell
            cell = p.d["raw_cell"]
            ntot = p.d["total"]["ntot"]
            print(f"! {i: 3d} from {fn}", file=f)
            print("UNIT_CELL_CONSTANTS= {:10.2f}{:10.2f}{:10.2f}{:10.2f}{:10.2f}{:10.2f} WEIGHT= {ntot}".format(*cell, ntot=ntot), file=f)

    print(f"Wrote file {fn}")


def cells_to_yaml(ps, fn="cells.yaml"):
    import yaml
    ds = []
    for i, p in enumerate(ps):
        i += 1
        d = {}
        d["directory"] = str(p.filename.parent)
        d["number"] = i
        d["unit_cell"] = p.d["cell"]
        d["raw_unit_cell"] = p.d["raw_cell"]
        d["space_group"] = p.d["spgr"]
        d["weight"] = p.d["total"]["ntot"]
        ds.append(d)

    yaml.dump(ds, open(fn, "w"))

    print(f"Wrote {i} cells to file {fn}")


def gather_DIALS_refl(ps, min_completeness=10.0, min_cchalf=90.0, gather=False):
    """Takes a list of `dials_parser` instances and gathers the
    corresponding relf files into the current directory.
    The data source and numbering scheme is summarized in the file `filelist.txt`.
    """
    fn = "filelist.txt"

    # gather DIALS_ascii and prepare filelist
    n = 0
    with open(fn, "w") as f:
        for i, p in enumerate(ps):
            i += 1

            completeness = p.d["total"]["completeness"]
            cchalf = p.d["total"]["cchalf"]

            if cchalf < min_cchalf:
                continue

            if completeness < min_completeness:
                continue

            src = p.filename.with_name("DIALS_ASCII.HKL")
            dst = f"{i:02d}_DIALS_ASCII.HKL"
            if gather:
                shutil.copy(src, dst)
                ascii_name = dst
            else:
                ascii_name = src

            dmax, dmin = p.d["res_range"]
            print(f" {i: 3d} {ascii_name} {dmax:8.2f} {dmin:8.2f}  # {p.filename}", file=f)
            n += 1

    print(f"Wrote {n} entries to file {fn} (completeness > {min_completeness}%, CC(1/2) > {min_cchalf}%)")


def lattice_to_space_group(lattice):
    return { 'aP':  1, 'mP':  3, 'mC':  5, 'mI':  5,
             'oP': 16, 'oC': 21, 'oI': 23, 'oF': 22,
             'tP': 75, 'tI': 79, 'hP':143, 'hR':146,
             'cP':195, 'cF':196, 'cI':197 }[lattice]


def evaluate_symmetry(ps):
    from collections import Counter

    c_score = Counter()
    c_freq = Counter()

    for p in ps:
        spgr = p.d["spgr"]
        weight = p.d["total"]["ntot"]
        d = spglib[spgr]
        lattice = d["lattice"]
        c_score[lattice] += weight
        c_freq[lattice] += 1

    print("\nMost likely lattice types:")
    n = 1
    for lattice, score in c_score.most_common(100):
        count = c_freq[lattice]
        spgr = lattice_to_space_group(lattice)
        print(f"{n:3} Lattice type `{lattice}` (spgr: {spgr:3}) was found {count:3} times (score: {score:7})")
        n += 1

    return lattice_to_space_group(c_score.most_common()[0][0])

def parse_xparm_for_uc(fn):
    with open(fn, "r") as f:
        f.readline()
        f.readline()
        f.readline()

        line = f.readline().split()
        uc = line[1:]
        uc = [float(i) for i in uc]
        return uc

def cells_to_yaml_xparm(uc, fn="cells_xparm.yaml"):
    import yaml
    ds = []
    for i, p in enumerate(uc):
        i += 1
        d = {}

        d["directory"] = str(Path(p[1]).parent.resolve())

        """get rotation range from dials.import.log"""
        DIALSinp = Path(p[1]).parent / "DIALS.INP"
        with open(DIALSinp, "r") as f:
            for line in f:
                if line.startswith("DATA_RANGE="):
                    datarange = list(map(float, line.strip("\n").split()[1:]))
                elif line.startswith("OSCILLATION_RANGE="):
                    osc_angle = float(line.strip("\n").split()[1])

        rr = osc_angle * (datarange[1] - datarange[0] + 1)

        d["number"] = i
        d["rotation_range"] = rr
        d["raw_unit_cell"] = p[0]
        d["space_group"] = "P1"
        d["weight"] = 1
        ds.append(d)

    yaml.dump(ds, open(fn, "w"))

    print(f"Wrote {i} cells to file {fn}")

def main():
    import argparse

    description = "Program to consolidate data from a large series of data sets from a serial crystallography experiment."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of log files or list of directories. If a list of directories is given "
                        "the program will find all CORRECT.LP files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("--match",
                        action="store", type=str, dest="match",
                        help="Include the log files only if they are in the given directories (i.e. --match SMV_reprocessed)")

    parser.add_argument("-g", "--gather",
                        action="store_true", dest="gather",
                        help="Gather refl files in local directory.")

    parser.add_argument("-j", "--job",
                        action="store", type=str, dest="job",
                        help="Type of ")

    parser.set_defaults(match=None, gather=False, job='index')

    options = parser.parse_args()

    match = options.match
    gather = options.gather
    args = options.args
    job = options.job

    if job == 'index':
        fns = parse_args_for_fns(args, name="dials.index.log", match=match)
        dials_all = []
        for fn in fns:
            try:
                p = dials_parser(fn, type='index')
            except UnboundLocalError as e:
                print(e)
                continue
            else:
                if p and p.d:
                    dials_all.append(p)
        for i, p in enumerate(dials_all):
            i += 1
            print(p.cell_info(sequence=i))


if __name__ == '__main__':
    main()
