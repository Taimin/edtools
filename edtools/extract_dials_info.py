from pathlib import Path
import os
from datetime import datetime
import shutil
import re
from .utils import volume, parse_args_for_fns
from .utils import space_group_lib
import copy
import pandas as pd
import numpy as np
import yaml

spglib = space_group_lib()
CWD = Path(os.getcwd())

class dials_parser:
    """docstring for dials_parser"""
    def __init__(self, filename, job='index'):
        self.ios_threshold = 0.8

        self.filename = Path(filename).resolve()
        self.d = self.parse(job)

    def parse(self, job):
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
                    d['percent'] = int(infos[2]) / (int(infos[2]) + int(infos[3]))
                elif line.startswith("    A = UB:"):
                    A_matrix = []
                    tmp = line.split('{')[-1].split(',')[:3]
                    tmp[2] = tmp[2].split('}')[0]
                    A_matrix.append([float(num) for num in tmp])
                    tmp = crystal[index+1].strip().split('{')[-1].split(',')[:3]
                    tmp[2] = tmp[2].split('}')[0]
                    A_matrix.append([float(num) for num in tmp])
                    tmp = crystal[index+2].strip().split('{')[-1].split(',')[:3]
                    tmp[2] = tmp[2].split('}')[0]
                    A_matrix.append([float(num) for num in tmp])
                    A_matrix = np.array(A_matrix).flatten()
            d["volume"] = volume(cell)
            d["cell"] = cell
            d["spgr"] = spgr
            d["model_num"] = model_num
            d['A_matrix'] = A_matrix
            d["fn"] = str(os.path.relpath(fn.parent, CWD))
            d_list.append(copy.deepcopy(d))
        
        return d_list

    def print_filename(self):
        print("#", self.filename)

    def cell_info(self, sequence=0):
        d = self.d
        i = sequence
        s_list = []
        for crystal in d:
            fn = self.filename
            model_num = crystal["model_num"]
            s = f"{i: 4d}: {fn.parents[0]} # Time: {datetime.fromtimestamp((os.path.getmtime(fn))).strftime('%Y-%m-%d %H:%M:%S')} #\n"
            s += "Spgr {} - Cell {:10.2f}{:10.2f}{:10.2f}{:10.2f}{:10.2f}{:10.2f} - Vol {:10.2f} Indexed: {} Unindexed: {}\n"\
                    .format(crystal["spgr"], *crystal["cell"], crystal["volume"], crystal["indexed"], crystal["unindexed"])
            s += f"# Model {model_num} # Indexed: {crystal['indexed']} # Unindexed: {crystal['unindexed']} \
                    # Percentage: {crystal['percent']:.1%} #\n"
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

def cells_to_yaml(ps, fn="cells.yaml"):
    import yaml
    ds = []
    for i, p in enumerate(ps):
        i += 1
        d = {}
        for crystal in p.d:
            d["directory"] = crystal["fn"]
            d["number"] = i
            d["unit_cell"] = crystal["cell"]
            d["space_group"] = crystal["spgr"]
            d["indexed"] = crystal["indexed"]
            d["percent"] = crystal["percent"]
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
                        help="Type of job performed")

    parser.add_argument("-sc", "--single_crystal",
                        action="store", type=bool, dest="single_crystal",
                        help="Only obtain datasets with only one lattice.")

    parser.add_argument("-t_i", "--thresh_indexed",
                        action="store", type=float, dest="thresh_indexed",
                        help="Threshold based on indexed number of spots.")

    parser.add_argument("-t_p", "--thresh_percent",
                        action="store", type=float, dest="thresh_percent",
                        help="Threshold based on percentage of indexed spots.")

    parser.set_defaults(match=None, gather=False, job='index', single_crystal=True, thresh_percent=None, thresh_indexed=None)

    options = parser.parse_args()

    match = options.match
    gather = options.gather
    args = options.args
    job = options.job
    single_crystal = options.single_crystal
    thresh_percent = options.thresh_percent
    thresh_indexed = options.thresh_indexed

    if job == 'index':
        if args:
            fns = []
            for arg in args:
                ds = yaml.load(open(arg, "r"), Loader=yaml.Loader)
                for d in ds:
                    fns.append(Path(d['directory']) / "dials.index.log")
            fns = [fn.resolve() for fn in fns]
        else:
            fns = parse_args_for_fns(args, name="dials.index.log", match=match)
        dials_all = []
        records = []
        for fn in fns:
            try:
                p = dials_parser(fn, job=job)
            except UnboundLocalError as e:
                print(e)
                continue
            else:
                if p and p.d:
                    if single_crystal:
                        if len(p.d) > 1:
                            continue
                    if thresh_indexed is not None:
                        if p.d[0]['indexed'] < thresh_indexed:
                            continue
                    if thresh_percent is not None:
                        if p.d[0]['percent'] < thresh_percent:
                            continue
                    dials_all.append(p)
                    for crystal in p.d:
                        records.append(crystal)

        for i, p in enumerate(dials_all):
            i += 1
            for info in p.cell_info(sequence=i):
                print(info)

        df = pd.DataFrame.from_dict(records)
        df.to_csv(CWD/'unit_cell.csv')

        # cells_to_cellparm(xdsall)
        cells_to_yaml(dials_all)

if __name__ == '__main__':
    main()
