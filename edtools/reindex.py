import numpy as np
from pathlib import Path
import shutil
import subprocess as sp

from instamatic import config

from .utils import parse_args_for_fns
from .extract_dials_info import dials_parser
import os

CWD = Path(os.getcwd())
DEVNULL = open(os.devnull, 'w')

def main():
    import argparse

    description = "Program to reindex and resolve ambiguities."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of files or list of directories. If a list of directories is given "
                        "the program will find all dials_process.bat files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("-m", "--match",
                        action="store", type=str, dest="match",
                        help="Include files only if they are in the given directories (i.e. --match dir)")

    parser.add_argument("-i", "--input_file",
                        action="store", type=str, dest="input_file",
                        help="A input file that list all the directories")

    parser.add_argument("-p", "--pets2",
                        action="store", type=bool, dest="pets2",
                        help="Whether use pets2 file to reindex")

    parser.add_argument("-r", "--reindex",
                        action="store", type=bool, dest="reindex",
                        help="Whether use reindex file from DIALS to reindex")

    parser.add_argument("-sg", "--space_group",
                        action="store", type=str, dest="space_group",
                        help="Specify space group in reindex procedure")

    parser.set_defaults(pets2=False, reindex=False, space_group=None)

    options = parser.parse_args()
    fns = options.args
    match = options.match
    input_file = options.input_file
    pets2 = options.pets2
    reindex = options.reindex
    space_group =options.space_group

    if input_file:
        fns_index = []
        fns_reindex = []
        fns_pets2 = []
        with open(input_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.split('/')[-1].split('_')
                folder = ['_'.join(line[0:3]), '_'.join(line[3:5]), line[5]]
                folder = '/'.join(folder)
                file = Path('./' + folder) / "dials.index.log"
                fns_index.append(file)
                if pets2:
                    file = Path('./' + folder).parent / 'pets_petsdata' / "pets.celllist"
                    fns_pets2.append(file)
                if reindex:
                    file = Path('./' + folder) / "dials.reindex.log"
                    fns_reindex.append(file)
    else:
        fns_index = parse_args_for_fns(fns, name="dials.index.log", match=match)
        if pets2:
            fns_pets2 = parse_args_for_fns(fns, name="pets.celllist", match=match)
        if reindex:
            fns_reindex = parse_args_for_fns(fns, name="dials.reindex.log", match=match)

    if reindex:
        with open(CWD / 'reindex.log', 'w') as f:
            pass
    
    if reindex:
        with open(CWD / "reindex.log", "a") as fr:
            for fn in fns_reindex:
                with open(fn, 'r') as f:
                    lines = f.readlines()
                    fn = os.path.relpath(fn.parent, CWD)
                    for line in lines:
                        if line.startswith('No reindexing required'):
                            print(f'{fn} a, b, c', file=fr)
                        elif line.startswith('Reindexing required with the twin operator'):
                            operator = line.split(':')[-1].strip()
                            print(f'{fn} {operator}', file=fr)
    
    if pets2:
        transform_a_b = np.array([[0,1,0],[1,0,0],[0,0,1]])
        transform_a_c = np.array([[0,0,1],[0,1,0],[1,0,0]])
        transform_b_c = np.array([[1,0,0],[0,0,1],[0,1,0]])
        with open(CWD / "reindex.log", "a") as fr:
            for fn_p, fn_i in zip(fns_pets2, fns_index):
                with open(fn_p, 'r') as f:
                    lines = f.readlines()
                    for index, line in enumerate(lines):
                        if line.startswith('ubmatrix'):
                            print(f'{fn_p}', file=fr)
                            ubmatrix = []
                            ubmatrix.append([float(num) for num in lines[index+1].split()])
                            ubmatrix.append([float(num) for num in lines[index+2].split()])
                            ubmatrix.append([float(num) for num in lines[index+3].split()])
                            ubmatrix = np.array(ubmatrix)
                            print(f'{ubmatrix}', file=fr)
                        elif line.startswith('cell  '):
                            uc = np.array([float(num) for num in line.split()[1:]])
                            print(f'{uc}', file=fr)
                with open(fn_i, 'r') as f:
                    crystals = dials_parser(fn_i, job='index')
                    fn_i = os.path.relpath(fn_i.parent, CWD)
                    for crystal in crystals.d:
                        print(f'{fn_i}', file=fr)
                        A_matrix = crystal['A_matrix']
                        print(f'{A_matrix}', file=fr)
                        uc = np.array(crystal['cell'])
                        print(f'{uc}', file=fr)
                    
                ubmatrix = np.abs(ubmatrix)
                A_matrix = np.abs(A_matrix)
                diff_a_c =  np.abs(A_matrix@transform_a_c - ubmatrix).sum()
                diff_a_c_b =  np.abs(A_matrix@transform_a_c@transform_a_b - ubmatrix).sum()
                diff_a_b =  np.abs(A_matrix@transform_a_b - ubmatrix).sum()
                diff = np.abs(A_matrix - ubmatrix).sum()
                if space_group is None:
                    cmd = f'dials.reindex.bat indexed.expt indexed.refl change_of_basis_op='
                else:
                    cmd = f'dials.reindex.bat indexed.expt indexed.refl space_group={space_group} change_of_basis_op='
                
                try:
                    #print(f'a<->c<->b: {diff_a_c_b}, a<->c: {diff_a_c}', file=fr)
                    print(f'a<->b: {diff_a_b}, same: {diff}', file=fr)
                    if diff_a_b < diff:
                        p = sp.Popen(cmd+'b,a,-c', cwd=str(CWD/fn_i), stdout=DEVNULL)
                    else:
                        p = sp.Popen(cmd+'a,b,c', cwd=str(CWD/fn_i), stdout=DEVNULL)
                    p.communicate()
                except Exception as e:
                    print("ERROR in subprocess call:", e)
                print('', file=fr)

