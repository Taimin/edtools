import numpy as np
from pathlib import Path
import shutil
import traceback
import pandas as pd


from .utils import parse_args_for_fns

def parse_shelx(fn, input_name):
	# Obtain data statistics after refinement
	res_file = f'{input_name}.res'
	lst_file = f'{input_name}.lst'

	with open(res_file, 'r') as f:
		pass

	with open(lst_file, 'r') as f:
		pass

def run_shelx(fn, input_name):
	# run shelxl on every dataset
	cmd = f"shelxl {input_name}.ins"
	try:
        p = sp.Popen(cmd, cwd=cwd, stdout=DEVNULL)
        p.communicate()
    except Exception as e:
        print("ERROR in subprocess call:", e)

def read_conf(fn):
	pass

def edit_shelx_input(fn, input_name, sc_factor=None):
	# Obtain data statistics after refinement
	pass


def main():
    import argparse

    description = "Program to refine structures using shelxl."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="Name of the res files.")

    parser.add_argument("-m", "--match",
                        action="store", type=str, dest="match",
                        help="Include the dials_process.bat files only if they are in the given directories (i.e. --match DIALS_reprocessed)")

    parser.add_argument("-s", "--save",
                        action="store", type=str, dest="save",
                        help="A directory to save the refined results")

    parser.add_argument("-c", "--conf",
                        action="store", type=str, dest="conf",
                        help="Specify a configuration file. Needed when extra jobs are specified")

    parser.add_argument("-chg", "--charge",
                        action="store", type=bool, dest="charge",
                        help="Use the scattering parameter specified in conf file to investigate the charge effect.")

    parser.add_argument("-r", "--refine",
                        action="store", type=bool, dest="refine",
                        help="Start refinement using shelxl.")
    
    parser.set_defaults()

    options = parser.parse_args()
    args = options.args
    match = options.match
    save = options.save
    conf = options.conf
    charge = options.charge
    refine = options.refine

    for arg in args:
    	fns = parse_args_for_fns(None, name=arg, match=match)
    
    if refine:
	    if charge:
	    	pass
	    	# read conf file

	    	# make combination for the scattering factor

	    	# change res file

	    	# run shelxl for each changed res file

	    	# add the refinement result to the pandas dataframe

	    else:
	    	# run shelxl for each res file

	    	# add the refinement result to the pandas dataframe