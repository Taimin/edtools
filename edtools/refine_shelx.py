import numpy as np
from pathlib import Path
import shutil
import traceback


from .utils import parse_args_for_fns

def parse_shelx():
	pass


def main():
    import argparse

    description = "Program to refine structures using shelxl."
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
