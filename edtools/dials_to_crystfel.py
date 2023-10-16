import numpy as np
from pathlib import Path
import shutil

from instamatic import config
from instamatic.formats import read_image
from instamatic.formats import write_adsc
from instamatic.formats import write_tiff

from .utils import parse_args_for_fns


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

    for fn in fns:

        print("\033[K", fn, end='\r')  # "\033[K" clears line


    print(f"\033[KUpdated {len(fns)} files")