import numpy as np
import scipy.ndimage as ndimage
from pathlib import Path
import shutil

from instamatic import config
from instamatic.formats import read_image, write_tiff
from instamatic.processing.PETS2_template import PETS2_template

from .utils import parse_args_for_fns
import traceback


def main():
    import argparse

    description = "Program to convert tiff file and update pets processing commands."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of files or list of directories. If a list of directories is given "
                        "the program will find all dials_process.bat files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("-m", "--match",
                        action="store", type=str, dest="match",
                        help="Include the files only if they are in the given directories (i.e. --match file)")

    parser.add_argument("-w_t", "--write_tif",
                        action="store", type=bool, dest="write_tif",
                        help="Write Tiffs")
