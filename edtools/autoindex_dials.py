import threading
import socket
import shutil
import sys, os
import traceback
from pathlib import Path
import yaml
import concurrent.futures
from .utils import parse_args_for_fns

import subprocess as sp

from .extract_dials_info import dials_parser

try:
    from instamatic import config
    HOST = config.settings.indexing_server_host
    PORT = config.settings.indexing_server_port
    BUFF = 1024
except ImportError:
    HOST, PORT = None, None

DEVNULL = open(os.devnull, 'w')
CWD = Path(os.getcwd())
rlock = threading.RLock()

def connect(payload: str) -> None:
    """Try to connect to `instamatic` indexing server

    Parameters
    ----------
    payload : str
        Directory where DIALS should be run.
    """
    payload = str(payload).encode()

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        with rlock:
            print("Sending job to server...", end=" ")
            s.connect((HOST, PORT))
            s.send(payload)
            data = s.recv(BUFF).decode()
            print(data)

        data = s.recv(BUFF).decode()

        with rlock:
            print(data)


def parse_dials_index(path: str, sequence: int=0) -> None:
    """Parse DIALS index output and print summary about indexing progress
    to the screen.
    
    Parameters
    ----------
    path : str
        Path in which dials.index has been run
    sequence : int
        Sequence number, needed for output and house-keeping
    """
    drc = Path(path)
    dials_index_log = drc / "dials.index.log"
    msg = None

    # rlock prevents messages getting mangled with 
    # simultaneous print statements from different threads
    with rlock:
        if dials_index_log.exists():
            try:
                p = dials_parser(dials_index_log, job='index')
            except UnboundLocalError:
                msg = f"{sequence: 4d}: {drc} -> Indexing completed but no cell reported..."
            else:
                msg = "\n"
                msg += "".join(p.cell_info(sequence=sequence))
                msg += "\n"
            print(msg)

        with open(CWD / "index_results.log", "a") as f:
            if msg is not None:
                print(msg, file=f)

def process_data(i, fn, job, restrain, find_spot, use_server, only_index):
    drc = fn.parent
    if use_server:
        connect(drc)
    else:
        cwd = str(drc)
        if restrain:
            shutil.copy(str(CWD/'restrain.phil'), str(drc/'restrain.phil'))
        else:
            with open(drc/'restrain.phil', 'w') as f:
                print('indexing {', file=f)
                print('}', file=f)
        if find_spot:
            shutil.copy(str(CWD/'find_spot.phil'), str(drc/'find_spot.phil'))
        else:
            with open(drc/'find_spot.phil', 'w') as f:
                print('spotfinder {', file=f)
                print('}', file=f)
        cmd = str(drc/"dials_process.bat")
        if only_index:
            cmd = 'dials.index.bat imported.expt strong.refl restrain.phil nproc=2'
        if os.path.exists(drc/'indexed.expt'): 
            os.remove(drc/'indexed.expt')
        if os.path.exists(drc/'indexed.refl'): 
            os.remove(drc/'indexed.refl')
        try:
            p = sp.Popen(cmd, cwd=cwd, stdout=DEVNULL)
            p.communicate()
        except Exception as e:
            print("ERROR in subprocess call:", e)
        try:
            if job == 'index':
                parse_dials_index(drc, sequence=i)
        except Exception:
            traceback.print_exc()

def main():
    import argparse

    description = "Program to automate the indexing of a large series of data sets from a serial crystallography experiment."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
        
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of dials_process.bat files or list of directories. If a list of directories is given "
                        "the program will find all dials_process.bat files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("-s","--server",
                        action="store_true", dest="use_server",
                        help="Use instamatic server for indexing")

    parser.add_argument("-m", "--match",
                        action="store", type=str, dest="match",
                        help="Include the dials_process.bat files only if they are in the given directories (i.e. --match SMV_reprocessed)")

    parser.add_argument("-u", "--unprocessed_only",
                        action="store_true", dest="unprocessed_only",
                        help="Run dials_process.bat only in unprocessed directories ")

    parser.add_argument("-j", "--job",
                        action="store", type=str, dest="job",
                        help="Type of ")

    parser.add_argument("-r", "--restrain",
                        action="store", type=bool, dest="restrain",
                        help="If True, copy the restrain.phil file for indexing to each dials directory")

    parser.add_argument("-f", "--find_spot",
                        action="store", type=bool, dest="find_spot",
                        help="If True, copy the find_spot.phil file to each dials directory")

    parser.add_argument("-min", "--min_spots",
                        action="store", type=int, dest="min_spots",
                        help="If True, copy the find_spot.phil file to each dials directory")

    parser.add_argument("-oi", "--only_index",
                        action="store", type=bool, dest="only_index",
                        help="Only index")

    parser.set_defaults(use_server=False,
                        match=None,
                        unprocessed_only=False,
                        job='index',
                        restrain=False,
                        find_spot=False,
                        min_spots=None,
                        only_index=False
                        )

    options = parser.parse_args()

    use_server = options.use_server
    match = options.match
    unprocessed_only = options.unprocessed_only
    job = options.job
    args = options.args
    restrain = options.restrain
    find_spot = options.find_spot
    min_spots = options.min_spots
    only_index = options.only_index

    if args:
        fns = []
        for arg in args:
            if arg.split('.')[-1] == "yaml":
                ds = yaml.load(open(arg, "r"), Loader=yaml.Loader)
                for d in ds:
                    fns.append(Path(d['directory']) / "dials_process.bat")
            elif arg.split('.')[-1] == "lst":
                with open(arg, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.split('/')[-1].split('_')
                        folder = ['_'.join(line[0:3]), '_'.join(line[3:5])]
                        folder = '/'.join(folder)
                        file = Path('./' + folder) / "dials_process.bat"
                        fns.append(file)
        fns = [fn.resolve() for fn in fns]
    else:
        fns = parse_args_for_fns(args, name="dials_process.bat", match=match)

    if unprocessed_only:
        fns = [fn for fn in fns if not fn.with_name("dials.index.log").exists()]
        print(f"Filtered directories which have already been processed, {len(fns)} left")

    if len(fns) == 0:
        return -1

    if (CWD / "index_results.log").is_file():
        shutil.copy(str(CWD/'index_results.log'), str(CWD/'~index_results.log'))
    with open(CWD / "index_results.log", "w") as f:
        pass

    if min_spots is not None:
        with open(CWD / 'find_spot.phil', 'w') as f:
            print('spotfinder {', file=f)
            print('  filter {', file=f)
            print(f'    min_spot_size = {min_spots}', file=f)
            print('    max_spot_size = 100', file=f)
            print('  }', file=f)
            print('  threshold {', file=f)
            print('    algorithm = *dispersion dispersion_extended radial_profile', file=f)
            print('    dispersion {', file=f)
            print('      gain = 1', file=f)
            print('      sigma_strong = 3', file=f)
            print('      global_threshold = 1', file=f)
            print('    }', file=f)
            print('  }', file=f)
            print('}', file=f)

    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        for i, fn in enumerate(fns):
            futures.append(executor.submit(process_data, i, fn, job, restrain, find_spot, use_server, only_index))
    concurrent.futures.wait(futures, return_when=concurrent.futures.ALL_COMPLETED)
        

if __name__ == '__main__':
    main()
