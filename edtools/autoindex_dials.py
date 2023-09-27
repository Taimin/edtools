import threading
import socket
import sys, os
import traceback
from pathlib import Path
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

rlock = threading.RLock()

def connect(payload: str) -> None:
    """Try to connect to `instamatic` indexing server

    Parameters
    ----------
    payload : str
        Directory where XDS should be run.
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


def parse_dials(path: str, sequence: int=0) -> None:
    """Parse XDS output (CORRECT.LP) and print summary about indexing progress
    to the screen.
    
    Parameters
    ----------
    path : str
        Path in which XDS has been run
    sequence : int
        Sequence number, needed for output and house-keeping
    """
    drc = Path(path)
    correct_lp = drc / "CORRECT.LP"

    lookBack = 160
    msg = None

    # rlock prevents messages getting mangled with 
    # simultaneous print statements from different threads
    with rlock:
        if correct_lp.exists():
        # if all files exist, try parsing CORRECT.LP
            try:
                p = dials_parser(correct_lp)
            except UnboundLocalError:
                msg = f"{sequence: 4d}: {drc} -> Indexing completed but no cell reported..."
            else:
                msg = "\n"
                msg += p.cell_info(sequence=sequence)
                msg += "\n"
                msg += p.info_header(hline=False)
                msg += p.integration_info(sequence=sequence)

            print(msg)
        else:
            for i, job in enumerate(XDSJOBS):
                error = None
                fn = (drc / job).with_suffix(".LP")
                if fn.exists():
                    with open(fn, "rb") as f:
                        f.seek(-lookBack, 2)
                        lines = f.readlines()

                    for line in lines:
                        if b"ERROR" in line:
                            error = line.decode()
                            error = error.split("!!!")[-1].strip()

                    if error:
                        msg = f"{sequence: 4d}: {drc} -> Error in {job}: {error}"
                        print(msg)
                        return
        with open(drc.parent.parent / "index_results.log", "a") as f:
            if msg is not None:
                print(msg, file=f)



def dials_process(path: str, sequence: int=0) -> None:
    """Run dials_process.bat at given path.
    
    Parameters
    ----------
    path : str
        Run dials_process.bat in this directory, expects dials_process.bat in this directory
    sequence : int
        Sequence number, needed for output and house-keeping
    """

    cmd = "dials_process.bat"

    cwd = str(path)

    try:
        p = sp.Popen(cmd, cwd=cwd, stdout=DEVNULL)
        p.wait()
    except Exception as e:
        print("ERROR in subprocess call:", e)

    try:
        parse_dials(path, sequence=sequence)
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

    parser.add_argument("-j", "--jobs",
                        action="store", type=int, dest="n_jobs",
                        help="Number of jobs to run in parallel")

    parser.set_defaults(use_server=False,
                        match=None,
                        unprocessed_only=False,
                        n_jobs=1,
                        )
    
    options = parser.parse_args()

    use_server = options.use_server
    match = options.match
    unprocessed_only = options.unprocessed_only
    n_jobs = options.n_jobs
    args = options.args

    fns = parse_args_for_fns(args, name="dials_process.bat", match=match)

    if unprocessed_only:
        fns = [fn for fn in fns if not fn.with_name("dials.index.log").exists()]
        print(f"Filtered directories which have already been processed, {len(fns)} left")


    drc = fns[0].parent
    with open(drc.parent.parent / "index_results.log", "w") as f:
        pass

    for i, fn in enumerate(fns):
        drc = fn.parent

        if use_server:
            connect(drc)
        else:
            dials_process(drc, i)
        futures.append(f)


if __name__ == '__main__':
    main()
