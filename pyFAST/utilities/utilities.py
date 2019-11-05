"""Utility functions for the regression testing suite."""


import os
import subprocess
import sys
from stat import ST_MODE


def ignore_baseline(directory, contents):
    itemFilter = ("linux-intel", "linux-gnu", "macos-gnu", "windows-intel")
    caught = []
    for c in contents:
        if c in itemFilter:
            caught.append(c)
    return caught


def validate_directory(directory: str, create: bool = True):
    """
    Validates if a directory exists exists.

    Parameters
    ----------
    directory : str
        Path to the file of concern.
    create : bool
        Indicator for creating a the directory if it isn't found.

    Raises
    ------
    FileExistsError
        Error is raised if the file can't be found.
    """
    if not os.path.isdir(directory):
        if not create:
            raise FileExistsError(f"{directory} is not a valid directory")
        os.makedirs(directory)


def validate_file(file_path: str):
    """
    Validates if a file exists.

    Parameters
    ----------
    file_path : str
        Path to the file of concern.

    Raises
    ------
    FileExistsError
        Error is raised if the file can't be found.
    """
    if not os.path.isfile(file_path):
        raise FileExistsError(f"{file_path} is not a valid file")


def validate_executable(file_path: str):
    """
    Validates if the passed file is a valid executable.

    Parameters
    ----------
    file_path : str
        Path to the exectuable file.

    Raises
    ------
    PermissionError
        Error is raised if the file has incorrect permissions.
    """
    permissions = oct(os.stat(file_path)[ST_MODE])[-1:]
    if permissions % 2 != 1:
        raise PermissionError(f"{file_path} does not have proper permissions")


def run_openfast_case(
    executable: str, in_file: str, verbose: bool = False, beamdyn: str = False
):
    """
    Runs an OpenFAST regression test case.

    Parameters
    ----------
    in_file : str
        Input file for the OpenFAST test case.
    executable : str
        File path to the OpenFAST executable.
    verbose : bool, optional
        Flag to include verbose output, by default False.
    """

    if beamdyn:
        case_dir = os.path.sep.join(in_file.split(os.path.sep)[:-1])
        os.chdir(case_dir)

    stdout = sys.stdout if verbose else open(os.devnull, "w")

    validate_file(in_file)
    validate_executable(executable)

    base = os.path.sep.join(in_file.split(os.path.sep)[-1].split(".")[:-1])
    parent = os.path.sep.join(in_file.split(os.path.sep)[:-1])
    log = os.path.join(parent, "".join((base, ".log")))

    command = f"{executable} {in_file} > {log}"
    print(command)
    code = subprocess.call(command, stdout=stdout, shell=True)
    print(f"COMPLTE with code: {code}", flush=True)
    return code
