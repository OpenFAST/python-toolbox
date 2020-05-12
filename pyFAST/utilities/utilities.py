"""Utility functions for the regression testing suite."""


import os
import sys
import subprocess
from stat import ST_MODE
from time import perf_counter


def ignore_baseline(_, contents):
    item_filter = ("linux-intel", "linux-gnu", "macos-gnu", "windows-intel")
    return [c for c in contents if c in item_filter]


def validate_directory(directory: str, create: bool = True):
    """
    Validates if a directory exists.

    Parameters
    ----------
    directory : str
        Path to the file of concern.
    create : bool
        Indicator for creating a the directory if it isn't found.

    Raises
    ------
    FileNotFoundError
        Error is raised if the file can't be found.
    """
    if not os.path.isdir(directory):
        if not create:
            raise FileNotFoundError(f"{directory} is not a valid directory")
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
    FileNotFoundError
        Error is raised if the file can't be found.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{file_path} is not a file")


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
    if int(permissions) % 2 != 1:
        raise PermissionError(f"{file_path} does not have proper permissions")


