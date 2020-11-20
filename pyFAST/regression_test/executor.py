
import os
import sys
import shutil
import platform
from typing import List, Tuple
from pathlib import Path
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
from time import perf_counter
import subprocess

from .utilities import (
    validate_file,
    validate_directory,
    validate_executable,
    ignore_baseline,
    load_output,
)
from .case_map import CASE_MAP


class Executor:
    """
    Base execution class for OpenFAST regression tests.

    This class first constructs the internal directories containing the test case input files.
    Then, the test cases are executed using multiple processors, if requested.
    The locally generated outputs are tested against their corresponding baselines.

    Attributes
    ----------
    """

    def __init__(
            self,
            cases: List[str],
            executable: List[str],
            openfast_root: str,
            compiler: str,
            system: str = None,
            no_execution: bool = False,
            verbose: bool = False,
            jobs: bool = -1,
    ):
        """
        Initialize the required inputs

        TODO:
        - There should exist a test pipeline that executes the full regression test process for each case
        - Thats what should be parallelized, not just the case execution
        - As is, the regression test results must wait for all tests to finish executing
        - We want to be able to bail if one test case fails the regression test but others havent finished

        - Choose to use str or Path for all paths

        Parameters
        ----------
        case : List[str]
            Test case name(s) as a list of strings.
        executable : List[str]
            Path(s) to the OpenFAST executable(s).
        openfast_root : str
            Path to OpenFAST repository.
        compiler : str
            System compiler id. Should be one of "intel" or "gnu".
        system : str
            Operating system version of results used for comparison, default
            None (machine's OS). Should be one of "windows", "linux", or
            "macos".
        no_execution : bool, default: False
            Flag to avoid executing the simulations, but proceed with the regression test.
        verbose : bool, default: False
            Flag to include system ouptut.
        jobs : int, default: -1
            Maximum number of parallel jobs to run:
             - -1: Number of nodes available minus 1
             - >0: Minimum of the number passed and the number of nodes available
        """

        # These path variables are used throughout Executor for brevity
        self.build_directory = os.path.join(openfast_root, "build")
        self.rtest_modules = os.path.join(openfast_root, "reg_tests", "r-test", "modules")
        self.rtest_openfast = os.path.join(openfast_root, "reg_tests", "r-test", "glue-codes", "openfast")
        self.local_test_location = os.path.join(self.build_directory, "reg_tests", "local_results")

        system = platform.system() if system is None else system.lower()
        if system == "darwin":
            system = "macos"

        # Cases should be a list of case names ["awt_yfree", "..."]
        # to run all cases, pass an empty list
        self.cases = CASE_MAP if not cases else {case: CASE_MAP[case] for case in cases}

        self.output_type = "-".join((system, compiler.lower()))
        self.verbose = verbose
        self.no_execution = no_execution
        self.jobs = jobs if jobs != 0 else -1

        for exe in executable:
            _exe = os.path.basename(os.path.normpath(exe))
            if "openfast" in _exe:
                self.of_executable = os.path.abspath(exe)
            elif "beamdyn_driver" in _exe:
                self.bd_executable = os.path.abspath(exe)

        self._validate_inputs()

        # Set the appropriate number of parallel jobs to run
        if self.jobs == -1:
            self.jobs = max(1, cpu_count() - 1)
        if self.jobs > 0:
            self.jobs = min(self.jobs, cpu_count())
        if self.jobs > len(self.cases):
            self.jobs = len(self.cases)

    def _validate_inputs(self):

        # Is the output type one of the supported combinations?
        _options = ("macos-gnu", "linux-intel", "linux-gnu", "windows-intel")
        if self.output_type not in _options:
            self.output_type = "macos-gnu"
            print(f"Defaulting to {self.output_type} for output type")

        # Are the required executables provided in a valid location with correct permissions?
        _required_executables = set([ CASE_MAP[c]["driver"] for c in self.cases ])
        if "openfast" in _required_executables:
            try:
                assert self.of_executable
            except AttributeError as error:
                raise AttributeError("An OpenFAST case was requested but no OpenFAST executable given.")
            validate_executable(self.of_executable)

        if "beamdyn" in _required_executables:
            try:
                assert self.bd_executable
            except AttributeError as error:
                raise AttributeError("A BeamDyn case was requested but no BeamDyn Driver executable given.")
            validate_executable(self.bd_executable)

        # Do the given directories exist?
        validate_directory(self.build_directory)

        #  Is the jobs flag within the supported range?
        if self.jobs < -1:
            raise ValueError("Invalid value given for 'jobs'")

    def _build_local_case_directories(self):
        """
        Copies the input data to the local directories where the tests will be run
        """
        for case in self.cases:
            case_info = CASE_MAP[case]

            if CASE_MAP[case]["driver"] == "openfast":
                # Copy the case files
                destination = os.path.join(self.local_test_location, "glue-codes",  case_info["driver"], case)
                if not os.path.isdir(destination):
                    source = os.path.join(self.rtest_openfast, case)
                    shutil.copytree(source, destination, ignore=ignore_baseline)
                
                # Copy the required turbine model
                destination = os.path.join(self.local_test_location, "glue-codes", case_info["driver"], case_info["turbine_directory"])
                if not os.path.isdir(destination):
                    source = os.path.join(self.rtest_openfast, case_info["turbine_directory"])
                    shutil.copytree(source, destination, ignore=ignore_baseline)

            else:
                # Copy the case files
                destination = os.path.join(self.local_test_location, "modules", case_info["driver"], case)
                if not os.path.isdir(destination):
                    source = os.path.join(self.rtest_modules, case_info["driver"], case)
                    shutil.copytree(source, destination)

    def _execute_case(
            self,
            executable: str,
            case_directory: str,
            input_file: str,
            ix: str,
            case: str,
            verbose: bool = False,
        ):
        """
        Runs an OpenFAST regression test case.

        Parameters
        ----------
        executable : str
            File path to the approptiate executable for this case.
        case_directory : str
            Directory containing the test case files.
        input_file : str
            Input file for the test case.
        ix : str
            String index/total of case being run.
        case : str
            Name of the case being run
        verbose : bool, optional
            Flag to include verbose output, by default False.
        """
        cwd = os.getcwd()
        os.chdir(case_directory)

        stdout = sys.stdout if verbose else open(os.devnull, "w")

        validate_file(input_file)
        executable = os.path.abspath(executable)
        validate_executable(executable)

        base = os.path.sep.join(input_file.split(os.path.sep)[-1].split(".")[:-1])
        parent = os.path.sep.join(input_file.split(os.path.sep)[:-1])
        log = os.path.join(parent, "".join((base, ".log")))

        command = f"{executable} {input_file} > {log}"
        print(f"{ix.rjust(6)} Start: {case}")
        if verbose:
            print(f"command: {command}")

        start = perf_counter()
        code = subprocess.call(command, stdout=stdout, shell=True)
        end = perf_counter()
        elapsed = f"{end - start:.2f}"
        status = "FAILED".rjust(8) if code != 0 else "complete"
        ix = ix.split("/")[0]
        message = (
            f"{ix.rjust(6)}   End: {case.ljust(40, '.')} {status} with code "
            f"{code}{elapsed.rjust(8)} seconds"
        )
        print(message, flush=True)

        os.chdir(cwd)

    def _run_case(self, index: str, case: str):
        """
        Runs a single OpenFAST test case

        Parameters
        ----------
        index : str
            String index as "i/n" for which case is being run.
        case : str
            Case name.
        """

        case_info = CASE_MAP[case]
        if CASE_MAP[case]["driver"] == "openfast":
            case_directory = os.path.join(self.local_test_location, "glue-codes", case_info["driver"], case)
        else:
            case_directory = os.path.join(self.local_test_location, "modules", case_info["driver"], case)

        if CASE_MAP[case]["driver"] == "beamdyn":
            exe = self.bd_executable
            case_input = "bd_driver.inp"
        elif CASE_MAP[case]["driver"] == "openfast":
            exe = self.of_executable
            case_input = "".join((case, ".fst"))
        else:
            raise ValueError

        self._execute_case(exe, case_directory, case_input, index, case, verbose=self.verbose)

    def _run_cases(self):
        """
        Runs all of the OpenFAST cases in parallel, if defined.
        """
        indeces = [f"{i}/{len(self.cases)}" for i in range(1, len(self.cases) + 1)]
        arguments = list(zip(indeces, self.cases))
        with Pool(self.jobs) as pool:
            pool.starmap(self._run_case, arguments)

    def run(self):
        """
        Function to build the references to ouput directories. If executing
        OpenFAST, then the directories are created if they don't already exist.
        """

        if not self.no_execution:
            self._build_local_case_directories()
            self._run_cases()
    
    def read_output_files(self) -> Tuple[list]:
        """
        Reads in the output files corresponding to `case` and returns the
        cases, baseline outputs, and locally produced outputs.

        Returns
        -------
        case_list : list
            List of valid non-linear regression test cases.
        baseline_list : List[tuple]
            List of valid non-linear regression test cases (attribute, data).
        test_list : List[tuple]
            List of valid non-linear regression test cases (attribute, data).
        """

        test_list = []
        baseline_list = []
        error_list = []

        # Get the baseline results
        for case in self.cases:

            case_info = CASE_MAP[case]
            if CASE_MAP[case]["driver"] == "openfast":
                results_file = os.path.join(self.rtest_openfast, case, self.output_type, CASE_MAP[case]["reference_output"])
            else: 
                results_file = os.path.join(self.rtest_modules, case_info["driver"], case, CASE_MAP[case]["reference_output"])

            try:
                validate_file(results_file)
            except FileNotFoundError as error:
                error_list.append( (case, "Baseline: " + str(error)) )
            else:
                data, info, _ = load_output(results_file)
                baseline_list.append( (data, info) )

        # Get the local test results
        for case in self.cases:

            case_info = CASE_MAP[case]
            if CASE_MAP[case]["driver"] == "openfast":
                results_file = os.path.join(self.local_test_location, "glue-codes", case_info["driver"], case, CASE_MAP[case]["reference_output"])
            else:
                results_file = os.path.join(self.local_test_location, "modules", case_info["driver"], case, CASE_MAP[case]["reference_output"])

            try:
                validate_file(results_file)
            except FileNotFoundError as error:
                error_list.append( (case, "Local: " + str(error)) )
            else:
                data, info, _ = load_output(results_file)
                test_list.append( (data, info) )

        if error_list:
            print("\n\x1b[1;31mErrors:\x1b[0m")
            for case, e in error_list:
                case = f"\x1b[1;31m{case}\x1b[0m"
                print(f"  {case}: {e}")
            print("")

        return baseline_list, test_list
