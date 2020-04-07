
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
    Finally, a summary HTML file is created with plots, if requested.

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
            tolerance: float = 1e-5,
            plot: int = 0,
            plot_path: str = None,
            no_execution: bool = False,
            verbose: bool = False,
            jobs: bool = -1,
    ):
        """
        Initialize the required inputs

        NOTE: Make the plotting a little more modular so that all are done in one grid?

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
        tolerance: float, default: 1e-5
            Error tolerance for pass/fail condition.
        plot : int, default: 0
            Flag to include plotting:
             - 0 (default): No plots will be produced
             - 1: All plots will be produced.
             - 2: Only the plots for failing cases will be produced.
            All plots will be output to <path_to_case_name>/results.html.
        plot_path : str, default None
            Path to save all case result summaries and their plots. If `None`,
            the local output directory is used.
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
        self.local_test_location = os.path.join(self.build_directory, "local_results")

        system = platform.system() if system is None else system.lower()
        if system == "darwin":
            system = "macos"

        # Cases should be a list of case names ["awt_yfree", "..."]
        # to run all cases, pass an empty list or leave the argument out
        if cases == 0:
            self.cases = CASE_MAP
        else:
            self.cases = {case: CASE_MAP[case] for case in cases}

        self.output_type = "-".join((system, compiler.lower()))
        self.verbose = verbose
        self.no_execution = no_execution
        self.tolerance = tolerance
        self.plot = plot
        self.plot_path = plot_path
        self.jobs = jobs if jobs != 0 else -1

        for exe in executable:
            _exe = os.path.basename(os.path.normpath(exe))
            if "openfast" in _exe:
                self.of_executable = Path(exe)
            elif "beamdyn_driver" in _exe:
                self.bd_executable = Path(exe)

        self._validate_inputs()

        # Set the appropriate number of parallel jobs to run
        if self.jobs == -1:
            self.jobs = max(1, cpu_count() - 1)
        if self.jobs > 0:
            self.jobs = min(self.jobs, cpu_count())
        if self.jobs > len(self.cases):
            self.jobs = len(self.cases)

        # Initialize these variables for use when collecting test directories and opening files
        self.inputs = []                  # Directory in r-test containing the test case inputs
        self.baseline_directories = []    # Directory in r-test containing the test case baseline solutions
        self.local_case_directories = []  # Direcory containing the locally generated test case inputs and outputs

    def _validate_inputs(self):

        # Is the output type one of the supported combinations?
        _options = ("macos-gnu", "linux-intel", "linux-gnu", "windows-intel")
        if self.output_type not in _options:
            self.output_type = "macos-gnu"
            print(f"Defaulting to {self.output_type} for output type")

        # Are the required executable provided?
        # TODO
        # _required_executables = set([ CASE_MAP[c]["driver"] for c in self.cases ])

        # Do the given executables exist with the correct permissions?
        validate_executable(self.bd_executable)
        validate_executable(self.of_executable)

        # Do the given directories exist?
        validate_directory(self.build_directory)

        # Is the plot flag within the supported range?
        _options = (0, 1, 2)
        if self.plot not in _options:
            raise ValueError(f"Input 'plot' must be one of {_options}")

    def _build_beamdyn_output_directories(self, _to_build):
        """
        Creates the local output directories for BeamDyn cases and intializes
        it with the input files.
        """
        for case in _to_build:
            ix = self.case.index(case)
            in_dir, test = self.inputs[ix], self.test_build[ix]
            for bd_file in ("bd_driver.inp", "bd_primary.inp", "beam_props.inp"):
                shutil.copy(os.path.join(in_dir, bd_file), test)

    def _check_5MW_dll_files(self):
        """
        Checks for the .dll libraries in the 5MW Baseline folder and creates
        them if they don't exist.
        """

        source = os.path.join(self.module, "5MW_Baseline", "ServoData")
        target = os.path.join(self.build, "local_results", "5MW_Baseline", "ServoData")
        if not os.path.isdir(target):
            os.makedirs(target)

        discon = "DISCON/build/DISCON.dll"
        discon_itibarge = "DISCON_ITI/build/DISCON_ITIBarge.dll"
        discon_oc3hywind = "DISCON_OC3/build/DISCON_OC3Hywind.dll"
        for f in (discon, discon_itibarge, discon_oc3hywind):
            to_copy = os.path.join(source, f)
            _check = os.path.join(target, f.split("/")[-1])
            if not os.path.isfile(_check):
                shutil.copy2(to_copy, target)

    def _build_5MW_directories(self):
        """Copies the 5MW Baseline folder"""

        # NOTE: revisit this piece

        source = os.path.join(self.module, "5MW_Baseline")
        target = os.path.join(self.build, "local_results", "5MW_Baseline")
        if not os.path.isdir(target):
            shutil.copytree(source, target)
        else:
            for name in os.listdir(source):
                if name == "ServoData":
                    continue
                _source = os.path.join(source, name)
                _target = os.path.join(target, name)
                if os.path.isdir(_source):
                    if not os.path.isdir(_target):
                        shutil.copytree(_source, _target)
                else:
                    shutil.copy2(_source, _target)
        #  Is the jobs flag within the supported range?
        if self.jobs < -1:
            raise ValueError("Invalid value given for 'jobs'")

    def _build_local_test_directory(self):
        """
        Copies the input data to the test build directory
        """
        for source, destination in zip(self.inputs, self.local_case_directories):
            if not os.path.isdir(destination):
                shutil.copytree(source, destination, ignore=ignore_baseline)
            else:
                for f in os.listdir(source):
                    if os.path.isfile(f):
                        shutil.copy2(os.path.join(source, f), destination)

    def _build_test_output_directories(self):
        """
        Creates the local output directories
        """

        _linear = ("Ideal_Beam", "WP_Baseline")
        _regression = ("AOC", "AWT27", "SWRT", "UAE_VI", "WP_Baseline")
        directories = []
        _to_build_beamdyn = []
        _to_build_5mw = [c for c in self.case if "5MW" in c]

        case_types = set(CASE_MAP[c] for c in self.case)
        if "linear" in case_types:
            directories.extend(_linear)

        if "regression" in case_types:
            directories.extend(_regression)

        if "beamdyn" in case_types:
            _to_build_beamdyn = [c for c in self.case if CASE_MAP[c] == "beamdyn"]

        for data in directories:
            _dir = os.path.join(self.build, "local_results", data)
            if not os.path.isdir(_dir):
                shutil.copytree(os.path.join(self.module, data), _dir)

        self._build_beamdyn_output_directories(_to_build_beamdyn)
        self._check_5MW_dll_files()
        self._build_5MW_directories()
        self._build_local_test_directory()

    def _build_directory_references(self):
        """
        Creates lists of the directories that will be used
        throughout the test.
        
        These include
        - Directories containing the input files for a single test case
        - Directories containing the baseline results; for some cases this is the same as the
            input directory, but other cases have designated baseline directories
        - Directories that will be created locally to run the case on the tested system; these
            are copied from the "input" directories above into a "build" directory on the tested
            system. The local results are here.
        """
        for i, case in enumerate(self.cases):
            # The driver is currently either "openfast" or a specific module
            if CASE_MAP[case]["driver"] == "openfast":
                self.inputs.append(os.path.join(self.rtest_openfast, case))
                self.baseline_directories.append(os.path.join(self.inputs[i], self.output_type))
            else:
                self.inputs.append(os.path.join(self.rtest_modules, CASE_MAP[case]["driver"], case))
                self.baseline_directories.append(os.path.join(self.rtest_modules, CASE_MAP[case]["driver"], case))
            self.local_case_directories.append(os.path.join(self.build_directory, "reg_tests", "local_results", CASE_MAP[case]["driver"], case))

    def _run_openfast_case(
            self,
            executable: str,
            in_file: str,
            ix: str,
            case: str,
            verbose: bool = False,
        ):
        """
        Runs an OpenFAST regression test case.

        Parameters
        ----------
        executable : str
            File path to the OpenFAST executable.
        in_file : str
            Input file for the OpenFAST test case.
        ix : str
            String index/total of case being run.
        case : str
            Name of the case being run
        verbose : bool, optional
            Flag to include verbose output, by default False.
        """
        cwd = os.getcwd()
        os.chdir(os.path.dirname(in_file))

        stdout = sys.stdout if verbose else open(os.devnull, "w")

        validate_file(in_file)
        executable = os.path.abspath(executable)
        validate_executable(executable)

        base = os.path.sep.join(in_file.split(os.path.sep)[-1].split(".")[:-1])
        parent = os.path.sep.join(in_file.split(os.path.sep)[:-1])
        log = os.path.join(parent, "".join((base, ".log")))

        command = f"{executable} {in_file} > {log}"
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

    def _run_single_case(self, index: str, case: str, test_build: str):
        """
        Runs a single OpenFAST test case

        Parameters
        ----------
        ix : str
            String index as "i/n" for which case is being run.
        case : str
            Case name.
        test_build : str
                Testing build directory.
        """
        if CASE_MAP[case]["driver"] == "beamdyn":
            exe = self.bd_executable
            case_input = os.path.join(test_build, "bd_driver.inp")
        elif CASE_MAP[case]["driver"] == "openfast":
            exe = self.of_executable
            case_input = os.path.join(test_build, "".join((case, ".fst")))
        else:
            raise ValueError

        self._run_openfast_case(exe, case_input, index, case, verbose=self.verbose)

    def _run_openfast_cases(self):
        """
        Runs all of the OpenFAST cases in parallel, if defined.
        """
        indeces = [f"{i}/{len(self.cases)}" for i in range(1, len(self.cases) + 1)]
        arguments = list(zip(indeces, self.cases, self.local_case_directories))
        with Pool(self.jobs) as pool:
            pool.starmap(self._run_single_case, arguments)

    def run(self):
        """
        Function to build the references to ouput directories. If executing
        OpenFAST, then the directories are created if they don't already exist.
        """

        self._build_directory_references()
        if not self.no_execution:
            self._build_test_output_directories()
            self._run_openfast_cases()
    
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
        for directory, case in zip(self.baseline_directories, self.cases):
            results_file = os.path.join(directory, CASE_MAP[case]["reference_output"])
            try:
                validate_file(results_file)
            except FileNotFoundError as error:
                error_list.append( (case, "Baseline: " + str(error)) )
            else:
                data, info, _ = load_output(results_file)
                baseline_list.append( (data, info) )

        # Get the local test results
        for directory, case in zip(self.local_case_directories, self.cases):
            results_file = os.path.join(directory, CASE_MAP[case]["reference_output"])
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
