"""Provides the Executor case"""


import os
import shutil
import platform
from typing import List, Tuple
from pathlib import Path
from functools import partial
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
import numpy as np

from .utilities import (
    load_output,
    validate_file,
    calculate_norms,
    ignore_baseline,
    validate_directory,
    validate_executable,
    pass_regression_test
)
from .case_map import CASE_MAP



class Executor:
    """
    Base execution class for OpenFAST regression tests. 

    Attributes
    ----------
    """

    def __init__(
            self,
            case: List[str],
            executable: List[str],
            openfast_root: str,
            compiler: str,
            system: str = None,
            tolerance: float = 1e-5,
            plot: int = 0,
            plot_path: str = None,
            execution: bool = False,
            verbose: bool = False,
            jobs: bool = -1,
    ):
        """
        Initialize the required inputs

        NOTE: Make the plotting a little more modular so that all are done in one grid?

        Parameters
        ----------
        case : List[str]
            Test case name(s) as a list of strings.
        executable : List[str]
            Path(s) to the OpenFAST executable(s). Should be no more than
            length 2 with one exe being for OpenFAST and the other for beamdyn.
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
        execution : bool, default: True
            Flag to run the test case(s). If `False`, ....
        verbose : bool, default: False
            Flag to include system ouptut.
        jobs : int, default: -1
            Maximum number of parallel jobs to run:
             - -1: 80% of maximum number of nodes available
             - >0: Minimum of the number passed and the number of nodes available
        """

        system = platform.system() if system is None else system.lower()
        if system == "darwin":
            system = "macos"

        if case == "all":
            self.case = [*CASE_MAP]
        else:
            self.case = case if isinstance(case, list) else [case]
        self.compiler = compiler
        self.output_type = "-".join((system, self.compiler.lower()))

        self.root = Path(openfast_root)
        self.build = os.path.join(self.root, "build")

        self.verbose = verbose
        self.execution = execution
        self.tolerance = tolerance
        self.plot = plot
        self.plot_path = plot_path
        self.jobs = jobs if jobs != 0 else -1

        self.rtest = os.path.join(self.root, "reg_tests", "r-test")
        self.module = os.path.join(self.rtest, "glue-codes", "openfast")

        for exe in executable:
            _exe = os.path.basename(os.path.normpath(exe))
            if "openfast" in _exe:
                self.of_executable = Path(exe)
            elif "beamdyn_driver" in _exe:
                self.bd_executable = Path(exe)

        self._validate_inputs()

        # Initialize these variables for use when collecting test directories
        # and opening files.
        self.inputs = []                  # Directory in r-test containing the test case inputs
        self.baseline_directories = []    # Directory in r-test containing the test case baseline solutions
        self.local_case_directories = []  # Direcory containing the locally generated test case inputs and outputs
        self.local_test_location = os.path.join(self.build, "local_results")

    def _validate_inputs(self):
        """Method to ensure inputs are valid."""

        _opts = ("macos-gnu", "linux-intel", "linux-gnu", "windows-intel")
        if self.output_type not in _opts:
            self.output_type = "macos-gnu"
            print(f"Defaulting to {self.output_type} for output type")

        if self.bd_executable != self.root:
            validate_executable(self.bd_executable)
        if self.of_executable != self.root:
            validate_executable(self.of_executable)

        validate_directory(self.build)

        _opts = (0, 1, 2)
        if self.plot not in _opts:
            raise ValueError(f"Input 'plot' must be one of {_opts}!")

        if self.jobs < -1:
            raise ValueError("Input 'jobs' cannot be negative!")
        if self.jobs == -1:
            self.jobs = int(np.ceil(cpu_count() * 0.8))
        if self.jobs > 0:
            self.jobs = min(self.jobs, cpu_count())
        if self.jobs > len(self.case):
            self.jobs = len(self.case)

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

    def _build_local_test_directory(self):
        """
        Copies the input data to the test build directory
        """
        for input_dir, test_dir in zip(self.inputs, self.local_case_directories):
            if not os.path.isdir(test_dir):
                shutil.copytree(input_dir, test_dir, ignore=ignore_baseline)
            else:
                for f in os.listdir(input_dir):
                    if os.path.isfile(f):
                        shutil.copy2(os.path.join(input_dir, f), test_dir)

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
        Builds the necessary directories
        """
        self.inputs = []
        self.outputs = []
        self.test_build = []
        for i, case in enumerate(self.case):
            if CASE_MAP[case] == "beamdyn":
                self.inputs.append(os.path.join(self.rtest, "modules", "beamdyn", case))
            else:
                self.inputs.append(os.path.join(self.module, case))
            self.outputs.append(os.path.join(self.inputs[i], self.output_type))
            self.test_build.append(os.path.join(self.build, "local_results", case))

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

    def _run_single_case(self, ix: str, case: str, test_build: str):
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

        beamdyn = CASE_MAP[case] == "beamdyn"

        if beamdyn:
            exe = self.bd_executable
            case_input = os.path.join(test_build, "bd_driver.inp")
        else:
            exe = self.of_executable
            case_input = os.path.join(test_build, "".join((case, ".fst")))

        self._run_openfast_case(exe, case_input, ix, case, verbose=self.verbose)

    def _run_openfast_cases(self):
        """
        Runs all of the openfast cases in parallel (if defined).
        """
        arguments = list(zip(self.ix, self.cases, self.local_case_directories))
        with Pool(self.jobs) as pool:
            pool.starmap(self._run_single_case, arguments)

    def run(self):
        """
        Function to build the references to ouput directories. If executing
        OpenFAST, then the directories are created if they don't already exist.
        """

        self._build_directory_references()
        self.n_cases = len(self.case)
        self.ix = [f"{i}/{self.n_cases}" for i in range(1, self.n_cases + 1)]
        if self.execution:
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

