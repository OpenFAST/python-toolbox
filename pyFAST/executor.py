"""Provides the Executor case"""


import os
import sys
import shutil
import platform
from typing import List, Tuple
from pathlib import Path
from multiprocessing import cpu_count
from multiprocessing.pool import Pool

import bokeh
import numpy as np
from pyfast.utilities import (
    load_output,
    validate_file,
    calculate_norms,
    ignore_baseline,
    run_openfast_case,
    validate_directory,
    validate_executable,
)

SYSTEM_MAP = {"Darwin": "macos", "Linux": "linux", "Windows": "windows"}

CASE_MAP = {
    "5MW_ITIBarge_DLL_WTurb_WavesIrr": "regression",
    "5MW_Land_BD_DLL_WTurb": "regression",
    "5MW_Land_BD_Linear": "linear",
    "5MW_Land_DLL_WTurb": "regression",
    "5MW_OC3Mnpl_DLL_WTurb_WavesIrr": "regression",
    "5MW_OC3Spar_DLL_WTurb_WavesIrr": "regression",
    "5MW_OC3Trpd_DLL_WSt_WavesReg": "regression",
    "5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth": "regression",
    "5MW_OC4Semi_WSt_WavesWN": "regression",
    "5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti": "regression",
    "AOC_WSt": "regression",
    "AOC_YFix_WSt": "regression",
    "AOC_YFree_WTurb": "regression",
    "AWT_WSt_StartUpShutDown": "regression",
    "AWT_WSt_StartUp_HighSpShutDown": "regression",
    "AWT_YFix_WSt": "regression",
    "AWT_YFree_WSt": "regression",
    "AWT_YFree_WTurb": "regression",
    "Ideal_Beam_Fixed_Free_Linear": "linear",
    "Ideal_Beam_Free_Free_Linear": "linear",
    "SWRT_YFree_VS_EDC01": "regression",
    "SWRT_YFree_VS_EDG01": "regression",
    "SWRT_YFree_VS_WTurb": "regression",
    "UAE_Dnwind_YRamp_WSt": "regression",
    "UAE_Upwind_Rigid_WRamp_PwrCurve": "regression",
    "WP_Stationary_Linear": "linear",
    "WP_VSP_ECD": "regression",
    "WP_VSP_WTurb": "regression",
    "WP_VSP_WTurb_PitchFail": "regression",
    "bd_5MW_dynamic": "beam_dyn",
    "bd_5MW_dynamic_gravity_Az00": "beam_dyn",
    "bd_5MW_dynamic_gravity_Az90": "beam_dyn",
    "bd_curved_beam": "beam_dyn",
    "bd_isotropic_rollup": "beam_dyn",
    "bd_static_cantilever_beam": "beam_dyn",
    "bd_static_twisted_with_k1": "beam_dyn",
}


def _get_linear_out_files():
    """
    .. note:: Not yet implemented but don't want it to fail.
    """
    return None, None


def _get_regression_out_files(case: str, out_dir: str, baseline_dir: str):
    """
    Reads the standard regression test output files for regression testing.

    Parameters
    ----------
    case : str
        Case name.
    out_dir : str
        Test build direcory.
    baseline_dir : str
        Target output directory.
    """

    case_file = "".join((case, ".outb"))
    local = os.path.join(out_dir, case_file)
    baseline = os.path.join(baseline_dir, case_file)

    for f in (local, baseline):
        validate_file(f)

    return local, baseline


def _get_beamdyn_out_files(case: str, out_dir: str, baseline_dir: str):
    """
    Reads the beamdyn test output files for regression testing.

    Parameters
    ----------
    case : str
        Case name. Not used.
    out_dir : str
        Test build direcory.
    baseline_dir : str
        Target output directory.
    """

    local = os.path.join(out_dir, "bd_driver.out")
    baseline = os.path.join(baseline_dir, "bd_driver.out")

    for f in (local, baseline):
        validate_file(f)

    return local, baseline


class Executor:
    """Base execution class for OpenFast

    Attributes
    ----------
    """

    FUNC_MAP = {
        "linear": _get_linear_out_files,
        "regression": _get_regression_out_files,
        "beamdyn": _get_beamdyn_out_files,
    }

    def __init__(
        self,
        case: str,
        executable: str,
        source: str,
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
        case : list(str, ...)
            Test case name(s) as a list of strings.
        exececutable : str
            Path to the OpenFAST executable.
        source : str
            Path to OpenFAST repository.
        compiler : str
            System compiler id. Should be one of "intel" or "gnu".
        system : str
            Operating system version of results used for comparison, default
            None (machine's OS). Should be one of "Windows", "Linux", or
            "Darwin".
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

        system = SYSTEM_MAP[platform.system() if system is None else system]

        if case == "all":
            self.case = [*CASE_MAP]
        else:
            self.case = case if isinstance(case, list) else [case]
        self.compiler = compiler
        self.output_type = "-".join((system, self.compiler.lower()))

        self.source = Path(source)
        self.executable = Path(executable)
        self.build = os.path.join(self.source, "build")

        self.verbose = verbose
        self.execution = execution
        self.tolerance = tolerance
        self.plot = plot
        self.plot_path = plot_path
        self.jobs = jobs if jobs != 0 else -1

        self.rtest = os.path.join(self.source, "reg_tests", "r-test")
        self.module = os.path.join(self.rtest, "glue-codes", "openfast")

        self._validate_inputs()

    def _validate_inputs(self):
        """Method to ensure inputs are valid."""

        _opts = ("macos-gnu", "linux-intel", "linux-gnu", "windows-intel")
        if self.output_type not in _opts:
            self.output_type = "macos-gnu"
            print(f"Defaulting to {self.output_type} for output type")

        validate_executable(self.executable)
        validate_directory(self.build)

        _opts = (0, 1, 2)
        if self.plot not in _opts:
            raise ValueError(f"Input 'plot' must be one of {_opts}!")

        if self.jobs < -1:
            raise ValueError("Input 'jobs' cannot be negative!")
        elif self.jobs == -1:
            self.jobs = int(np.ceil(cpu_count() * 0.8))
        elif self.jobs > 0:
            self.jobs = min(self.jobs, cpu_count())
        if self.jobs > len(self.case):
            self.jobs = len(self.case)

    def _build_beamdyn_output_directories(self):
        """
        Creates the local output directories for BeamDyn cases and intializes
        it with the input files.
        """

        for case, in_dir, test in zip(self.case, self.inputs, self.test_build):
            if CASE_MAP[case] != "beamdyn":
                continue
            for bd_file in ("bd_driver.inp", "bd_primary.inp", "beam_props.inp"):
                shutil.copy(os.path.join(in_dir, bd_file), test)

    def _build_5MW_directories(self):
        """Copies the 5MW Baseline folder"""

        target = os.path.join(self.build, "5MW_Baseline")
        source = os.path.join(self.module, "5MW_Baseline")
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

    def _build_test_directory(self):
        """Copies the input data to the test build directory"""

        for input_dir, test_dir in zip(self.inputs, self.test_build):
            if not os.path.isdir(test_dir):
                shutil.copytree(input_dir, test_dir, ignore=ignore_baseline)

    def _build_output_directories(self, case_type: List[str]):
        """Creates the local output directories"""

        _linear = ("Ideal_Beam", "WP_Baseline")
        _regression = ("AOC", "AWT27", "SWRT", "UAE_VI", "WP_Baseline")
        directories = []

        if "linear" in case_type:
            directories.extend(_linear)

        if "regression" in case_type:
            directories.extend(_regression)

        if "beamdyn" in case_type:
            self._build_beamdyn_output_directories()
            if directories == []:
                return

        for data in directories:
            _dir = os.path.join(self.build, data)
            if not os.path.isdir(_dir):
                shutil.copytree(os.path.join(self.module, data), _dir)

        self._build_5MW_directories()
        self._build_test_directory()

    def _build_directories(self):
        """Builds the necessary directories"""

        self.inputs = []
        self.outputs = []
        self.test_build = []
        case_types = set()
        for i, case in enumerate(self.case):
            self.inputs.append(os.path.join(self.module, case))
            self.outputs.append(os.path.join(self.inputs[i], self.output_type))
            self.test_build.append(os.path.join(self.build, case))
            case_types.add(CASE_MAP[case])

    def _run_single_case(self, args: List[Tuple[str, str]]):
        """
        Runs a single OpenFAST test case

        Parameters
        ----------
        args : List[Tuple[str, str]]
            ix : str
                String index as "i/n" for which case is being run.
            case : str
                Case name.
            test_build : str
                Testing build directory.
        """

        ix, case, test_build = args
        beamdyn = CASE_MAP[case] == "beamdyn"

        if self.verbose:
            print(f"{ix.rjust(4)}: {case} running")

        if beamdyn:
            case_input = os.path.join(test_build, "bd_driver.inp")
        else:
            case_input = os.path.join(test_build, "".join((case, ".fst")))

        code = run_openfast_case(
            self.executable, case_input, verbose=self.verbose, beamdyn=beamdyn
        )
        if code != 0:
            sys.exit("")

        if self.verbose:
            print(f"{ix.rjust(4)}: {case} complete")

    @staticmethod
    def _run_case_logging(ix: List[str], case: List[str]):
        """
        Prints the cases to be run.

        Parameters
        ----------
        ix : List[str]
            Index with format i/n.
        case : List[str]
            List of cases.
        """

        print("    Cases to be run:\n")
        for i, c in zip(ix, case):
            print(f"{i.rjust(6)}: {c}")
        print()

    def _run_openfast_cases(self):
        """
        Runs all of the openfast cases in parallel (if defined).
        """

        n_cases = len(self.case)
        ix = [f"{i}/{n_cases}" for i in range(1, n_cases + 1)]
        if self.verbose:
            self._run_case_logging(ix, self.case)

        arguments = list(zip(ix, self.case, self.test_build))
        with Pool(self.jobs) as pool:
            pool.map(self._run_single_case, arguments)

    def run(self):
        """
        Function to build the references to ouput directories. If executing
        OpenFAST, then the directories are created if they don't already exist.
        """

        self._build_directories()
        if self.execution:
            case_types = list(set(CASE_MAP[c] for c in self.case))
            self._build_output_directories(case_types)
            self._run_openfast_cases()

    def read_out_files(self):
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

        case_list = []
        test_list = []
        baseline_list = []
        for case, out_dir, target in zip(self.case, self.test_build, self.outputs):
            # Process the files
            _func = self.FUNC_MAP[CASE_MAP[case]]
            local, baseline = _func(case, out_dir, target)

            # Check for linear case
            if local is None and baseline is None:
                continue

            # Extract the data
            test_data, test_info, _ = load_output(local)
            baseline_data, baseline_info, _ = load_output(baseline)
            case_list.append(case)
            test_list.append((test_data, test_info))
            baseline_list.append((baseline_data, baseline_info))

        return case_list, baseline_list, test_list

    @staticmethod
    def test_norm(
        case_list: List[str],
        baseline_list: List[Tuple[np.ndarray, list]],
        test_list: List[Tuple[np.ndarray, list]],
        norm_list: List[str] = [
            "max_norm",
            "max_norm_over_range",
            "l2_norm",
            "relative_l2_norm",
        ],
    ) -> List[np.ndarray]:
        """
        Computes the norms for each of the valid test cases.

        Parameters
        ----------
        case_list : List[str]
            List of valid cases where a norm can be computed.
        baseline_list : List[Tuple[np.ndarray, list]]
            Tuples of baseline data and info corresponding to `case_list.
        test_list : List[Tuple[np.ndarray, list]]
            Tuples of test data and info correpsonding to `case_list`.
        norm_list : List[str], optional
            List of norms to be computed, by default ["max_norm","max_norm_over_range","l2_norm","relative_l2_norm",]

        Returns
        -------
        norm_results : List[np.ndarray]
            List of norm results corresponding to `case_list`. Each array will
            have shape [len(attributes), len(norm_list)]
        """

        norm_results = []
        for case in case_list:
            norm_results.append(calculate_norms(baseline, test, norm_list))

        return norm_results

    def create_results_summary(
        self, case_list: List[str], norm_results: List[np.ndarray]
    ):
        """
        Creates the results summary html file for each case in `case_list`.

        Parameters
        ----------
        case_list : List[str]
            List of cases to produce html summaries.
        norm_results : List[np.ndarray]
            Computed norms for each of the cases in `case_list`.
        """

        if self.plot_path is None:
            self.plot_path = [self.test_build[self.case.index(c)] for c in case_list]

        #### LEFT OFF HERE ####
        create_case_summary(
            path, case, results, results_max, tolerance, plots, results_columns
        )
