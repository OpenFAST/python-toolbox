"""Provides the Executor case"""


import argparse
import os
import platform
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from pathlib import Path
from typing import List

import bokeh
import numpy

from pyFAST.utilities import (ignore_baseline, load_output, run_openfast_case,
                              validate_directory, validate_executable,
                              validate_file)

CASE_MAP = {
    '5MW_ITIBarge_DLL_WTurb_WavesIrr': 'regression',
    '5MW_Land_BD_DLL_WTurb': 'regression',
    '5MW_Land_BD_Linear': 'linear',
    '5MW_Land_DLL_WTurb': 'regression',
    '5MW_OC3Mnpl_DLL_WTurb_WavesIrr': 'regression',
    '5MW_OC3Spar_DLL_WTurb_WavesIrr': 'regression',
    '5MW_OC3Trpd_DLL_WSt_WavesReg': 'regression',
    '5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth': 'regression',
    '5MW_OC4Semi_WSt_WavesWN': 'regression',
    '5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti': 'regression',
    'AOC_WSt': 'regression',
    'AOC_YFix_WSt': 'regression',
    'AOC_YFree_WTurb': 'regression',
    'AWT_WSt_StartUpShutDown': 'regression',
    'AWT_WSt_StartUp_HighSpShutDown': 'regression',
    'AWT_YFix_WSt': 'regression',
    'AWT_YFree_WSt': 'regression',
    'AWT_YFree_WTurb': 'regression',
    'Ideal_Beam_Fixed_Free_Linear': 'linear',
    'Ideal_Beam_Free_Free_Linear': 'linear',
    'SWRT_YFree_VS_EDC01': 'regression',
    'SWRT_YFree_VS_EDG01': 'regression',
    'SWRT_YFree_VS_WTurb': 'regression',
    'UAE_Dnwind_YRamp_WSt': 'regression',
    'UAE_Upwind_Rigid_WRamp_PwrCurve': 'regression',
    'WP_Stationary_Linear': 'linear',
    'WP_VSP_ECD': 'regression',
    'WP_VSP_WTurb': 'regression',
    'WP_VSP_WTurb_PitchFail': 'regression',
    'bd_5MW_dynamic': 'beam_dyn',
    'bd_5MW_dynamic_gravity_Az00': 'beam_dyn',
    'bd_5MW_dynamic_gravity_Az90': 'beam_dyn',
    'bd_curved_beam': 'beam_dyn',
    'bd_isotropic_rollup': 'beam_dyn',
    'bd_static_cantilever_beam': 'beam_dyn',
    'bd_static_twisted_with_k1': 'beam_dyn'
 }


class Executor:
    """Base execution class for OpenFast

    Attributes
    ----------
    """

    def __init__(
        self,
        case,
        exec,
        source,
        compiler,
        tolerance=1e-5,
        plot=0,
        execution=False,
        verbose=False,
        jobs=1
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
        tolerance: float, default: 1e-5
            Error tolerance for pass/fail condition.
        plot : int, default: 0
            Flag to include plotting:
             - 0 (default): No plots will be produced
             - 1: All plots will be produced.
             - 2: Only the plots for failing cases will be produced.
            All plots will be output to <path_to_case_name>/results.html.
        execution : bool, default: True
            Flag to run the test case(s). If `False`, ....
        verbose : bool, default: False
            Flag to include system ouptut.
        jobs : int, default: 1
            Maximum number of parallel jobs to run:
             - -1: 80% of maximum number of nodes available
             - >0: Minimum of the number passed and the number of nodes available
        """

        system_map = {"Darwin": "macos", "Linux": "linux", "Windows": "windows"}

        self.case = case
        self.source = Path(source)
        self.executable = Path(executable)
        self.output_type = "-".join((system_map[platform.system()], compiler.lower()))
        self.build = os.path.join(source, "build")
        self.verbose = verbose
        self.execution = execution
        self.tolerance = tolerance
        self.plot = plot
        self.jobs = jobs if jobs != 0 else -1

        self.rtest = os.path.join(sourceDirectory, "reg_tests", "r-test")
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

        _opts = ("gnu", "intel")
        if self.compiler not in _opts:
            raise ValueError(f"Input 'compiler' must be one of {_opts}!")

        _opts = (0, 1, 2)
        if self.plot not in _opts:
            raise ValueError(f"Input 'plot' must be one of {_opts}!")

        if self.jobs < -1:
            raise ValueError("Input 'jobs' cannot be negative!")
        elif self.jobs == -1:
            self.jobs = int(np.ceil(cpu_count() * 0.8))
        elif self.jobs > 0:
            self.jobs = min(self.jobs, cpu_count())

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
                _source = os.path.join(souce, name)
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

    def build_output_directories(self, case_type: List[str]):
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

    def build_directories(self):
        """Builds the necessary directories"""

        self.inputs = []
        self.outputs = []
        self.test_build = []
        case_types = set()
        for i, case in enumerate(self.case):
            self.inputs.append(os.path.join(self.module, case))
            self.outputs.append(os.path.join(self.inputs[i], self.output_type))
            self.test_build.append(os.path.join(self.build, case))
            case_types.add(CASE_MAP["case"])

        self._build_output_directories()


    def _run_single_case(self, case: str, input_file: str, test_build: str):
        """
        Runs a single OpenFAST test case

        Parameters
        ----------
        case : str
            Case name.
        input_file : str
            Input data for test case.
        test_build : str
            Testing build directory.
        """
        beamdyn = True if CASE_MAP[case] == "beamdyn" else False

        if beamdyn:
            case_input = os.path.join(test_build, "bd_driver.inp")
        else:
            case_input = os.path.join(test_build, "".join((case, ".fst")))

        code = run_openfast_case(self.executable, case_input, verbose=self.verbose, beamdyn=beamdyn)
        if code != 0:
            sys.exit("")

    def run_openfast_cases(self):
        """
        Runs all of the openfast cases in parallel (if defined).
        """

        arguments = list(zip(self.case, self.inputs, self.test_build))
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            pool.map(self._run_single_case, arguments)

    def run(self):
        self.build_directories()
        if self.execution:
            self.run_openfast_cases()

    @staticmethod
    def _get_linear_out_files(self):
        """
        .. note:: Not yet implemented but don't want it to fail.
        """
        return None, None

    @staticmethod
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

    @staticmethod
    def _get_beamdyn_out_files(case: str, out_dir: str, baseline_dir: str):
        """
        Reads the beamdyn test output files for regression testing.

        Parameters
        ----------
        case : str
            Case name.
        out_dir : str
            Test build direcory.
        baseline_dir : str
            Target output directory.
        """

        local = os.path.join(out_dir, "bd_driver.out")
        baseline = os.path.join(baseline_dir, "bd_driver.out")

        for f in (local, baseline):
            validate_file(f)

    def read_out_files(self):

        FUNC_MAP = {
            "linear": self._get_linear_out_files,
            "regression": self._get_regression_out_files,
            "beamdyn": self._get_beamdyn_out_files
        }
        linear = (None, None)
        case_list = []
        test_list = []
        baseline_list = []
        for case, out_dir, target in zip(self.case, self.test_build, self.outputs):
            # Process the files
            _func = FUNC_MAP[CASE_MAP[case]]
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

    def




# regression_tests = {
#     "regression": {
#         "AWT_YFix_WSt": "openfast;elastodyn;aerodyn14;servodyn",
#         "AWT_WSt_StartUp_HighSpShutDown": "openfast;elastodyn;aerodyn15;servodyn"
#         "AWT_YFree_WSt": "openfast;elastodyn;aerodyn15;servodyn",
#         "AWT_YFree_WTurb": "openfast;elastodyn;aerodyn14;servodyn",
#         "AWT_WSt_StartUpShutDown": "openfast;elastodyn;aerodyn15;servodyn",
#         "AOC_WSt": "openfast;elastodyn;aerodyn14;servodyn",
#         "AOC_YFree_WTurb": "openfast;elastodyn;aerodyn15;servodyn",
#         "AOC_YFix_WSt": "openfast;elastodyn;aerodyn15;servodyn",
#         "UAE_Dnwind_YRamp_WSt": "openfast;elastodyn;aerodyn14;servodyn",
#         "UAE_Upwind_Rigid_WRamp_PwrCurve": "openfast;elastodyn;aerodyn15;servodyn",
#         "WP_VSP_WTurb_PitchFail": "openfast;elastodyn;aerodyn14;servodyn",
#         "WP_VSP_ECD": "openfast;elastodyn;aerodyn15;servodyn",
#         "WP_VSP_WTurb": "openfast;elastodyn;aerodyn15;servodyn",
#         "SWRT_YFree_VS_EDG01": "openfast;elastodyn;aerodyn15;servodyn",
#         "SWRT_YFree_VS_EDC01": "openfast;elastodyn;aerodyn15;servodyn",
#         "SWRT_YFree_VS_WTurb": "openfast;elastodyn;aerodyn14;servodyn",
#         "5MW_Land_DLL_WTurb": "openfast;elastodyn;aerodyn15;servodyn",
#         "5MW_OC3Mnpl_DLL_WTurb_WavesIrr": "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn",
#         "5MW_OC3Trpd_DLL_WSt_WavesReg": "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn",
#         "5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth": "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn",
#         "5MW_ITIBarge_DLL_WTurb_WavesIrr": "openfast;elastodyn;aerodyn14;servodyn;hydrodyn;map",
#         "5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti": "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map",
#         "5MW_OC3Spar_DLL_WTurb_WavesIrr": "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map",
#         "5MW_OC4Semi_WSt_WavesWN": "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn",
#         "5MW_Land_BD_DLL_WTurb": "openfast;beamdyn;aerodyn15;servodyn",
#     }
#     "linear":{
#         "WP_Stationary_Linear": "openfast;linear;elastodyn;aerodyn15",
#         "Ideal_Beam_Fixed_Free_Linear": "openfast;linear;beamdyn",
#         "Ideal_Beam_Free_Free_Linear": "openfast;linear;beamdyn",
#         "5MW_Land_BD_Linear": "openfast;linear;beamdyn;servodyn",
#     }
#     "beam_dyn": {
#         "bd_5MW_dynamic": "beamdyn;dynamic",
#         "bd_5MW_dynamic_gravity_Az00": "beamdyn;dynamic",
#         "bd_5MW_dynamic_gravity_Az90": "beamdyn;dynamic",
#         "bd_curved_beam": "beamdyn;static",
#         "bd_isotropic_rollup": "beamdyn;static",
#         "bd_static_cantilever_beam": "beamdyn;static",
#         "bd_static_twisted_with_k1": "beamdyn;static",
#     }
# }
