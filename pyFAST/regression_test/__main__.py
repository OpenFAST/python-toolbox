"""CLI functionality for pyFAST"""

import sys
import argparse
from typing import List
import re

from pyFAST import (
    CASE_MAP,
    Executor,
    RegressionTester,
    SummaryHandler
)


def match_cases(case_regex: str) -> List[str]:
    """
    Match provided case with available test cases.

    Parameters
    ----------
    case_regex : str
        Case to match against available test cases.

    Returns
    -------
    List[str]
        List of cases that have `case_regex` as a case-sensitive substring.
    """

    prog = re.compile(case_regex)
    return [c for c in CASE_MAP if prog.match(c)]


def main():
    """
    Runs the pyFAST suite.
    """

    parser = argparse.ArgumentParser(
        description="Executes the requested cases and quantifies the difference in local and baseline results.",
        prog="pyFAST"
    )

    parser.add_argument(
        "-l",
        "--list-cases",
        dest="list",
        action="store_true",
        help=('Display the available test cases.'),
    )
    _list_flag_given = '-l' in sys.argv or '--list' in sys.argv

    parser.add_argument(
        "-r",
        "--regex-cases",
        dest="cases",
        type=str,
        nargs="+",
        default=[],
        help=(
            'Case sensitive regex of case names to run from the list of available cases.'
        ),
    )
    parser.add_argument(
        "-e",
        dest="executable",
        type=str,
        nargs="+",
        required=not _list_flag_given,
        help="Test-produced data file.",
    )
    parser.add_argument(
        "--openfast-root",
        dest="openfast_root",
        type=str,
        required=not _list_flag_given,
        help="Path to the OpenFAST repository",
    )
    parser.add_argument(
        "--system",
        dest="system",
        type=str.lower,
        choices=["macos", "linux", "windows"],
        help="Operating system to use for baseline results.",
    )
    parser.add_argument(
        "-c",
        dest="compiler",
        type=str.lower,
        required=not _list_flag_given,
        choices=["intel", "gnu"],
        help="Compiler ID for the system",
    )
    parser.add_argument(
        "-t",
        "--tol",
        dest="tolerance",
        type=float,
        default=1e-5,
        help="Tolerance for determing the failure of a norm",
    )
    parser.add_argument(
        "-p",
        "--plot",
        dest="plot",
        choices=[0, 1, 2],
        type=int,
        default=0,
        help="0: no plots; 1: all plots; 2: plot failure cases only",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        dest="jobs",
        type=int,
        default=-1,
        help="Number of cases to run in parallel. Use -1 for 80 pct of available cores",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Use flag to use verbose output.",
    )
    parser.add_argument(
        "-n",
        "--no-execution",
        dest="no_execution",
        action="store_true",
        help="Avoid executing the OpenFAST case simulations. The simulation results should already exist and the results comparison will be performed.",
    )
    # parser.add_argument(
    #     "--norm",
    #     dest="norm_list",
    #     type=str.lower,
    #     nargs="+",
    #     default=["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"],
    #     choices=["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"],
    #     help="The norm(s) to be computed.",
    # )
    # parser.add_argument(
    #     "--test-norm",
    #     dest="test_norm",
    #     type=str.lower,
    #     nargs="+",
    #     default=["relative_l2_norm"],
    #     choices=["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"],
    #     help="Norm(s) used to determine if the test(s) pass. Must be a normed passed to `-norm`.",
    # )

    args = parser.parse_args()

    if args.list:
        # Display all available cases and return
        print("OpenFAST cases:")
        for case in [key for key, value in CASE_MAP.items() if value["driver"] == "openfast"]:
            print(f"\t{case}")
        print("")
        print("BeamDyn module cases:")
        for case in [key for key, value in CASE_MAP.items() if value["driver"] == "beamdyn"]:
            print(f"\t{case}")
        return

    # Get the cases given by the reg ex, if requested
    cases = args.cases if not args.cases else [c for el in args.cases for c in match_cases(el)]
    cases = list(set(cases))

    executor = Executor(
        cases,
        args.executable,
        args.openfast_root,
        args.compiler,
        system=args.system,
        no_execution=args.no_execution,
        verbose=args.verbose,
        jobs=args.jobs,
    )

    # Run cases
    executor.run()

    # Gather the outputs
    baseline, test = executor.read_output_files()

    # Run the regression test
    reg_test = RegressionTester(args.tolerance)
    ix = [f"{i}/{len(executor.cases)}" for i in range(1, len(executor.cases) + 1)]
    norm_res, pass_fail_list, norm_list = reg_test.test_norm(
        ix,
        cases,
        baseline,
        test
    )

    # Extract the attributes metadata and the data
    attributes = [
        list(zip(info["attribute_names"], info["attribute_units"]))
        for _, info in baseline
    ]
    baseline_data = [data for data, _ in baseline]
    test_data = [data for data, _ in test]

    # Create the regression test summaries
    summary = SummaryHandler(args.plot, executor.local_test_location)
    plots = summary.retrieve_plot_html(baseline_data, test_data, attributes, pass_fail_list)
    summary.create_results_summary(cases, attributes, norm_res, norm_list, plots, args.tolerance)


if __name__ == "__main__":

    main()
