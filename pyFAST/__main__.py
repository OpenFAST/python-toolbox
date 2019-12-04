"""CLI functionality for pyfast"""

import argparse
from typing import List

from pyfast import CASE_LIST, Executor


def match_cases(case_regex: str) -> List[str]:
    """Match provided case with available test cases.

    Parameters
    ----------
    case_regex : str
        Case to match against available test cases.

    Returns
    -------
    List[str]
        List of cases that have `case_regex` as a non-case-sensitive substring.
    """

    case_regex = case_regex.lower()
    return [c for c in CASE_LIST if case_regex in c.lower()]


def main():
    cases = (
        "all"
        if "all" in args.case
        else [c for el in args.case for c in match_cases(el)]
    )

    execution = not args.no_execution
    reg_test = Executor(
        cases,
        args.executable,
        args.source,
        args.compiler,
        tolerance=args.tolerance,
        plot=args.plot,
        execution=execution,
        verbose=args.verbose,
        jobs=args.jobs,
    )

    # Run openFAST cases
    reg_test.run()

    # Gather the outputs
    ix, cases, baseline, test = reg_test.read_out_files()

    # Run the regression test
    norm_res, pass_fail_list, norm_list = reg_test.test_norm(ix, cases, baseline, test)

    # Extract the attributes metadata and the data
    attributes = [
        list(zip(info["attribute_names"], info["attribute_units"]))
        for _, info in baseline
    ]
    baseline_data = [data for data, _ in baseline]
    test_data = [data for data, _ in test]

    # Create the case summary for each case
    plots = reg_test.retrieve_plot_html(
        cases, baseline_data, test_data, attributes, pass_fail_list
    )
    reg_test.create_results_summary(cases, attributes, norm_res, norm_list, plots)


parser = argparse.ArgumentParser(
    description="Calculates the regression test norms.", prog="pyfast"
)

parser.add_argument(
    "-r",
    "--regex-case",
    dest="case",
    type=str,
    nargs="+",
    default=["all"],
    required=False,
    help=(
        '"Regex" case names where the text. Looks to see if the provided'
        "string is contained in any of the valid cases. "
        "Note: not case sensitive"
    ),
)
parser.add_argument(
    "-e",
    dest="executable",
    type=str,
    nargs="+",
    required=True,
    help="Test-produced data file.",
)

parser.add_argument(
    "-s", dest="source", type=str, required=True, help="Path to the OpenFAST repository"
)
parser.add_argument(
    "-c",
    dest="compiler",
    type=str,
    required=True,
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
    help="Doesn't run the actual test case(s). The data should already exist.",
)
# Doesn't run openfast, uses pre-compiled results

parser.add_argument(
    "--norm",
    dest="norm_list",
    type=str,
    nargs="+",
    choices=["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"],
    help="The norm(s) to be computed.",
)

# Parse the arguments, find the cases to be run, and initialize the openFAST
# execution class
args = parser.parse_args()

main()
