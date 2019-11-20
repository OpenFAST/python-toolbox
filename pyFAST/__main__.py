"""CLI functionality for pyfast"""

import os
import sys
import argparse

from pyfast import Executor, utilities
from case_list import CASE_LIST


def match_cases(case_regex):
    case_regex = case_regex.lower()
    if case_regex == "all":
        return case_regex
    return [c for c in CASE_LIST if case_regex in c.lower()]


parser = argparse.ArgumentParser(
    description="Calculates the regression test norms.", prog="pyfast"
)

parser.add_argument(
    "-r",
    "--regex-case",
    dest="case",
    type=str,
    nargs="+",
    required=True,
    help='"Regex" case names where the text. Looks to see if the provided string is contained in any of the valid cases. Note: not case sensitive',
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
    default=0,
    help="0: no plots; 1: all plots; 2: plot failure cases only",
)
parser.add_argument(
    "-j",
    "--jobs",
    dest="jobs",
    type=int,
    default=-1,
    help="Number of cases to run in parallel. Use -1 for 80% of available cores",
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
cases = [c for el in args.case for c in match_cases(el)]

reg_test = Executor(
    cases,
    args.executable,
    args.source,
    args.compiler,
    tolerance=args.tolerance,
    plot=args.plot,
    execution=~args.no_execution,
    verbose=args.verbose,
    jobs=args.jobs,
)

# Run openFAST cases
reg_test.run()

# Gather the outputs
cases, baseline, test = reg_test.read_out_files()

# Check for pass/fail

# Plot the results
