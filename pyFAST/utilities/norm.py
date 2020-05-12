"""Provides the norm functions for the regression testing."""


import sys
import argparse
from typing import List

import numpy as np
import numpy.linalg as LA

from .fast_io import load_output
from .utilities import validate_file


def diff(baseline: np.ndarray, test: np.ndarray, abs_val: bool = True) -> np.ndarray:
    """
    Computes the difference between 2 arrays.

    Parameters
    ----------
    baseline : np.ndarray
        Baseline data for comparison.
    test : np.ndarray
        Test produced data.
    abs_val : bool, optional
        Indicator to use the absolute value, by default True

    Returns
    -------
    np.ndarray
        (Absolute) difference between two arrays.
    """

    diff = test - baseline
    return np.absolute(diff) if abs_val else diff


def pass_regression_test(norm: np.ndarray, tolerance: float) -> bool:
    """
    Indicator for if a norm passes the regression test tolerance condition.

    Parameters
    ----------
    norm : np.ndarray
        Output from one of the norms.
    tolerance : float
        Value that all values in `norm` must be less than.

    Returns
    -------
    bool
        Indicator for if the regression test norm passed.
    """
    return norm.max() < tolerance


# The Norms


def max_norm(baseline: np.ndarray, test: np.ndarray) -> np.ndarray:
    """
    Compute the max norm of the difference between baseline and test results.

    Parameters
    ----------
    baseline : np.ndarray
        Baseline data.
    test : np.ndarray
        Test-produced data.
    axis : int, optional
        Axis for computing the norm, by default 0.

    Returns
    -------
    float
        Max norm of the differene betwen baseline and test data.
    """

    return LA.norm(diff(baseline, test, True), np.inf, axis=0)


def l2_norm(baseline: np.ndarray, test: np.ndarray, abs_val: bool = True) -> np.ndarray:
    """
    Compute the L2 norm of the difference between baseline and test results.

    Parameters
    ----------
    baseline : np.ndarray
        Baseline data.
    test : np.ndarray
        Test-produced data.
    axis : int, optional
        Axis for computing the norm, by default 0.

    Returns
    -------
    float
        L2 norm of the differene betwen baseline and test data.
    """

    return LA.norm(diff(baseline, test, abs_val), 2, axis=0)


def relative_l2_norm(baseline: np.ndarray, test: np.ndarray) -> np.ndarray:
    """
    Compute the relative L2 norm of the difference between baseline and test
    results.

    Parameters
    ----------
    baseline : np.ndarray
        Baseline data.
    test : np.ndarray
        Test-produced data.
    abs_val : bool, optional
        Indicator to use absolute value, by default True.
    axis : int, optional
        Axis for computing the norm, by default 0.

    Returns
    -------
    float
        Relative L2 norm of the differene betwen baseline and test data.
    """

    norm_diff = l2_norm(baseline, test, abs_val=False)
    norm_baseline = LA.norm(baseline, 2, axis=0)

    # Replace zeros with a small value before division
    norm_baseline[norm_baseline == 0] = 1e-16

    norm = norm_diff.copy()
    ix_no_diff = norm_baseline >= 1
    norm[ix_no_diff] = norm_diff[ix_no_diff] / norm_baseline[ix_no_diff]
    return norm


def max_norm_over_range(baseline: np.ndarray, test: np.ndarray) -> np.ndarray:
    """
    Compute the maximum norm of the difference between baseline and test
    results over the range of the baseline data.

    Parameters
    ----------
    baseline : np.ndarray
        Baseline data.
    test : np.ndarray
        Test-produced data.
    axis : int, optional
        Axis for computing the norm, by default 0.

    Returns
    -------
    float
        Maximum norm of the differene betwen baseline and test data.
    """

    ranges = diff(baseline.max(axis=0), baseline.min(axis=0), True)
    norm = max_norm(baseline, test)

    ix_no_diff = ranges >= 1
    norm[ix_no_diff] = norm[ix_no_diff] / ranges[ix_no_diff]
    return norm


def calculate_norms(
    baseline: np.ndarray,
    test: np.ndarray,
    norms: List[str] = [
        "max_norm",
        "max_norm_over_range",
        "l2_norm",
        "relative_l2_norm",
    ],
) -> np.ndarray:
    """
    Compute all the listed norm of the difference between baseline and test
    results over the range of the baseline data.

    Parameters
    ----------
    baseline : np.ndarray
        Baseline data.
    test : np.ndarray
        Test-produced data.
    norms : List[str]
        List of norms to compute.
        Default: ["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"]

    Returns
    -------
    np.ndarray
        Maximum norm of the differene betwen baseline and test data.
    """

    norm_results = np.hstack(
        list(NORM_MAP[norm_type](baseline, test).reshape(-1, 1) for norm_type in norms)
    )
    return norm_results


NORM_MAP = {
    "max_norm": max_norm,
    "max_norm_over_range": max_norm_over_range,
    "l2_norm": l2_norm,
    "relative_l2_norm": relative_l2_norm,
}


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculates the regression test norms."
    )
    parser.add_argument(
        "baseline",
        metavar="Baseline",
        type=str,
        nargs=1,
        help="Baseline data file for comparison.",
    )
    parser.add_argument(
        "test", metavar="Test", type=str, nargs=1, help="Test-produced data file."
    )
    parser.add_argument(
        "norms",
        metavar="[<norm-to-compute>]",
        type=str,
        nargs="+",
        choices=["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"],
        help="The norms to be computed.",
    )
    parser.add_argument(
        "-t",
        "-tol",
        dest="tolerance",
        default=1e-5,
        nargs=1,
        metavar="Tolerance",
        help="Tolerance level for pass/fail condition.",
    )
    parser.add_argument(
        "-a",
        "-abs",
        dest="abs_val",
        action="store_false",
        metavar="Absolute-value-flag",
        help="Indicate absolute value of differences.",
    )
    parser.add_argument(
        "-x",
        "-axis",
        dest="axis",
        default=0,
        metavar="Axis",
        type=int,
        nargs=1,
        help="Axis to compute the norm over.",
    )

    args = parser.parse_args()

    validate_file(args.baseline)
    validate_file(args.test)

    baseline, *_ = load_output(args.baseline)
    test, *_ = load_output(args.test)

    results = calculate_norms(baseline, test, args.norms, args.abs_val, args.axis)

    norms_pass = np.array(
        [pass_regression_test(norm, args.tolerance) for norm in results]
    )
    if norms_pass.sum() == norms_pass.shape[0]:
        print(f"All norms pass with a {args.tolerance:.5f} tolerance.")
        sys.exit(0)
    else:
        fail_norms = args.norms[~norms_pass]
        print(f"{fail_norms} did not pass all cases within {args.tolerance}")
        sys.exit(1)
