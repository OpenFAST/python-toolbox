
from typing import List, Tuple
import numpy as np
from functools import partial
from multiprocessing.pool import Pool
from .utilities.norm import calculate_norms, pass_regression_test


class RegressionTester:

    def __init__(self, tolerance):
        self.jobs = 1
        self.tolerance = tolerance

    def test_norm(
        self,
        ix_list: List[str],
        case_list: List[str],
        baseline_list: List[Tuple[np.ndarray, list]],
        test_list: List[Tuple[np.ndarray, list]],
        norm_list: List[str] = ["max_norm", "max_norm_over_range", "l2_norm", "relative_l2_norm"],
        test_norm_condition: List[str] = ["relative_l2_norm"],  # flag in __main__.py
    ) -> Tuple[List[np.ndarray], List[bool], List[str]]:
        """
        Computes the norms for each of the valid test cases.

        Parameters
        ----------
        ix_list : List[str]
            List of indices for cases to be run in the form of "i/N".
        case_list : List[str]
            List of valid cases where a norm can be computed.
        baseline_list : List[Tuple[np.ndarray, list]]
            Tuples of baseline data and info corresponding to `case_list.
        test_list : List[Tuple[np.ndarray, list]]
            Tuples of test data and info correpsonding to `case_list`.
        norm_list : List[str], optional
            List of norms to be computed, by default:
            ["max_norm","max_norm_over_range","l2_norm","relative_l2_norm"]
        test_norm_condition : List[str]
            Defines which norm(s) to use for the pass/fail condition, by
            default ["relative_l2_norm"].
        jobs : int
            Number of parallel jobs to compute the norms and test them.

        Returns
        -------
        norm_results : List[np.ndarray]
            List of norm results corresponding to `case_list`. Each array will
            have shape [len(attributes), len(norm_list)]
        pass_fail_list : List[bool]
            A list of indicators for if the case passed the regression test.
        norm_list : List[str]
            `norm_list`.
        """

        # Test to make sure that the test norm is included in the computed norms
        if not set(test_norm_condition).issubset(norm_list):
            message = (
                f"test_norm_condition: {test_norm_condition} should be contained in "
                f"norm_list: {norm_list}."
            )
            raise ValueError(message)

        norm_results = []
        arguments = [[b[0], t[0]] for b, t in zip(baseline_list, test_list)]
        partial_norms = partial(calculate_norms, norms=norm_list)
        with Pool(self.jobs) as pool:
            norm_results = list(pool.starmap(partial_norms, arguments))

        norm_ix = [norm_list.index(norm) for norm in test_norm_condition]
        arguments = [(norm[:, norm_ix], self.tolerance) for norm in norm_results]
        with Pool(self.jobs) as pool:
            pass_fail_list = list(pool.starmap(pass_regression_test, arguments))

        n_fail = len(pass_fail_list) - sum(pass_fail_list)
        fail = []
        for ix, case, _pass in zip(ix_list, case_list, pass_fail_list):
            ix = ix.split("/")[0]
            if not _pass:
                fail.append((ix, case))
            pf = "pass" if _pass else "\x1b[1;31mFAIL\x1b[0m"
            print(f"{ix.rjust(6)}  TEST: {case.ljust(40, '.')} {pf}")

        print(f"\n{str(n_fail).rjust(6)} cases \x1b[1;31mFAILED\x1b[0m")
        for ix, case in fail:
            ix = ix.split("/")[0]
            message = f"\x1b[1;31m{ix.rjust(8)} {case}\x1b[0m"
            print(message)

        return norm_results, pass_fail_list, norm_list    

