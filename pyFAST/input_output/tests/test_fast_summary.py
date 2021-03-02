import unittest
import os
import numpy as np

from .helpers_for_test import MyDir, reading_test 
import pyFAST
from pyFAST.input_output import FASTSummaryFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTSum*.*', FASTSummaryFile)

    def test_FASTSum(self):
        f = FASTSummaryFile(os.path.join(MyDir, 'FASTSum_Pendulum.SD.sum.yaml'))
        np.testing.assert_almost_equal(f['CB_frequencies'].ravel(),[2.571561E-02,5.154897E+00,3.448768E+01,3.639185E+01,9.826435E+01], 5)
        pass


if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
