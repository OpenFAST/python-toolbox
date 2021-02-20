import unittest
import os
import numpy as np
from .helpers_for_test import MyDir, reading_test 
import pyFAST
from pyFAST.input_output import FASTLinearizationFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTLin*.*', FASTLinearizationFile)

    def test_FASTLin(self):
        F=FASTLinearizationFile(os.path.join(MyDir,'FASTLin.lin'))
        self.assertAlmostEqual(F['A'][3,1], 3.91159454E-04 )
        self.assertAlmostEqual(F['u'][7]   ,4.00176055E+04)

        F=FASTLinearizationFile(os.path.join(MyDir,'FASTLin_EDM.lin'))
        dfs=F.toDataFrame()
        M=dfs['M']
        self.assertAlmostEqual(M['7_TwFADOF1']['7_TwFADOF1'],0.436753E+06)
        self.assertAlmostEqual(M['13_GeAz']['13_GeAz']     , 0.437026E+08)

if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
