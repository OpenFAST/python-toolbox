import unittest
import os
import numpy as np
from .helpers_for_test import MyDir, reading_test 

import pyFAST
from pyFAST.input_output import FASTOutputFile

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTOut*.*', FASTOutputFile)

    def DF(self,FN):
        """ Reads a file and return a dataframe """ 
        return FASTOutputFile(os.path.join(MyDir,FN)).toDataFrame()
 
    def test_FASTOut(self):
        self.assertEqual(self.DF('FASTOut.out').values[-1,1],1036)
 
    def test_FASTOutBin(self):
        M = self.DF('FASTOutBin.outb')
        self.assertAlmostEqual(M['GenPwr_[kW]'].values[-1],40.57663190807828)

if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
