import unittest
import glob
import pyFAST
from pyFAST.input_output import FASTLinearizationFile
import os
import numpy as np
MyDir=os.path.join(os.path.dirname(__file__),'example_files')

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        nError=0
        for f in glob.glob(os.path.join(MyDir,'FASTLin*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                obj = FASTLinearizationFile(f)
                s=type(obj).__name__.replace('file','')[:20]
                if DEBUG:
                    print('[ OK ] {:30s}\t{:20s}'.format(os.path.basename(f)[:30],s))
            except:
                nError += 1
                if DEBUG:
                    print('[FAIL] {:30s}\tException occurred'.format(os.path.basename(f)[:30]))
                raise 
        if nError>0:
            raise Exception('Some tests failed')

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
