import unittest
import glob
import pyFAST
from pyFAST.input_output import FASTOutputFile
import os
import numpy as np
MyDir=os.path.join(os.path.dirname(__file__),'example_files')

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        nError=0
        for f in glob.glob(os.path.join(MyDir,'FASTOut*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                obj = FASTOutputFile(f)
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
