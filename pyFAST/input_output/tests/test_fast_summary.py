import unittest
import glob
import pyFAST
from pyFAST.input_output import FASTSummaryFile
import os
import numpy as np
MyDir=os.path.join(os.path.dirname(__file__),'example_files')

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        nError=0
        for f in glob.glob(os.path.join(MyDir,'FASTSum*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                obj = FASTSummaryFile(f)
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

    def test_FASTSum(self):
        f = FASTSummaryFile(os.path.join(MyDir, 'FASTSum_Pendulum.SD.sum.yaml'))
        np.testing.assert_almost_equal(f['CB_frequencies'].ravel(),[2.571561E-02,5.154897E+00,3.448768E+01,3.639185E+01,9.826435E+01], 5)
        pass


if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
