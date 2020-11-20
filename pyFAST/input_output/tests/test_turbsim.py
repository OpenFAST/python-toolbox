import unittest
import glob
import pyFAST
from pyFAST.input_output import TurbSimFile
import os
import numpy as np
MyDir=os.path.join(os.path.dirname(__file__),'example_files')

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        nError=0
        for f in glob.glob(os.path.join(MyDir,'TurbSim*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                obj = TurbSimFile(f)
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

    def test_TurbSim(self):
        # --- Test without tower
        F = TurbSimFile(os.path.join(MyDir,'TurbSimNoTwr.bts'))
        F.write(      os.path.join(MyDir,'TurbSimNoTwr_TMP.bts'))
        F2= TurbSimFile(os.path.join(MyDir,'TurbSimNoTwr_TMP.bts'))
        os.remove(    os.path.join(MyDir,'TurbSimNoTwr_TMP.bts'))
        np.testing.assert_almost_equal(F['u'][0,:,:,:],F2['u'][0,:,:,:],4)
        np.testing.assert_almost_equal(F['u'][1,:,:,:],F2['u'][1,:,:,:],4)
        np.testing.assert_almost_equal(F['u'][2,:,:,:],F2['u'][2,:,:,:],4)
        # --- Test with tower
        F = TurbSimFile(os.path.join(MyDir,'TurbSimWithTwr.bts'))
        np.testing.assert_almost_equal(F['u'][2,-1,1,3], 0.508036, 5)
        np.testing.assert_almost_equal(F['u'][0, 4,2,0], 7.4867466, 5)
        np.testing.assert_almost_equal(F['uTwr'][0, 4, :], [6.1509, 6.4063, 8.9555, 7.6943], 4)
        F.write(      os.path.join(MyDir,'TurbSimWithTwr_TMP.bts'))
        F2= TurbSimFile(os.path.join(MyDir,'TurbSimWithTwr_TMP.bts'))
        os.remove(    os.path.join(MyDir,'TurbSimWithTwr_TMP.bts'))
        np.testing.assert_almost_equal(F['u'][0,:,:,:],F2['u'][0,:,:,:],3)
        np.testing.assert_almost_equal(F['u'][1,:,:,:],F2['u'][1,:,:,:],3)
        np.testing.assert_almost_equal(F['u'][2,:,:,:],F2['u'][2,:,:,:],3)

if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
