import unittest
import glob
import pyFAST
from pyFAST.input_output import FASTInputFile
import os
import numpy as np

MyDir=os.path.join(os.path.dirname(__file__),'example_files')

class Test(unittest.TestCase):
 
    def test_001_read_all(self, DEBUG=True):
        nError=0
        for f in glob.glob(os.path.join(MyDir,'FASTIn*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                obj = FASTInputFile(f)
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

    def test_FASTIn(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_BD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['PitchK'],2.0e+07)
        self.assertEqual(F['MemberGeom'][-1,2],61.5)
        self.assertEqual(F['MemberGeom'][-2,3],0.023000)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['RotSpeed'],0.2)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED_bld.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['BldEdgSh(6)'],-0.6952)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ED_twr.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['AdjFASt'],1)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_AD15.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertTrue(F['TipLoss'])

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_ExtPtfm_SubSef.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['StiffnessMatrix'][2,2],1.96653266e+09)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_HD.dat'))
        #F.test_ascii(bCompareWritesOnly=True,bDelete=True) # TODO
        self.assertEqual(F['RdtnDT'],0.0125)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_IF_NoHead.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['Z0'],0.03)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_SbD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['Joints'][0,3],-100)
        self.assertEqual(int(F['Members'][0,1]),1)
        self.assertEqual(int(F['Members'][0,2]),2)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_SD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['PitManRat(1)'],2)

if __name__ == '__main__':
    #Test().test_FASTIn()
    unittest.main()
