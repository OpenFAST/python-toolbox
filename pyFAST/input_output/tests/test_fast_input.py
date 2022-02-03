import unittest
import os
import numpy as np

from .helpers_for_test import MyDir, reading_test 

import pyFAST
from pyFAST.input_output import FASTInputFile
from pyFAST.input_output.fast_wind_file import FASTWndFile


class Test(unittest.TestCase):
 
    def test_001_read_all(self, DEBUG=True):
        reading_test('FASTIn*.*', FASTInputFile)

    def test_FASTIn(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_BD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(F['PitchK'],2.0e+07)
        self.assertEqual(F['MemberGeom'][-1,2],61.5)
        self.assertEqual(F['MemberGeom'][-2,3],0.023000)

        F=FASTInputFile(os.path.join(MyDir,'FASTIn_BD_bld.dat'))
        F.test_ascii(bCompareWritesOnly=False,bDelete=True)
        self.assertEqual(F['DampingCoeffs'][0][0],0.01)
        # TODO BeamDyn Blade properties are not really "user friendly"
        self.assertEqual(F['BeamProperties']['span'][1],1.0)
        self.assertEqual(F['BeamProperties']['K'][1][0,0],1.8e+08) # K11 @ section 2
        self.assertEqual(F['BeamProperties']['M'][1][0,0],1.2) # M11 @ section 2

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
        
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_MD.dat'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)
        self.assertEqual(float(F['LineTypes'][0,1]),0.02)

    def test_FASTWnd(self):
        F=FASTWndFile(os.path.join(MyDir,'FASTWnd.wnd'))
        F.test_ascii(bCompareWritesOnly=True,bDelete=True)

    def test_FASTInGraph(self):
        F=FASTInputFile(os.path.join(MyDir,'FASTIn_HD.dat'))
        #graph = F.toGraph()
        #print(graph)
        #self.assertEqual(len(graph.Nodes), 4)
        #self.assertEqual(len(graph.Elements), 3)
# 
        #F=FASTInputFile(os.path.join(MyDir,'FASTIn_SbD.dat'))
        #print(F)
        #graph = F.toGraph()
#         self.assertEqual(len(graph.Nodes), 2)
#         self.assertEqual(len(graph.Elements), 1)

if __name__ == '__main__':
    #Test().test_FASTIn()
    unittest.main()
