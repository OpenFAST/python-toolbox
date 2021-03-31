import unittest
import os
import numpy as np
import pyFAST
# import pyFAST.lin import TurbSimFile
import pyFAST.linearization.mbc.mbc3 as mbc


MyDir=os.path.join(os.path.dirname(__file__))

class Test(unittest.TestCase):

    def mbc3_standstill(self, lin_file):
        # Script Parameters
        BladeLen     = 40.04                # Blade length, used to tune relative modal energy [m]
        TowerLen     = 55.59                # Tower length, used to tune relative modal energy [m]

        # Derived parameters
        lin_files = np.array([lin_file])

        # Performing MBC (NOTE: not stricly necessary without rotation)
        mbc_data, matData, FAST_linData = mbc.fx_mbc3(lin_files, verbose=False)
        CD = mbc.campbell_diagram_data(mbc_data,BladeLen,TowerLen)

        nModesMax = np.min([len(CD['Modes']),10])
        Freq = np.array([CD['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
        Damp = np.array([CD['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])
        LogDec = Damp*100*2*np.pi
        return Freq, Damp, LogDec

    def test_mbc3_standstill(self):
        lin_file     = os.path.join(MyDir,'../../../../data/example_files/Standstill.1.lin') # Linearization file
        Freq, Damp, LogDec = self.mbc3_standstill(lin_file)
        np.testing.assert_almost_equal(Freq[:3]  ,[0.427, 0.450, 0.669], 3)
        np.testing.assert_almost_equal(LogDec[:3],[1.9505,2.1309,5.0649], 4)

    def test_mbc3_standstill_old(self):
        lin_file = os.path.join(MyDir, '../../../../data/example_files/Standstill_old.1.lin') # Linearization file
        Freq, Damp, LogDec = self.mbc3_standstill(lin_file)
        np.testing.assert_almost_equal(Freq[:3]  ,[0.427, 0.449, 0.667], 3)
        np.testing.assert_almost_equal(LogDec[:3],[1.9497,2.1162,5.0113], 4)



if __name__ == '__main__':
#     Test().test_000_debug()
    unittest.main()
