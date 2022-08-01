import unittest
import os
import numpy as np    

from pyFAST.fastfarm import *

MyDir=os.path.dirname(__file__)

class Test(unittest.TestCase):

    def test_box_extent(self):
        # Test the turbulence box extent function

        # --- TurbSim Box paarameters
        yBox  = np.arange(-187.5,188  ,5  )
        zBox  = np.arange(1     ,282  ,5  )
        tBox  = np.arange(0     ,4.901,0.1)
        meanU = 5.968230471458111
        # --- Parameters for TurbSim Extent
        D             = 77.0                       # Turbine diameter (m)
        hubHeight     = 78.045                     # Hub Height (m)
        extent_X_high = 1.2                        # x-extent of high res box in diamter around turbine location
        extent_Y_high = 1.2                        # y-extent of high res box in diamter around turbine location
        chord_max     = 5                          # maximum blade chord (m). Turbine specific.
        Cmeander      = 1.9                        # Meandering constant (-)
        # --- Layout
        xWT = [0.0, 265.]    # x positions of turbines
        yWT = [0.0, 50.0] # y postitions of turbines
        zWT = [0.0, 0.0 ]   # z postitions of turbines

        # --- Determine Box extent for FAST>Farm
        FFTS = fastFarmBoxExtent(yBox, zBox, tBox, meanU, hubHeight, D, xWT, yWT, Cmeander=Cmeander, chord_max=chord_max, extent_X=extent_X_high, extent_YZ=extent_Y_high)

        # --- Test values
        #print(FFTS)
        np.testing.assert_almost_equal(FFTS['DT_Low']  , 2.5  , 5)
        np.testing.assert_almost_equal(FFTS['DT_High'] , 0.1  , 5)
        np.testing.assert_almost_equal(FFTS['NX_Low']  , 212  , 5)
        np.testing.assert_almost_equal(FFTS['NY_Low']  , 75   , 5)
        np.testing.assert_almost_equal(FFTS['NZ_Low']  , 56   , 5)
        np.testing.assert_almost_equal(FFTS['X0_Low']  , -154 , 5)
        np.testing.assert_almost_equal(FFTS['Y0_Low']  , -188 , 5)
        np.testing.assert_almost_equal(FFTS['Z0_Low']  , 1    , 5)
        np.testing.assert_almost_equal(FFTS['dX_Low']  , 4.87 , 5)
        np.testing.assert_almost_equal(FFTS['dY_Low']  , 5.0  , 5)
        np.testing.assert_almost_equal(FFTS['dZ_Low']  , 5.0  , 5)
        np.testing.assert_almost_equal(FFTS['NX_High'] , 20   , 5)
        np.testing.assert_almost_equal(FFTS['NY_High'] , 19   , 5)
        np.testing.assert_almost_equal(FFTS['NZ_High'] , 26   , 5)
        np.testing.assert_almost_equal(FFTS['dX_High'] , 4.87 , 5)
        np.testing.assert_almost_equal(FFTS['dY_High'] , 5    , 5)
        np.testing.assert_almost_equal(FFTS['dZ_High'] , 5    , 5)
        np.testing.assert_almost_equal(FFTS['X0_High'] , [-43.83, 219.15], 5)
        np.testing.assert_almost_equal(FFTS['Y0_High'] , [-45   , 5], 5)

        # --- Write Fast Farm file with layout and Low and High res extent
        templateFSTF = os.path.join(MyDir, '../examples/SampleFiles/TestCase.fstf')      # template file used for FastFarm input file, need to exist
        outputFSTF   = os.path.join(MyDir, '../examples/SampleFiles/_TestCase_mod.fstf') # new file that will be written
        writeFastFarm(outputFSTF, templateFSTF, xWT, yWT, zWT, FFTS=FFTS)

if __name__ == '__main__':
    unittest.main()

if __name__ == '__main__':
    unittest.main()

