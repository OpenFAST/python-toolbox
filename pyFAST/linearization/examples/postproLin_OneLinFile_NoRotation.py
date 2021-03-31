""" 
Standalone script to post-process one linearization file from OpenFAST.
NOTE: this should not be used if the rotor is turning.
This script would typically be used for a standstill analysis (no rotation),
or to compute modes of when only isolated degrees of freedom are turned on (and no rotation).

NOTE: should match the script found in the matlab-toolbox
"""
import os
import numpy as np
import pyFAST.linearization.mbc.mbc3 as mbc # TODO will be moved in the future

MyDir = os.path.dirname(__file__)

## Script Parameters
BladeLen     = 40.04                # Blade length, used to tune relative modal energy [m]
TowerLen     = 55.59                # Tower length, used to tune relative modal energy [m]
lin_file     = os.path.join(MyDir,'../../../data/example_files/Standstill.1.lin') # Linearization file
#lin_file     = '../../../data/example_files/Standstill_old.1.lin' # Linearization file
nModesMax    = 10                   # Maximum number of modes to be shown
nCharMaxDesc = 50                   # Maximum number of characters for description written to screen

## Derived parameters
lin_files = np.array([lin_file])

# Performing MBC (NOTE: not stricly necessary without rotation)
mbc_data, matData, FAST_linData = mbc.fx_mbc3(lin_files, verbose=False)
CD = mbc.campbell_diagram_data(mbc_data,BladeLen,TowerLen)

# Outputs to screen
nModesMax = np.min([len(CD['Modes']),nModesMax])
Freq = np.array([CD['Modes'][i]['NaturalFreq_Hz'] for i in np.arange(nModesMax)])
Damp = np.array([CD['Modes'][i]['DampingRatio']   for i in np.arange(nModesMax)])
print('Mode, NatFreq_[Hz], Damp_Ratio_[-], LogDec._[%], Mode_content_[-]')
for i in np.arange(nModesMax):
    Mode = CD['Modes'][i]
    # Extracting description the best we can
    Desc = mbc.extractShortModeDescription(Mode)
    print('{:3d} ,{:12.3f}, {:8.5f}       , {:7.4f},  {:s}'.format(i+1,Mode['NaturalFreq_Hz'],Mode['DampingRatio'],Mode['DampingRatio']*100*2*np.pi, Desc[:min(nCharMaxDesc,len(Desc))]))

if __name__=='__main__':
    pass

if __name__=='__test__':
    np.testing.assert_almost_equal(Freq[:3]  ,[0.427, 0.450, 0.669], 3)
    np.testing.assert_almost_equal(Damp[:3]*np.pi*2*100,[1.9505,2.1309,5.0649], 4)
