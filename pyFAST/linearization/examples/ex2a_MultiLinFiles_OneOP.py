""" 
Script to postprocess linearization multiple lin files from OpenFAST.
 - the multiple lin files are obtained for the same periodic operating point
 - the different lin files are at differnt azimuthal positions
 - MBC3 is performed


"""
import os
import glob
import numpy as np
import pyFAST.linearization.mbc.mbc3 as mbc # TODO will be moved in the future
import pyFAST.linearization as lin

scriptDir = os.path.dirname(__file__)

## Script Parameters
simDir       = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Land_Lin_Rotating/') # Simulation directory
basename     = os.path.join(simDir,'./Main') # Basename, lin files will be assumed to be base_name.i.lin
nLin         = 36       # Number of linearization file for this operating point
BladeLen     = 63       # Blade length, used to tune relative modal energy [m]
TowerLen     = 87.6     # Tower length, used to tune relative modal energy [m]
nModesMax    = 10       # Maximum number of modes to be shown
nCharMaxDesc = 50       # Maximum number of characters for description written to screen


# List of lin files for this operating point
#lin_files = np.array([basename+'.{}.lin'.format(i+1) for i in range(nLin)])
lin_files = glob.glob(basename+'.*.lin')

# TODO Simplify interface
mbc_data, matData, FAST_linData = mbc.fx_mbc3(lin_files, verbose=False)
# mbc_data, matData, FAST_linData = lin.MBC_OF(basename+'.fst')
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
    np.testing.assert_almost_equal(Freq[:4],     [0.588,  0.722 , 0.842, 0.937],3)
    np.testing.assert_almost_equal(Damp[:4]*100, [63.106, 52.529, 44.01, 1.634],3)
