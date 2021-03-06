""" 
Script to post-process several linearization files generated by different OpenFAST simulations.

Each OpenFAST simulation is considered to be at a different operating point (OP).
Typically, this is run for several wind speed/RPM.

A Campbell diagram is plotted, showing the frequencies and damping of each modes for each operating point.

An attempt to identify the turbine modes is done by the script, but a manual sorting is usually needed.
This is done by opening the csv file generated (Campbell_ModesID.csv), and changing the indices. 

The "plot call" at the end of the script can then be repeated with the updated csv file.


"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pyFAST.linearization.linearization as lin

MyDir = os.path.dirname(__file__)

# Script Parameters
BladeLen     = 61.5  # Blade length, used to tune relative modal energy [m] NOTE: not needed if fst files exists
TowerLen     = 87.6  # Tower length, used to tune relative modal energy [m] idem
fstFiles = glob.glob(os.path.join(MyDir,'../../../data/linearization_outputs/*.fst')) # list of fst files where linearization were run, lin file will be looked for
#fstFiles = glob.glob('../../../data/NREL5MW/__5MW_Land_Lin_Trim/*.fst') # list of fst files where linearization were run, lin file will be looked for
fstFiles.sort() # not necessary

# Find lin files, perform MBC, and try to identify modes. A csv file is written with the mode IDs.
OP, Freq, Damp, UnMapped, ModeData, modeID_file = lin.postproCampbell(fstFiles, BladeLen, TowerLen)

# Edit the mode ID file manually to better identify/distribute the modes
print('[TODO] Edit this file manually: ',modeID_file)

# Plot Campbell
fig, axes, figName =  lin.plotCampbellDataFile(modeID_file, 'ws', ylim=None)


if __name__=='__main__':
    plt.show()

if __name__=='__test__':
    # Something weird is happening on github action, order is different, 
    np.testing.assert_almost_equal(Freq['1st_Tower_FA'][:2], [0.324446, 0.331407],3)
    np.testing.assert_almost_equal(Damp['1st_Tower_FA'][:2], [0.00352, 0.06034],4)
    np.testing.assert_almost_equal(OP['WS_[m/s]'], [0, 3],2)
    np.testing.assert_almost_equal(OP['RotSpeed_[rpm]'], [0, 6.972],2)

