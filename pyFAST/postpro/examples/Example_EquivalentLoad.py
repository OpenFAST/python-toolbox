""" 
- Open and OpenFAST binary file
- Convert it to a pandas dataframe
- Compute damage equivalent load for a given Wohler exponent
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from pyFAST.input_output import FASTOutputFile
from pyFAST.postpro import equivalent_load

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)

# Read an openFAST binary
fastoutFilename = os.path.join(scriptDir, '../../../data/example_files/fastout_allnodes.outb')
df = FASTOutputFile(fastoutFilename).toDataFrame()


# Compute equivalent load for one signal and Wohler slope
T = df['Time_[s]'].values[-1] # number of 1Hz load cycles (time series length in second)
m = 10 # Wohler slope
Leq = equivalent_load(df['RootMyc1_[kN-m]'], Teq=T,  m=1)
print('Leq ',Leq)
# Leq = equivalent_load(df['RootMyc1_[kN-m]'], Teq=T,  m=1, method='fatpack') # requires package fatpack


if __name__ == '__main__':
    plt.show()
if __name__ == '__test__':
    np.testing.assert_almost_equal(Leq , 284.30398, 3)
