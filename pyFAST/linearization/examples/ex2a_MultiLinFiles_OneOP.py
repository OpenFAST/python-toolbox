""" 
Script to postprocess linearization multiple lin files from OpenFAST.
 - the multiple lin files are obtained for the same periodic operating point
 - the different lin files are at differnt azimuthal positions
 - MBC3 is performed
"""
import os
import glob
import numpy as np
import pyFAST.linearization as lin

scriptDir = os.path.dirname(__file__)

## Script Parameters
simDir      = os.path.join(scriptDir,'../../../data/NREL5MW/5MW_Land_Lin_Rotating/') # Simulation directory
fstFile     = os.path.join(simDir,'./Main.fst') # fstFile, lin files will be assumed to be basename.i.lin

# Get Campbell Diagram Data for one Operating Point (CDDOP) given an OpenFAST (OF) input file
# Perform MBC transformation based on all lin files found next to .fst file
CDDOP, MBCOP = lin.getCDDOP(fstFile) # alternatively provide a list of lin files

# Outputs to screen
Freq,Damp = lin.printCDDOP(CDDOP, nModesMax=10, nCharMaxDesc=50)


if __name__=='__main__':
    pass

if __name__=='__test__':
    np.testing.assert_almost_equal(Freq[:4],     [0.588,  0.722 , 0.842, 0.937],3)
    np.testing.assert_almost_equal(Damp[:4]*100, [63.106, 52.529, 44.01, 1.634],3)
