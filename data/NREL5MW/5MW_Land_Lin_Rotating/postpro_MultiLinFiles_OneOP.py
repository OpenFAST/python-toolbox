""" 
Script to postprocess linearization files from OpenFAST for one operating point.

Adapted from:
     https://github.com/OpenFAST/python-toolbox/blob/dev/openfast_toolbox/linearization/examples/ex2a_MultiLinFiles_OneOP.py

"""
import os
import glob
import numpy as np
import openfast_toolbox
import openfast_toolbox.linearization as lin

## Script Parameters
fstFile = './Main.fst' # Main .fst file, .lin files will be sought for with same basename

# Get Campbell Diagram Data for one Operating Point (CDDOP) given an OpenFAST input file
# Perform MBC transformation based on all lin files found next to .fst file
CDDOP, MBCOP = lin.getCDDOP(fstFile, writeModes=True, verbose=True) 

# Outputs to screen
Freq,Damp = lin.printCDDOP(CDDOP, nModesMax=10, nCharMaxDesc=50)
