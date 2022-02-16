""" 
Setup a FAST.Farm fstf input file using a TurbSim box:
  - Low Res grid extent is based on the TurbSim box data
  - High Res grid extent and location is set based 



NOTE: to run this script you first need to generate the TurbSim box from 
      the file SampleFiles/TestCase.inp, as follows:
         turbsim.exe TestCase.inp


"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Local packages
from pyFAST.fastfarm.FFarmCaseCreation import FFarmCaseCreation, WriteFFarmFile
from pyFAST.fastfarm.fastfarm import plotFastFarmSetup

# --- Parameters required for this script
FSTPath       = '../turbineModel/Test18'        # Path to FST files,  do not include .fst extention
BTSFilename   = 'SampleFiles/TestCase.bts'       # TurbSim Box to be used in FAST.Farm simulation, need to exist.
OldFSTFFile   = 'SampleFiles/TestCase.fstf'     # template file used for FastFarm input file, need to exist
NewTSTFFile   = 'SampleFiles/_TestCase_mod.fstf' # new file that will be written
D             = 77.0                            # Turbine diameter (m)
HubHt         = 78.045                          # Hub Height (m)
high_extent_X = 1.2                             # x-extent of high res box in diamter around turbine location
high_extent_Y = 1.2                             # y-extent of high res box in diamter around turbine location
xlocs = [0.0, 265.643]  # x positions of turbines
ylocs = [0.0, 50.0   ]  # y postitions of turbines
# xlocs = [0.0      , 265.643, 653.506,  871.276,  653.901]  # x positions of turbines
# ylocs = [-377.410 , 0.0    ,  -5.378,  28.494 ,  342.903]  # y postitions of turbines

# --- "Optional" inputs
cmax     = 5   # maximum blade chord (m). Turbine specific.
Cmeander = 1.9 # Meandering constant (-)

# --- Setup FASTFarm Case creation class
Case = FFarmCaseCreation(D, HubHt, x=xlocs, y=ylocs, cmax = cmax)

# --- Setup low res and high res based on BTS file
Case.setupFromBTS(BTSFilename, high_extent_X, high_extent_Y, Cmeander=Cmeander)

# --- Write FFarm Input File
Case.writeFFarmFile(OldFSTFFile, NewTSTFFile, FSTPath, NewFile=False)

# --- Visualize low&high extent and turbine positions
fig, ax = Case.plotSetup()    # Visualize based in "Case" data
plotFastFarmSetup(NewTSTFFile) # Visualize based on what is in the new fstf file

plt.show()
