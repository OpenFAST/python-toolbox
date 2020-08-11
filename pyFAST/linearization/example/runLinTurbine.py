import numpy as np
import pandas as pd
import os
import pyFAST.linearization.linearization as lin
import pyFAST.linearization.LinearModel as lin_mod
import pyFAST.case_generation.case_gen as case_gen
import pyFAST.input_output.fast_output_file as fo
import pyFAST.case_generation.runner as runner
from scipy.io import loadmat

import matplotlib.pyplot as plt


if __name__ == '__main__':

    # Load matlab systems, before using Emmanuel's tools
    m = loadmat('/Users/dzalkind/Tools/matlab-toolbox/Simulations/SaveData/LinearModels/PitTwr.mat')

    linTurb = lin_mod.LinearTurbineModel(fromMat=True,matDict = m)

    # load disturbance from file
    fast_out = fo.FASTOutputFile('/Users/dzalkind/Tools/matlab-toolbox/Simulations/SaveData/072720_183300.out')
    fast_out._read()

    wind_ind = fast_out.info['attribute_names'].index('RtVAvgxh')
    u_h      = fast_out.data[:,wind_ind]

    tt       = fast_out.data[:,fast_out.info['attribute_names'].index('Time')]

    

    linTurb.solve(tt,u_h)
    
    print('here')



