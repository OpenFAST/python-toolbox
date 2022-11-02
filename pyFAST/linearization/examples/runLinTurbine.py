""" Example to run a set of OpenFAST simulations with linearizations.
NOTE: for a more streamlined process to generate a Cmpbell diagram look at the example runCampbell.py

This script uses a reference directory which contains a reference input file (templateFstFile)
1) The reference directory is copied to a working directory (`simulationFolder`).
2) All the fast input files are generated in this directory based on a list of dictionaries (`PARAMS`).
For each dictionary in this list:
   - The keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `TMax`.
     These should correspond to the variables used in the FAST inputs files.
   - The values are the values corresponding to this parameter
For instance:
     PARAMS[0]['DT']                    = 0.01
     PARAMS[0]['EDFile|RotSpeed']       = 5
     PARAMS[0]['InflowFile|HWindSpeed'] = 10

3) The simulations are run, successively distributed on `nCores` CPUs.

4) The linearization file are postprocessed using MBC and the frequencies written to screen

"""
import numpy as np
import pandas as pd
import os
import pyFAST.linearization.linearization as lin
import pyFAST.case_generation.case_gen as case_gen
import pyFAST.case_generation.runner as runner
import scipy as sp

import matplotlib.pyplot as plt

def run_linearization():
    # --- Parameters for this script
    simulationFolder    = '../../../data/NREL5MW/_5MW_Land_Lin_Trim2/'  # Output folder for input files and linearization (will be created)
    templateFstFile     = '../../../data/NREL5MW/5MW_Land_Lin_Templates/Main_5MW_Land_Lin.fst'  # Main file, used as a template
    fastExe             = '../../../data/openfast2.5_x64.exe' # Path to a FAST exe (and dll) 

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS = [14,16,18]
    BaseDict = {'TMax': 600}
    # Set some default options 
    BaseDict = case_gen.paramsLinearTrim(BaseDict)   # Run linear trim case
    print(BaseDict)
    PARAMS=[]
    for i,wsp in enumerate(WS): 
        p=BaseDict.copy()
        p['CompServo']              = 1
        p['InflowFile|HWindSpeed']  = wsp
        p['InflowFile|WindType']    = 1 # Setting steady wind
        
        # Set DOFs
        p['EDFile|GenDOF']          = 'True'
        p['EDFile|FlapDOF1']        = 'True'
        p['EDFile|TwFADOF1']        = 'True'

        # NREL-5MW rated generator speed, torque, control params
        p['ServoFile|VS_RtGnSp']    = 1173 * .9
        p['ServoFile|VS_RtTq']      = 47402
        p['ServoFile|VS_Rgn2K']     = 0.0226
        p['ServoFile|VS_SlPc']      = 10.

        # Trim solution will converge to this rotor speed
        p['EDFile|RotSpeed']        = 12.1

        # Set number of linearizations
        p['CalcSteady']             = True
        p['NLinTimes']              = 12
        p['OutFmt']                 = '"ES20.12E3"'  # Important for decent resolution
        p['LinInputs']              = 0
        p['LinOutputs']             = 0
        p['TrimCase']               = 3
        p['TrimGain']               = 0.001

        p['__name__']='{:03d}_ws{:04.1f}'.format(i,p['InflowFile|HWindSpeed'])
        PARAMS.append(p)
        i=i+1
    # --- Generating all files in a output directory
    refDir    = os.path.dirname(templateFstFile)
    main_file = os.path.basename(templateFstFile)
    fastfiles=case_gen.templateReplace(PARAMS, refDir, outputDir=simulationFolder, removeRefSubFiles=True, main_file=main_file)
    print(fastfiles)

    # --- Creating a batch script just in case
    runner.writeBatch(os.path.join(simulationFolder,'_RUN_ALL.bat'),fastfiles,fastExe=fastExe)
    # --- Running the simulations
    runner.run_fastfiles(fastfiles, fastExe=fastExe, parallel=True, showOutputs=True, nCores=1)

    # --- Simple Postprocessing
    # (averaging each signal over the last period for each simulation)
    #outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
    # avg_results = postpro.averagePostPro(outFiles, avgMethod='periods', avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    # avg_results.drop('Time_[s]',axis=1, inplace=True)

    return fastfiles



if __name__ == '__main__':

    # 1. Run linearizations
    fastfiles = run_linearization()
    # 2. Do MBC
    #MBC = lin.run_pyMBC(fastfiles)
    MBC = lin.getMBCOPs(fastfiles)

    # Check natural frequencies against matlab
    for mbc in MBC:
        eigs = sp.linalg.eig(mbc['AvgA'])[0]
        nat_freq_hz = (np.abs(eigs)/2/np.pi)
        nat_freq_hz.sort()
        print('Natural Freqs. (hz): {}'.format(nat_freq_hz))

if __name__=='__test__':
    pass # this example needs an openfast binary
