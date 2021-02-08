import numpy as np
import pandas as pd
import os
import pyFAST.linearization.linearization as lin
# import pyFAST.linearization.LinearModel as lin_mod
import pyFAST.case_generation.case_gen as case_gen
import pyFAST.input_output.fast_output_file as fo
import sys
import pyFAST.case_generation.runner as runner
# from pCrunch.Analysis import Loads_Analysis
import yaml
import scipy as sp



import matplotlib.pyplot as plt

def run_linearization():
    """ Example to run a set of OpenFAST simulations (linearizations)

    This script uses a reference directory (`ref_dir`) which contains a reference input file (.fst)
    1) The reference directory is copied to a working directory (`out_dir`).
    2) All the fast input files are generated in this directory based on a list of dictionaries (`PARAMS`).
    For each dictionary in this list:
       - The keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `FAST|TMax`.
         These should correspond to the variables used in the FAST inputs files.
       - The values are the values corresponding to this parameter
    For instance:
         PARAMS[0]['DT']                    = 0.01
         PARAMS[0]['EDFile|RotSpeed']       = 5
         PARAMS[0]['InflowFile|HWindSpeed'] = 10

    3) The simulations are run, successively distributed on `nCores` CPUs.
    4) The output files are read, and averaged based on a method (e.g. average over a set of periods,
        see averagePostPro in postpro for the different averaging methods).
       A pandas DataFrame is returned

    """
    # --- Parameters for this script
    of_dir           = '/Users/dzalkind/Tools/openfast-dev'  # openfast dir, using dev branch as of Jan-14
    this_dir         = os.path.dirname(__file__)
    ref_dir          = '/Users/dzalkind/Tools/ROSCO_toolbox/Test_Cases/NREL-5MW'   # Folder where the fast input files are located (will be copied)
    out_dir          = os.path.join(this_dir,'NREL-5MW_Linear/')     # Output folder (will be created)
    main_file        = 'NREL-5MW.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = os.path.join(of_dir,'install/bin/openfast') # Location of a FAST exe (and dll)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    WS = [14,16,18]
    BaseDict = {'TMax': 600}
    BaseDict = case_gen.paramsLinearTrim(BaseDict)   # Run linear trim case

    PARAMS=[]
    for i,wsp in enumerate(WS): 
        p=BaseDict.copy()
        p['InflowFile|HWindSpeed']  = wsp
        p['InflowFile|WindType']    = 1 # Setting steady wind
        
        # Set DOFs
        p['EDFile|GenDOF']          = 'True'
        p['EDFile|FlapDOF1']         = 'True'
        p['EDFile|TwFADOF1']         = 'True'

        # NREL-5MW rated generator speed, torque, control params
        p['ServoFile|VS_RtGnSp']    = 1173 * .9
        p['ServoFile|VS_RtTq']      = 47402
        p['ServoFile|VS_Rgn2K']     = 0.0226
        p['ServoFile|VS_SlPc']      = 10.

        # Trim solution will converge to this rotor speed
        p['EDFile|RotSpeed']        = 12.1

        # Set number of linearizations
        p['NLinTimes']              = 12

        p['__name__']='{:03d}_ws{:04.1f}'.format(i,p['InflowFile|HWindSpeed'])
        PARAMS.append(p)
        i=i+1
    # --- Generating all files in a output directory
    fastfiles=case_gen.templateReplace(PARAMS,ref_dir,outputDir=out_dir,removeRefSubFiles=True,main_file=main_file, oneSimPerDir=False)
    print(fastfiles)

    # --- Creating a batch script just in case
    runner.writeBatch(os.path.join(out_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    runner.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,showOutputs=True,nCores=4)

    # --- Simple Postprocessing
    # (averaging each signal over the last period for each simulation)
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
    # avg_results = postpro.averagePostPro(outFiles, avgMethod='periods', avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    # avg_results.drop('Time_[s]',axis=1, inplace=True)

    return outFiles



if __name__ == '__main__':

    # 1. Run linearizations
    outfiles = run_linearization()


    # 2. Do MBC
    MBC = lin.run_pyMBC(outfiles)

    # Check natural frequencies against matlab
    for mbc in MBC:
        eigs = sp.linalg.eig(mbc['AvgA'])[0]
        nat_freq_hz = (np.abs(eigs)/2/np.pi)
        nat_freq_hz.sort()
        print('Natural Freqs. (hz): {}'.format(nat_freq_hz))

