import numpy as np
import pandas as pd
import os
import pyFAST.linearization.linearization as lin
import pyFAST.case_generation.case_gen as case_gen
import pyFAST.case_generation.runner as runner

import matplotlib.pyplot as plt


if __name__=='__main__':

    # --- Parameters to generate linearization input files
    templateFstFile     = 'C:/Work/FAST/matlab-toolbox/_ExampleData/5MW_Land_Lin_Templates/Main_5MW_Land_Lin.fst'  # Main file, used as a template
    simulationFolder    = 'C:/Work/FAST/matlab-toolbox/_ExampleData/5MW_Land_Lin/'  # Output folder for input files and linearization (will be created)
    operatingPointsFile = 'C:/Work/FAST/matlab-toolbox/Campbell/example/LinearizationPoints_NoServo.csv'
    tStart           = 500  # Time after which linearization is done (need to reach periodic steady state)
    nPerPeriod       = 12   # Number of linearization per revolution
    
    # --- Parameters to run OpenFAST
    runFast = True
    fastExe = 'C:/Work/FAST/matlab-toolbox/_ExampleData/openfast2.3_x64s.exe' # Path to a FAST exe (and dll) 

    # --- Parameters for MBC postpro
    runMBC = True
    toolboxDir = 'C:/Work/FAST/matlab-toolbox/'           # path to matlab-toolbox
    matlabExe  = 'C:/Bin/Octave-4.4.1/bin/octave-cli.exe' # path the matlab or octave exe

    # --- Step 1: Write OpenFAST inputs files for each operating points 
    baseDict={'DT':0.01} # Example of how inputs can be overriden (see case_gen.py templateReplace)
    FSTfilenames= lin.writeLinearizationFiles(templateFstFile, simulationFolder, operatingPointsFile, nPerPeriod=nPerPeriod, baseDict=baseDict, tStart=tStart)

    # Create a batch script (optional)
    runner.writeBatch(os.path.join(simulationFolder,'_RUN_ALL.bat'), FSTfilenames, fastExe=fastExe)

    # --- Step 2: run OpenFAST 
    if runFast:
        runner.run_fastfiles(FSTfilenames, fastExe=fastExe, parallel=True, showOutputs=False, nCores=4)

    # --- Step 3: Run MBC, identify modes and generate XLS or CSV file
    if runMBC:
        lin.postproLinearization(simulationFolder, operatingPointsFile, toolboxDir, matlabExe)

    # --- Step 4: Campbell diagram plot
    csvFile = os.path.join(simulationFolder, 'Campbell_ModesID.csv') # <<< TODO Change me if manual identification is done
    fig, axes, figName = lin.plotCampbellDataFile(csvFile, ws_or_rpm='ws', ylim=[0,4])
    #  fig.savefig(figName+'.png')
    plt.show()


