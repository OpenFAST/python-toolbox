import numpy as np
import pandas as pd
import os
import pyFAST.linearization.linearization as lin
import pyFAST.linearization.LinearModel as lin_mod
import pyFAST.case_generation.case_gen as case_gen
import pyFAST.input_output.fast_output_file as fo
import pyFAST.input_output.rosco_input_file as ri
import sys
from ROSCO_toolbox import utilities as ROSCO_utilities
from ROSCO_toolbox import controller as ROSCO_controller
from ROSCO_toolbox import turbine as ROSCO_turbine
import pyFAST.case_generation.runner as runner
from pCrunch.Analysis import Loads_Analysis
import yaml



import matplotlib.pyplot as plt



if __name__ == '__main__':

    # Load matlab systems, before using Emmanuel's tools
    if True:
        linTurb = lin_mod.LinearTurbineModel('/Users/dzalkind/Tools/matlab-toolbox/Simulations/SaveData/LinearModels/PitTwr.mat', \
            fromMat=True)
    else:
        lin_file_dir = '/Users/dzalkind/Tools/SaveData/TrimTest/LinearTwrPit_Tol1en5'
        linTurb = lin_mod.LinearTurbineModel(lin_file_dir,reduceStates=False)


    # load disturbance from file
    fast_out = fo.FASTOutputFile('/Users/dzalkind/Tools/matlab-toolbox/Simulations/SaveData/072720_183300.out')
    fast_out._read()

    wind_ind = fast_out.info['attribute_names'].index('RtVAvgxh')
    u_h      = fast_out.data[:, wind_ind]

    tt       = fast_out.data[:, fast_out.info['attribute_names'].index('Time')]

    

    # linTurb.solve(tt,u_h)

    # load controller file
    # f = ri.ROSCOInputFile('/Users/dzalkind/Tools/SaveData/TrimTest/LinearTwrPit_Tol1en5/lin_19_DISCON.IN')

    # a = f['LoggingLevel']

    if False:
        fp = ROSCO_utilities.FileProcessing()
        f = fp.read_DISCON('/Users/dzalkind/Tools/SaveData/Float_Test/UM_DLC0_100_DISCON.IN')

        linCont = lin_mod.LinearControlModel([],fromDISCON_IN=True,DISCON_file=f)
    else:
        # Load yaml file 
        parameter_filename = '/Users/dzalkind/Tools/ROSCO_toolbox/Tune_Cases/IEA15MW.yaml'
        inps = yaml.safe_load(open(parameter_filename))
        path_params         = inps['path_params']
        turbine_params      = inps['turbine_params']
        controller_params   = inps['controller_params']

        # Instantiate turbine, controller, and file processing classes
        turbine         = ROSCO_turbine.Turbine(turbine_params)
        controller      = ROSCO_controller.Controller(controller_params)

        # Load turbine data from OpenFAST and rotor performance text file
        turbine.load_from_fast(path_params['FAST_InputFile'],path_params['FAST_directory'],dev_branch=True,rot_source='txt',txt_filename=path_params['rotor_performance_filename'])

        # Tune controller 
        controller.tune_controller(turbine)

        controller.turbine = turbine

        linCont = lin_mod.LinearControlModel(controller)


    Lin_OutList, Lin_OutData, P_cl = linTurb.solve(tt,u_h,Plot=False,open_loop=False,controller=linCont)

    Non_OutList = fast_out.info['attribute_names']
    Non_OutData = fast_out.data



    # comparison plot
    if False:
        comp_channels = ['RtVAvgxh','GenSpeed','TwrBsMyt','PtfmPitch']
        ax = [None] * len(comp_channels)
        plt.figure(2)

        for iPlot in range(0,len(comp_channels)):
            ax[iPlot] = plt.subplot(len(comp_channels),1,iPlot+1)
            try:
                ax[iPlot].plot(tt,Non_OutData[:,Non_OutList.index(comp_channels[iPlot])])
            except:
                print(comp_channels[iPlot] + ' is not in OpenFAST OutList')

            try:
                ax[iPlot].plot(tt,Lin_OutData[:,Lin_OutList.index(comp_channels[iPlot])])
            except:
                print(comp_channels[iPlot] + ' is not in Linearization OutList')
            ax[iPlot].set_ylabel(comp_channels[iPlot])
            ax[iPlot].grid(True)
            if not iPlot == (len(comp_channels) - 1):
                ax[iPlot].set_xticklabels([])

        plt.show()


    # sweep controller.pc_omega and run linearization
    fast_data_lin = []

    if True:
        ww = np.linspace(.05,.45,6)

        for om in ww:
            # Tune controller 
            controller.omega_pc = om
            controller.tune_controller(turbine)

            controller.turbine = turbine

            # update linear controller
            linCont = lin_mod.LinearControlModel(controller)

            # solve 
            Lin_OutList, Lin_OutData, P_cl = linTurb.solve(tt,u_h,Plot=False,open_loop=False,controller=linCont)

            # convert into pCrunch form
            fd_lin  = {}
            fd_lin['TwrBsMyt'] = Lin_OutData[:,Lin_OutList.index('TwrBsMyt')]
            fd_lin['meta']    = {}
            fd_lin['meta']['name'] = 'omega: ' + str(om)
            fast_data_lin.append(fd_lin)



            comp_channels = ['RtVAvgxh','GenSpeed','TwrBsMyt','PtfmPitch']
            ax = [None] * len(comp_channels)
            plt.figure(3)

            for iPlot in range(0,len(comp_channels)):
                ax[iPlot] = plt.subplot(len(comp_channels),1,iPlot+1)

                try:
                    ax[iPlot].plot(tt,Lin_OutData[:,Lin_OutList.index(comp_channels[iPlot])])
                except:
                    print(comp_channels[iPlot] + ' is not in Linearization OutList')
                ax[iPlot].set_ylabel(comp_channels[iPlot])
                ax[iPlot].grid(True)
                if not iPlot == (len(comp_channels) - 1):
                    ax[iPlot].set_xticklabels([])


            

        plt.show()
        

    # Try some post processing using pCrunch
    chan_info = ('TwrBsMyt',4)

    la = Loads_Analysis()
    fDEL = la.get_DEL(fast_data_lin,chan_info)

    fd_non  = {}
    fd_non['TwrBsMyt'] = Non_OutData[:,Non_OutList.index('TwrBsMyt')]
    fd_non['meta']    = {}
    fd_non['meta']['name'] = 'nonlinear'

    fDEL_nl = la.get_DEL([fd_non],chan_info)
    print('here')




