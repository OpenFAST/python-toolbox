""" 
Setup a FAST.Farm suite of cases based on input parameters.

The extent of the high res and low res domain are setup according to the guidelines:
    https://openfast.readthedocs.io/en/dev/source/user/fast.farm/ModelGuidance.html

NOTE: If driving FAST.Farm using TurbSim inflow, the resulting boxes are necessary to
      build the final FAST.Farm case and are not provided as part of this repository. 
      If driving FAST.Farm using LES inflow, the VTK boxes are not necessary to exist.

"""

from pyFAST.fastfarm.FASTFarmCaseCreation import FFCaseCreation

def main():

    # -----------------------------------------------------------------------------
    # USER INPUT: Modify these
    #             For the d{t,s}_{high,low}_les paramters, use AMRWindSimulation.py
    # -----------------------------------------------------------------------------

    # ----------- Case absolute path
    path = 'path to case'
    
    # ----------- General hard-coded parameters
    cmax     = 5      # maximum blade chord (m)
    fmax     = 1   # maximum excitation frequency (Hz) 10/6
    Cmeander = 1.9    # Meandering constant (-)

    # ----------- Wind farm
    D = 240
    zhub = 135

    wts  = {
        0 :{'x':-884.53,   'y':-1175.67,    'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        1 :{'x':-1131.09,   'y':16.58,       'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        2 :{'x':-1048.27,   'y':899.75,      'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        3 :{'x':150.84,  'y':-909.41,     'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        4 :{'x':-139.48,  'y':205.46,      'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        5 :{'x':1191.84,  'y':517.81,      'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        6 :{'x':1161.26,  'y':-671.35,     'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        7 :{'x':1076.48,  'y':-167.62,     'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander},
        8 :{'x':1200.86,  'y':1181.68,     'z':0.0,  'D':D,  'zhub':zhub,  'cmax':cmax,  'fmax':fmax,  'Cmeander':Cmeander}
            }

    refTurb_rot = 0
    
    # ----------- Additional variables
    tmax = 1500    # Total windfield time
    tsim= 10         # Simulation time multiple of tmax
    
    nSeeds = 1      # Number of different seeds
    zbot = 10        # Bottom of your domain
    mod_wake = 1    # Wake model. 1: Polar, 2: Curl, 3: Cartesian
    
    # ----------- Desired sweeps
    vhub       = [10]
    shear      = [0.2]
    TIvalue    = [6]
    # inflow_deg = [0]
    # inflow_deg = [0, 90, 180, 270]
    inflow_deg = [0] # if they are not in order from smaller to bigger the folders produced will be wrong ned to change this in the future
    
    # ----------- Turbine parameters
    # Set the yaw of each turbine for wind dir. One row for each wind direction.
    yaw_init = [ [0,0,0,0,0,0,0,0,0] ,
                [0,0,0,0,0,0,0,0,0]]
    
    # ----------- Low- and high-res boxes parameters
    # Should match LES if comparisons are to be made; otherwise, set desired values
    # For an automatic computation of such parameters, omit them from the call to FFCaseCreation
    # High-res boxes settings
    dt_high_les = 0.25                # sampling frequency of high-res files
    ds_high_les = 5               # dx, dy, dz that you want these high-res files at
    extent_high = 3               # high-res box extent in y and x for each turbine, in D.
    # Low-res boxes settings
    dt_low_les  = 4                  # sampling frequency of low-res files
    ds_low_les  = 35               # dx, dy, dz of low-res files
    extent_low  = [1.5, 10,  1.5, 1.5, 0.55]   # extent in xmin, xmax, ymin, ymax, zmax, in D
    # extent_low  = [2, 10,  2,2, 2]   # extent in xmin, xmax, ymin, ymax, zmax, in D
    
    
    # ----------- Execution parameters
    ffbin = 'path to fast farm exe'
    TurbSimbin='path to turbsim exe'
    Mannbin='path to mann box exe'
    
    # ----------- LES parameters. This variable will dictate whether it is a TurbSim-driven or LES-driven case
    # LESpath = '/full/path/to/the/LES/case'
    LESpath = 'MannBox' # set as None if TurbSim-driven is desired
    # LESpath = 'TurbSim' # set as None if TurbSim-driven is desired
    
    floating=1
    
    
    # -----------------------------------------------------------------------------
    # ----------- Template files

    templatePath            = 'path to template'
    
    # Put 'unused' to any input that is not applicable to your case
    # Files should be in templatePath
    EDfilename              = 'Activefloat_ElastoDyn.T'
    SEDfilename             = 'unused'
    HDfilename              = 'Activefloat_HydroDyn.T'
    MDfilename              = 'Activefloat_MoorDyn.T'
    SrvDfilename            = 'Activefloat_ServoDyn.T'
    ADfilename              = 'Activefloat_IEA-15-240-RWT_AeroDyn15.dat'
    ADskfilename            = None
    SubDfilename            = None
    IWfilename              = 'Activefloat_IEA-15-240-RWT_InflowFile.dat'
    BDfilepath              = None
    bladefilename           = 'Activefloat_IEA-15-240-RWT_ElastoDyn_blade.dat'
    towerfilename           = 'Activefloat_IEA-15-240-RWT_ElastoDyn_tower.dat'
    turbfilename            = 'Activefloat.T'
    libdisconfilepath       = 'DISCON.dll'
    controllerInputfilename = 'DISCON.IN'
    coeffTablefilename      = 'Cp_Ct_Cq.IEA15MW.txt'
    FFfilename              = 'Activefloat_Farm.fstf'
    
    # TurbSim setups
    turbsimLowfilepath      = '\SampleFiles\template_Low_InflowXX_SeedY.inp'
    turbsimHighfilepath     = '\SampleFiles\template_HighT1_InflowXX_SeedY.inp'
    
    # SLURM scripts
    slurm_TS_high           = '\SampleFiles/runAllHighBox.sh'
    slurm_TS_low            = '/runAllLowBox.sh'
    slurm_FF_single         = '/runFASTFarm_cond0_case0_seed0.sh'


    # -----------------------------------------------------------------------------
    # END OF USER INPUT
    # -----------------------------------------------------------------------------


    # Initial setup
    case = FFCaseCreation(path, wts, tmax, tsim,zbot, vhub, shear,
                          TIvalue, inflow_deg, dt_high_les, ds_high_les, extent_high,
                          dt_low_les, ds_low_les, extent_low, ffbin,TurbSimbin,Mannbin, mod_wake, LESpath=LESpath,floating=floating,
                          verbose=1)

    case.setTemplateFilename(templatePath, EDfilename, SEDfilename, HDfilename, SrvDfilename, ADfilename,
                             ADskfilename, SubDfilename, IWfilename, BDfilepath, bladefilename, towerfilename,
                             turbfilename, libdisconfilepath, controllerInputfilename, coeffTablefilename,
                             turbsimLowfilepath, turbsimHighfilepath, FFfilename,MDfilename)

    # Get domain paramters
    case.getDomainParameters()

    # Organize file structure
    if LESpath =='MannBox':
        case.Create_Mannbox(33.6, 3.9)

    if LESpath =='TurbSim':
        case.TS_low_setup()
        case.TS_low_bat_file()
        case.TS_low_bat_execute()


        case.TS_high_setup()
        case.TS_high_bat_file()
        case.TS_high_bat_execute()

    # Final setup
    case.FF_setup()
    case.copyTurbineFilesForEachCase()
    case.FF_bat_file()
    # case.FF_execute_bat()



if __name__ == '__main__':
    # This example cannot be fully run.
    pass
