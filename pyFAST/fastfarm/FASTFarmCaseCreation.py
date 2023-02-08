import pandas as pd
import numpy as np
import os, sys, shutil
import subprocess
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from pyFAST.input_output import FASTInputFile, FASTOutputFile, TurbSimFile, VTKFile
from pyFAST.fastfarm import writeFastFarm, fastFarmTurbSimExtent, plotFastFarmSetup
from pyFAST.fastfarm.TurbSimCaseCreation import TSCaseCreation, writeTimeSeriesFile

def cosd(t): return np.cos(np.deg2rad(t))
def sind(t): return np.sin(np.deg2rad(t))


class FFCaseCreation:


    def __init__(self, path, wts, cmax, fmax, Cmeander, tmax, zbot, vhub, shear, TIvalue, inflow_deg, dt_high_les=None, ds_high_les=None, extent_high=None, dt_low_les=None, ds_low_les=None, extent_low=None, ffbin=None, yaw_init=None, ADmodel=None, EDmodel=None, nSeeds=6, LESpath=None, sweepWakeSteering=False, sweepYawMisalignment=False, seedValues=None, refTurb_rot=0, verbose=0):
        '''
        
        ffbin: str
            Full path of the FAST.Farm binary to be executed
        refTurb_rot: int
            Index of reference turbine which the rotation of the farm will occur. Default is 0, the first one.
        '''

        self.path        = path
        self.wts         = wts
        self.cmax        = cmax
        self.fmax        = fmax
        self.Cmeander    = Cmeander
        self.tmax        = tmax
        self.zbot        = zbot
        self.vhub        = vhub
        self.shear       = shear
        self.TIvalue     = TIvalue
        self.inflow_deg  = inflow_deg
        self.dt_high_les = dt_high_les
        self.ds_high_les = ds_high_les
        self.dt_low_les  = dt_low_les
        self.ds_low_les  = ds_low_les
        self.extent_low  = extent_low
        self.extent_high = extent_high
        self.ffbin       = ffbin
        self.yaw_init    = yaw_init
        self.ADmodel     = ADmodel
        self.EDmodel     = EDmodel
        self.nSeeds      = nSeeds
        self.LESpath     = LESpath
        self.sweepWS     = sweepWakeSteering
        self.sweepYM     = sweepYawMisalignment
        self.seedValues  = seedValues
        self.refTurb_rot = refTurb_rot
        self.verbose     = verbose
                                        
                                        
        if self.verbose>0: print(f'Checking inputs...', end='\r')
        self._checkInputs()       
        if self.verbose>0: print(f'Checking inputs... Done.')


        if self.verbose>0: print(f'Setting rotor parameters...', end='\r')
        self._setRotorParameters()
        if self.verbose>0: print(f'Setting rotor parameters... Done.')
                                        
                                        
        if self.verbose>0: print(f'Creating auxiliary arrays for all conditions and cases...', end='\r')
        self.createAuxArrays()          
        if self.verbose>0: print(f'Creating auxiliary arrays for all conditions and cases... Done.')
                                        

        if self.verbose>0: print(f'Determining box parameters...', end='\r')
        self._determineBoxParamters()       
        if self.verbose>0: print(f'Determining box paramters... Done.')


        if self.verbose>0: print(f'Creating directory structure and copying files...', end='\r')
        self._create_dir_structure()
        if self.verbose>0: print(f'Creating directory structure and copying files... Done.')


        #if self.verbose: print(f'Setting up FAST.Farm based on {self.inflowStr} inflow...', end='\r')
        #if self.verbose: print(f'Setting up FAST.Farm based on {self.inflowStr} inflow... Done.'


    def _checkInputs(self):
  
        # Create case path is doesn't exist
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        # Check the wind turbine dict
        if not  isinstance(self.wts,dict):
            raise ValueError (f'`wts` needs to be a dictionary with the following entries for each turbine: x, y, z, D, zhub.')
        self.nTurbines = len(self.wts)
        self.D         = self.wts[0]['D']
        self.zhub      = self.wts[0]['zhub']

        # Check values of each turbine
        for t in range(self.nTurbines):
            t_x    = self.wts[t]['x']
            t_y    = self.wts[t]['y']
            t_z    = self.wts[t]['z']
            t_D    = self.wts[t]['D']
            t_zhub = self.wts[t]['zhub']
            if t_D != self.D:
                raise ValueError(f'Different turbines are not currently supported. Turbine {t+1} has a different diamenter.')
            if t_zhub != self.zhub:
                raise ValueError(f'Different turbines are not currently supported. Turbine {t+1} has a different hub height.')
            if not isinstance(t_x,(float,int)):
                raise ValueError (f'The `x` value for the turbine {t+1} should be an integer or float. Received {t_x}.')
            if not isinstance(t_y,(float,int)):
                raise ValueError (f'The `y` value for the turbine {t+1} should be an integer or float. Received {t_y}.')
            if not isinstance(t_z,(float,int)):
                raise ValueError (f'The `z` value for the turbine {t+1} should be an integer or float. Received {t_z}.')
            if not isinstance(t_D,(float,int)):
                raise ValueError (f'The `D` value for the turbine {t+1} should be an integer or float. Received {t_D}.')
            if not isinstance(t_zhub,(float,int)):
                raise ValueError (f'The `zhub` value for the turbine {t+1} should be an integer or float. Received {t_zhub}.')
  
        # Check general variables
        if self.cmax     <= 0: raise ValueError('cmax cannot be negative')
        if self.fmax     <= 0: raise ValueError('fmax cannot be negative')
        if self.Cmeander <= 0: raise ValueError('Cmeander cannot be negative')
        if self.tmax     <= 0: raise ValueError('A positive tmax should be requested')
        if self.zbot     <= 0: raise ValueError('zbot should be greater than 0 (recommended 1)')
  
        # Ensure quantities are list
        self.vhub       = [self.vhub]       if isinstance(self.vhub,(float,int))       else self.vhub
        self.shear      = [self.shear]      if isinstance(self.shear,(float,int))      else self.shear
        self.TIvalue    = [self.TIvalue]    if isinstance(self.TIvalue,(float,int))    else self.TIvalue
        self.inflow_deg = [self.inflow_deg] if isinstance(self.inflow_deg,(float,int)) else self.inflow_deg

        # Fill turbine parameters arrays if not given
        if self.yaw_init is None:
            yaw = np.ones((1,self.nTurbines))*0
            self.yaw_init = np.repeat(yaw, len(self.inflow_deg), axis=0)      

        # Check TI values if given in percent
        for t in self.TIvalue:
            if t<0:  raise ValueError(f'TI cannot be negative. Received {t}.')
            if t<1:  raise ValueError(f'TI should be given in percentage (e.g. "10" for a 10% TI). Received {t}.')

        # Check the ds and dt for the high- and low-res boxes
        if None not in (self.dt_high_les, self.dt_low_les, self.ds_high_les, self.ds_low_les, self.extent_high, self.extent_low):
            # All values for ds, dt, and extent were given
            if not (np.array(self.extent_low)>=0).all():
                raise ValueError(f'The array for low-res box extents should be given with positive values')
            if self.dt_low_les%(self.dt_high_les-1e-15) > 1e-14:
                raise ValueError(f'The temporal resolution dT_Low should be a multiple of dT_High')
            if self.dt_low_les < self.dt_high_les:
                raise ValueError(f'The temporal resolution dT_High should not be greater than dT_Low on the LES side')
            if self.ds_low_les < self.ds_high_les:
                raise ValueError(f'The grid resolution dS_High should not be greater than dS_Low on the LES side')
            if not isinstance(self.extent_high, (float,int)):
                raise ValueError(f'The extent_high should be a scalar')
            if self.extent_high<=0:
                raise ValueError(f'The extent of high boxes should be positive')
        else:
            # At least one of the ds, dt, extent are None. If that is the case, all of them should be None
            if None in (self.dt_high_les, self.dt_low_les, self.ds_high_les, self.ds_low_les, self.extent_high, self.extent_low):
                raise ValueError (f'The dt, ds, and extent for high and low res boxes need to either be fully' \
                                  f'given or computed. Some were given, but some were not.')


        # Check the FAST.Farm binary
        if self.ffbin is None:
            import shutil
            self.ffbin = shutil.which('FAST.Farm')
            if verbose>1:
                print(f'FAST.Farm binary not given. Using the following from your $PATH: {self.ffbin}.')
        if not os.path.isfile(self.ffbin):
            raise ValueError (f'The FAST.Farm binary given does not appear to exist')


        # Check turbine conditions arrays for consistency
        if len(self.inflow_deg) != len(self.yaw_init):
            raise ValueError(f'One row for each inflow angle should be given in yaw_init. '\
                             f'Currently {len(self.inflow_deg)} inflow angle(s) and {len(self.yaw_init)} yaw entrie(s)')
        if self.ADmodel is None:
            self.ADmodel = np.tile(['ADyn'],(1,self.nTurbines))
        if self.EDmodel is None:
            self.EDmodel = np.tile(['FED'],(1,self.nTurbines))
        if np.shape(self.ADmodel) != np.shape(self.EDmodel):
            raise ValueError('Every case should have the aerodynamic and elastic model selected. The number of cases '\
                             '(lines) in `ADmodel` and `EDmodel` should be the same')
        if self.nTurbines != np.shape(self.ADmodel)[1]:
            raise ValueError(f'The number of turbines in wts ({len(self.wts)}) should match the number of turbines '\
                             f'in the ADmodel and EDmodel arrays ({np.shape(self.ADmodel)[1]})')
  
        # Check on seed parameters
        if not isinstance(self.nSeeds,int):
            raise ValueError(f'An integer number of seeds should be requested. Got {self.nSeeds}.')
        if self.seedValues is None:
            self.seedValues = [2318573, 122299, 123456, 389432, -432443, 9849898]
        if len(self.seedValues) != self.nSeeds:
            raise ValueError(f'Number of seeds is {self.nSeeds} but {len(self.seedValues)} seed values were given. '
                             f'Adjust the seedValues array accordingly')
  
        # Check LES parameters
        if self.LESpath is None:
            self.inflowStr = 'TurbSim'
        else:
            if not os.path.isdir(self.LESpath):
                raise ValueError (f'The path {self.LESpath} does not exist')
            self.inflowStr = 'LES'

        # Check the reference turbine for rotation
        if self.refTurb_rot >= self.nTurbines:
            raise ValueError(f'The index for the reference turbine for the farm to be rotated around is greater than the number of turbines')

        # Set aux variable
        self.templateFilesCreatedBool = False
        self.TSlowBoxFilesCreatedBool = False


    def _create_dir_structure(self):     
        # Create directory structure CondXX_*/CaseYY_*/Seed_Z/TurbSim; and CondXX_*/Seed_Y
        # Also saves the dir structure on array to write on SLURM script in the future.
        condDirList = []
        caseDirList_ = []
        for cond in range(self.nConditions):
            # Recover information about current condition for directory naming purposes
            Vhub_    = self.allCond['vhub'      ].isel(cond=cond).values
            shear_   = self.allCond['shear'     ].isel(cond=cond).values
            tivalue_ = self.allCond['TIvalue'   ].isel(cond=cond).values
            
            # Set current path name string
            condStr = f'Cond{cond:02d}_v{Vhub_:04.1f}_PL{shear_}_TI{tivalue_}'
            condDirList.append(condStr)
            condPath = os.path.join(self.path, condStr)
            
            for case in range(self.nCases):
                # Recover information about current case for directory naming purposes
                inflow_deg_   = self.allCases['inflow_deg'     ].sel(case=case).values
                wakeSteering_ = self.allCases['wakeSteering'   ].sel(case=case).values
                misalignment_ = self.allCases['misalignment'   ].sel(case=case).values
                nADyn_        = self.allCases['nFullAeroDyn'   ].sel(case=case).values
                nFED_         = self.allCases['nFulllElastoDyn'].sel(case=case).values
                yawCase_      = self.allCases['yawCase'        ].sel(case=case).values
            
                # Set current path name string. The case is of the following form: Case00_wdirp10_WSfalse_YMfalse_12fED_12ADyn
                caseStr = f"Case{case:02d}_wdir{f'{int(inflow_deg_):+03d}'.replace('+','p').replace('-','m')}"\
                          f"_WS{str(wakeSteering_).lower()}_YM{str(misalignment_).lower()}"\
                          f"_{nFED_}fED_{nADyn_}ADyn"
                # If sweeping on yaw, then add yaw case to dir name
                if len(np.unique(self.allCases.yawCase)) > 1:
                    caseStr = caseStr + f"_yawCase{yawCase_}"
                
                caseDirList_.append(caseStr)
                casePath = os.path.join(condPath, caseStr)
                if not os.path.exists(casePath):  os.makedirs(casePath)
                
                for seed in range(self.nSeeds):
                    seedPath = os.path.join(casePath, f'Seed_{seed}')
                    if not os.path.exists(seedPath):  os.makedirs(seedPath)
                    turbsimPath = os.path.join(seedPath, 'TurbSim')
                    if not os.path.exists(turbsimPath):  os.makedirs(turbsimPath)
                    
            # The following loop creates the turbsim files for low box. That should really only happen if inflowStr is `TurbSim`, as I have commented out below.
            # However, as of now I need to copy some files to obtain the limits of the turbsim box to be able to shift the coordinate system. Later on i delete these dirs
            # if inflowStr == 'TurbSim':
            #     for seed in range(nSeeds):
            #         seedPath = os.path.join(condPath, f'Seed_{seed}')
            #         if not os.path.exists(seedPath):  os.makedirs(seedPath)
            for seed in range(self.nSeeds):
                seedPath = os.path.join(condPath, f'Seed_{seed}')
                if not os.path.exists(seedPath):  os.makedirs(seedPath)
                
        # Get rid of duplicate entries due to the nature of the loop (equiv to only getting the first nCases entries)
        self.condDirList = condDirList
        self.caseDirList = sorted(list(set(caseDirList_)))
        assert self.caseDirList==caseDirList_[:self.nCases]



    def copyTurbineFilesForEachCase(self, writeFiles=True):

        if not self.templateFilesCreatedBool:
            raise SyntaxError('Template files not set. Call `setTemplateFilename` before calling this function.')

        # Loops on all conditions/cases creating DISCON and *Dyn files
        for cond in range(self.nConditions):
            for case in range(self.nCases):
                currPath = os.path.join(self.path, self.condDirList[cond], self.caseDirList[case])
        
                # Recover info about the current CondXX_*/CaseYY_*
                Vhub_ = self.allCond.sel(cond=cond)['vhub'].values
        
                # Update parameters to be changed in the *Dyn files
                self.HydroDynFile['WaveHs']     = self.bins.sel(wspd=Vhub_, method='nearest').WaveHs.values
                self.HydroDynFile['WaveTp']     = self.bins.sel(wspd=Vhub_, method='nearest').WaveTp.values
                self.HydroDynFile['WvHiCOffD']  = 2.0*np.pi/self.HydroDynFile['WaveTp']   # check with JJ
                self.HydroDynFile['WvLowCOffS'] = 2.0*np.pi/self.HydroDynFile['WaveTp']
                # !!! JJ: do i need to change the first 3 values of HD? can they be 'default'?
        
                self.ElastoDynFile['RotSpeed']   = self.bins.sel(wspd=Vhub_, method='nearest').RotSpeed.values
                self.ElastoDynFile['BlPitch(1)'] = self.bins.sel(wspd=Vhub_, method='nearest').BlPitch.values
                self.ElastoDynFile['BlPitch(2)'] = self.bins.sel(wspd=Vhub_, method='nearest').BlPitch.values
                self.ElastoDynFile['BlPitch(3)'] = self.bins.sel(wspd=Vhub_, method='nearest').BlPitch.values
        
                self.SElastoDynFile['RotSpeed']  = self.bins.sel(wspd=Vhub_, method='nearest').RotSpeed.values
                self.SElastoDynFile['BlPitch']   = self.bins.sel(wspd=Vhub_, method='nearest').BlPitch.values
        
                # Write updated DISCON and *Dyn files. 
                if writeFiles:
                    self.HydroDynFile.write(os.path.join(currPath, self.HDfilename))
                    status = os.system(f'cp {self.templatePath}/{self.controllerInputfilename} {currPath}/{self.controllerInputfilename}')
                
                # Depending on the controller, it might need to be in the same level as the fstf input file. The ideal solution would be give the full
                # path to the controller input, but we have no control over the compilation process and it is likely that a very long string with the full
                # path will get cut. So we need to give the relative path. We give the path as the current one, so here we create a link to ensure it will
                # work regardless of how the controller was compiled. There is no hard in having this extra link even if it's not needed.
                notepath = os.getcwd();  os.chdir(self.path)
                for seed in range(self.nSeeds):
                    try:
                        src = os.path.join('../', self.controllerInputfilename)
                        dst = os.path.join(self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', self.controllerInputfilename)
                        if writeFiles:
                            os.symlink(src, dst)                
                    except FileExistsError:
                        pass
                os.chdir(notepath)
        
                # Write InflowWind files
                self.InflowWindFile['WindType']       = 3
                self.InflowWindFile['PropagationDir'] = 0
                self.InflowWindFile['Filename_BTS']   = '"./TurbSim"'
                if writeFiles:
                    self.InflowWindFile.write( os.path.join(currPath,f'InflowWind.dat'))
        
                for t in range(self.nTurbines):
                    # Recover info about the current turbine in CondXX_*/CaseYY_
                    yaw_deg_     = self.allCases.sel(case=case, turbine=t)['yaw'].values
                    yaw_mis_deg_ = self.allCases.sel(case=case, turbine=t)['yawmis'].values
                    ADmodel_     = self.allCases.sel(case=case, turbine=t)['ADmodel'].values
                    EDmodel_     = self.allCases.sel(case=case, turbine=t)['EDmodel'].values
        
                    # Quickly check that yaw misaligned value is zero if case does not contain yaw misalignment
                    if self.allCases.sel(case=case, turbine=t)['misalignment'].values:
                        assert yaw_mis_deg_ != 0
                    else:
                        assert yaw_mis_deg_ == 0
        
                    if EDmodel_ == 'FED':
                        # Update each turbine's ElastoDyn
                        self.ElastoDynFile['NacYaw']   = yaw_deg_ + yaw_mis_deg_
                        # check if this change of EDfile to a variable works as it should with quotes and stuff
                        self.ElastoDynFile['BldFile1'] = self.ElastoDynFile['BldFile2'] = self.ElastoDynFile['BldFile3'] = f'"{self.bladefilename}"'
                        self.ElastoDynFile['TwrFile'] = f'"{self.towerfilename}"'
                        self.ElastoDynFile['Azimuth']  = round(np.random.uniform(low=0, high=360)) # start at a random value
                        if writeFiles:
                            self.ElastoDynFile.write(os.path.join(currPath,f'{self.EDfilename}{t+1}_mod.dat'))
        
                    elif EDmodel_ == 'SED':
                        # Update each turbine's Simplified ElastoDyn
                        self.SElastoDynFile['BlPitch']  = self.bins.sel(wspd=Vhub_, method='nearest').BlPitch.values
                        self.SElastoDynFile['RotSpeed'] = self.bins.sel(wspd=Vhub_, method='nearest').RotSpeed.values
                        self.SElastoDynFile['NacYaw']   = yaw_deg_ + yaw_mis_deg_
                        if writeFiles:
                            self.SElastoDynFile.write(os.path.join(currPath,f'{self.SEDfilename}{t+1}_mod.dat'))
        
                    # Update each turbine's ServoDyn
                    self.ServoDynFile['YawNeut']      = yaw_deg_ + yaw_mis_deg_
                    self.ServoDynFile['DLL_FileName'] = f'"{self.DLLfilepath}{t+1}.so"'
                    if writeFiles:
                        self.ServoDynFile.write( os.path.join(currPath,f'{self.SrvDfilename}{t+1}_mod.dat'))
        
                    # Update each turbine's OpenFAST input
                    self.turbineFile['TMax']         = self.tmax
                    self.turbineFile['CompInflow']   = 1  # 1: InflowWind;     2: OpenFoam (fully coupled; not VTK input to FF)
        
                    if EDmodel_ == 'FED':
                        self.turbineFile['CompElast']    = 1  # 1: full ElastoDyn; 2: full ElastoDyn + BeamDyn;  3: Simplified ElastoDyn
                        self.turbineFile['CompSub']      = 1
                        self.turbineFile['CompHydro']    = 1
                        self.turbineFile['EDFile']       = f'"./{self.EDfilename}{t+1}_mod.dat"'
                    elif EDmodel_ == 'SED':
                        self.turbineFile['CompElast']    = 3  # 1: full ElastoDyn; 2: full ElastoDyn + BeamDyn;  3: Simplified ElastoDyn
                        self.turbineFile['CompSub']      = 0  # need to be disabled with SED
                        self.turbineFile['CompHydro']    = 0  # need to be disabled with SED
                        self.turbineFile['IntMethod']    = 3
                        self.turbineFile['EDFile']       = f'"./{self.SEDfilename}{t+1}_mod.dat"'
                    self.turbineFile['BDBldFile(1)'] = f'"{self.BDfilepath}"'
                    self.turbineFile['BDBldFile(2)'] = f'"{self.BDfilepath}"'
                    self.turbineFile['BDBldFile(3)'] = f'"{self.BDfilepath}"'
                    self.turbineFile['InflowFile']   = f'"./{self.IWfilename}"'
                    if ADmodel_ == 'ADyn':
                        self.turbineFile['CompAero']     = 2  # 1: AeroDyn v14;    2: AeroDyn v15;   3: AeroDisk
                        self.turbineFile['AeroFile']     = f'"{self.ADfilepath}"'
                    elif ADmodel_ == 'ADsk':
                        # from Andy's email on 2022-11-01 So if you use AeroDisk with ElastoDyn, just set the blade DOFs to false for now.
                        self.turbineFile['CompAero']     = 3  # 1: AeroDyn v14;    2: AeroDyn v15;   3: AeroDisk
                        self.turbineFile['AeroFile']     = f'"{self.ADskfilepath}"'
                        if writeFiles:
                            status = os.system(f'cp {self.coeffTablefilepath} {os.path.join(currPath,self.coeffTablefilename)}')
                    self.turbineFile['ServoFile']    = f'"./{self.SrvDfilename}{t+1}_mod.dat"'
                    self.turbineFile['HydroFile']    = f'"./{self.HDfilename}"'
                    self.turbineFile['SubFile']      = f'"{self.SubDfilepath}"'
                    self.turbineFile['MooringFile']  = f'"unused"'
                    self.turbineFile['IceFile']      = f'"unused"'
                    self.turbineFile['TStart']       = 0 # start saving openfast output from time 0 (to see transient)
                    self.turbineFile['OutFileFmt']   = 3 # 1: .out; 2: .outb; 3: both
        
                    if writeFiles:
                        self.turbineFile.write( os.path.join(currPath,f'{self.turbfilename}{t+1}.fst'))


                                        
    def setTemplateFilename(self, templatePath=None, EDfilename=None, SEDfilename=None, HDfilename=None, SrvDfilename=None, ADfilename=None, ADskfilename=None, SubDfilename=None, IWfilename=None, BDfilepath=None, bladefilename=None, towerfilename=None, turbfilename=None, libdisconfilepath=None, controllerInputfilename=None, coeffTablefilename=None, turbsimLowfilepath=None, turbsimHighfilepath=None, FFfilename=None):
        '''                             
                                        
        *filename: str                  
        The filename of the current O   penFAST submodule, no complete path. Assumes it is
        inside `templatePath`           
        *filepath: str                  
        Complete path of the file. Ma   y or may not be inside `templatePath`

        '''

        self.EDfilename    = "unused";  self.EDfilepath    = "unused" 
        self.SEDfilename   = "unused";  self.SEDfilepath   = "unused" 
        self.HDfilename    = "unused";  self.HDfilepath    = "unused" 
        self.SrvDfilename  = "unused";  self.SrvDfilepath  = "unused" 
        self.ADfilename    = "unused";  self.ADfilepath    = "unused" 
        self.ADskfilename  = "unused";  self.ADskfilepath  = "unused" 
        self.SubDfilename  = "unused";  self.SubDfilepath  = "unused" 
        self.IWfilename    = "unused";  self.IWfilepath    = "unused" 
        self.BDfilepath    = "unused";  self.BDfilename    = "unused" 
        self.bladefilename = "unused";  self.bladefilepath = "unused" 
        self.towerfilename = "unused";  self.towerfilepath = "unused" 


        if templatePath is None:
            print(f'--- WARNING: No template files given. Complete setup will not be possible')
            return

        if not os.path.isdir(templatePath):
            raise ValueError (f'Template path {templatePath} does not seem to exist.')

        self.templatePath = templatePath


        def checkIfExists(f):
            if f == 'unused':
                return
            if not os.path.isfile(f):
                raise ValueError (f'File {f} does not exist.')

        if EDfilename is not None:
            if EDfilename != 'unused' and not EDfilename.endswith('.T'):
                raise ValueError (f'Name the template ED file "*.T.dat" and give "*.T" as `EDfilename`')
            self.EDfilepath = os.path.join(self.templatePath,f"{EDfilename}.dat")
            checkIfExists(self.EDfilepath)
            self.EDfilename = EDfilename

        if SEDfilename is not None:
            if SEDfilename != 'unused' and not SEDfilename.endswith('.T'):
                raise ValueError (f'Name the template SED file "*.T.dat" and give "*.T" as `SEDfilename`')
            self.SEDfilepath = os.path.join(self.templatePath,f"{SEDfilename}.dat")
            checkIfExists(self.SEDfilepath)
            self.SEDfilename = SEDfilename

        if HDfilename is not None:
            if HDfilename != 'unused' and not HDfilename.endswith('.dat'):
                raise ValueError (f'The HydroDyn filename should end in `.dat`.')
            self.HDfilepath = os.path.join(self.templatePath,HDfilename)
            checkIfExists(self.HDfilepath)
            self.HDfilename = HDfilename

        if SrvDfilename is not None:
            if not EDfilename.endswith('.T'):
                raise ValueError (f'Name the template ServoDyn file "*.T.dat" and give "*.T" as `SrvDfilename`')
            self.SrvDfilepath = os.path.join(self.templatePath,f"{SrvDfilename}.dat")
            checkIfExists(self.SrvDfilepath)
            self.SrvDfilename = SrvDfilename

        if ADfilename is not None:
            if ADfilename != 'unused' and not ADfilename.endswith('.dat'):
                raise ValueError (f'The AeroDyn filename should end in `.dat`.')
            self.ADfilepath = os.path.join(self.templatePath,ADfilename)
            checkIfExists(self.ADfilepath)
            self.ADfilename = ADfilename

        if ADskfilename is not None:
            if ADskfilename != 'unused' and not ADskfilename.endswith('.dat'):
                raise ValueError (f'The AeroDisk filename should end in `.dat`.')
            self.ADskfilepath = os.path.join(self.templatePath,ADskfilename)
            checkIfExists(self.ADskfilepath)
            self.ADskfilename = ADskfilename

        if SubDfilename is not None:
            if not SubDfilename.endswith('.dat'):
                raise ValueError (f'The SubDyn filename should end in `.dat`.')
            self.SubDfilepath = os.path.join(self.templatePath,SubDfilename)
            checkIfExists(self.SubDfilepath)
            self.SubDfilename = SubDfilename

        if IWfilename is not None:
            if IWfilename != 'unused' and not IWfilename.endswith('.dat'):
                raise ValueError (f'The InflowWind filename should end in `.dat`.')
            self.IWfilepath = os.path.join(self.templatePath,IWfilename)
            checkIfExists(self.IWfilepath)
            self.IWfilename = IWfilename

        if BDfilepath is not None:
            if BDfilepath != 'unused' and not BDfilepath.endswith('.dat'):
                raise ValueError (f'The BeamDyn filename should end in `.dat`.')
            self.BDfilepath = BDfilepath
            checkIfExists(self.BDfilepath)

        if bladefilename is not None:
            if not bladefilename.endswith('.dat'):
                raise ValueError (f'The blade filename should end in `.dat`.')
            self.bladefilepath = os.path.join(self.templatePath,bladefilename)
            checkIfExists(self.bladefilepath)
            self.bladefilename = bladefilename

        if towerfilename is not None:
            if not towerfilename.endswith('.dat'):
                raise ValueError (f'The tower filename should end in `.dat`.')
            self.towerfilepath = os.path.join(self.templatePath,towerfilename)
            checkIfExists(self.towerfilepath)
            self.towerfilename = towerfilename

        if turbfilename is not None:
            if not turbfilename.endswith('.T'):
                raise ValueError (f'Name the template turbine file "*.T.fst" and give "*.T" as `turbfilename`')
            self.turbfilepath = os.path.join(self.templatePath,f"{turbfilename}.fst")
            checkIfExists(self.turbfilepath)
            self.turbfilename = turbfilename

        if libdisconfilepath is not None:
            if not libdisconfilepath.endswith('.so'):
                raise ValueError (f'The libdiscon `libdisconfilepath` file should end in "*.so"')
            self.libdisconfilepath = libdisconfilepath
            checkIfExists(self.libdisconfilepath)
            self._create_copy_libdiscon()

        if controllerInputfilename is not None:
            if not controllerInputfilename.endswith('.IN'):
                print(f'--- WARNING: The controller input file typically ends in "*.IN". Currently {controllerInputfilename}. Double check.')
            self.controllerInputfilepath = os.path.join(self.templatePath, controllerInputfilename)
            checkIfExists(self.controllerInputfilepath)
            self.controllerInputfilename = controllerInputfilename

        if coeffTablefilename is not None:
            if not coeffTablefilename.endswith('.csv'):
                raise ValueError (f'The performance table `coeffTablefilename` file should end in "*.csv"')
            self.coeffTablefilepath = os.path.join(templatePath, coeffTablefilename)
            checkIfExists(self.coeffTablefilepath)
            self.coeffTablefilename = coeffTablefilename

        if turbsimLowfilepath is not None:
            if not turbsimLowfilepath.endswith('.inp'):
                raise ValueError (f'TurbSim file input for low-res box `turbsimLowfilepath` should end in ".inp".')
            self.turbsimLowfilepath = turbsimLowfilepath
            checkIfExists(self.turbsimLowfilepath)
        
        if turbsimHighfilepath is not None:
            if not turbsimHighfilepath.endswith('.inp'):
                raise ValueError (f'TurbSim file input for high-res box `turbsimHighfilepath` should end in ".inp".')
            self.turbsimHighfilepath = turbsimHighfilepath
            checkIfExists(self.turbsimHighfilepath)
        
        if FFfilename is not None:
            if not FFfilename.endswith('.fstf'):
                raise ValueError (f'FAST.Farm input file `FFfilename` should end in ".fstf".')
            self.FFfilepath = os.path.join(self.templatePath,FFfilename)
            checkIfExists(self.FFfilepath)
            self.FFfilename = FFfilename
        
        self._open_template_files()

        self.templateFilesCreatedBool = True

        return


    def _create_copy_libdiscon(self):
        # Make copies of libdiscon for each turbine if they don't exist
        copied = False
        for t in range(self.nTurbines):
            libdisconfilename = os.path.splitext(os.path.basename(self.libdisconfilepath))[0]
            currLibdiscon = os.path.join(os.path.dirname(self.libdisconfilepath), f'{libdisconfilename}.T{t+1}.so')
            self.DLLfilepath = os.path.join(os.path.dirname(self.libdisconfilepath), f'{libdisconfilename}.T')
            if not os.path.isfile(currLibdiscon):
                if self.verbose>0: print(f'    Creating a copy of the controller {libdisconfilename}.so in {currLibdiscon}')
                status = os.system(f'cp {self.libdisconfilepath} {currLibdiscon}')
                copied=True
        
        if copied == False and self.verbose>0:
            print(f'    Copies of the controller {libdisconfilename}.T[1-{self.nTurbines}].so already exists in {os.path.dirname(self.libdisconfilepath)}. Skipped step.')


    def _open_template_files(self):

        # Open template files
        def _check_and_open(f):
            if f != 'unused':
              return FASTInputFile(f)

        self.ElastoDynFile   = _check_and_open(self.EDfilepath)  
        self.SElastoDynFile  = _check_and_open(self.SEDfilepath) 
        self.HydroDynFile    = _check_and_open(self.HDfilepath)  
        self.ServoDynFile    = _check_and_open(self.SrvDfilepath)
        self.AeroDiskFile    = _check_and_open(self.ADskfilepath)
        self.turbineFile     = _check_and_open(self.turbfilepath)
        self.InflowWindFile  = _check_and_open(self.IWfilepath)  
            



    def print_template_files(self):
        raise NotImplementedError (f'Placeholder. Not implemented.')


    def createAuxArrays(self):

        self._rotate_wts()
        self._create_all_cond()
        self._create_all_cases()



    def _create_all_cond(self):

        if len(self.vhub)==len(self.shear) and len(self.shear)==len(self.TIvalue):
            self.nConditions = len(self.vhub)

            if self.verbose>1: print(f'The length of vhub, shear, and TI are the same. Assuming each position is a condition.')
            if self.verbose>0: print(f'Creating {self.nConditions} conditions')

            self.allCond = xr.Dataset({'vhub':    (['cond'], self.vhub   ),
                                       'shear':   (['cond'], self.shear  ),
                                       'TIvalue': (['cond'], self.TIvalue)},
                                       coords={'cond': np.arange(self.nConditions)} )
        
        else:
            import itertools
            self.nConditions = len(self.vhub) * len(self.shear) * len(self.TIvalue)

            if self.verbose>1: print(f'The length of vhub, shear, and TI are different. Assuming sweep on each of them.')
            if self.verbose>0: print(f'Creating {self.nConditions} conditions')
    
            # Repeat arrays as necessary to build xarray Dataset
            combination = np.vstack(list(itertools.product(self.vhub,self.shear,self.TIvalue)))

            self.allCond = xr.Dataset({'vhub':    (['cond'], combination[:,0]),
                                       'shear':   (['cond'], combination[:,1]),
                                       'TIvalue': (['cond'], combination[:,2])},
                                       coords={'cond': np.arange(self.nConditions)} )
  

    def _create_all_cases(self):
        # Generate the different "cases" (inflow angle, and misalignment and wakesteer bools). She reads from ../NewFFParams_Shell.csv
        # If misalignment true, then the actual yaw is yaw[turb]=np.random.uniform(low=-8.0, high=8.0). CaseSetup/ParameterManipulation.py:166
        
        # Calculate the total number of cases given sweeps requested. Multipliers for wake steering, yaw misalignment, and reduced-order models
        nWindDir = len(np.unique(self.inflow_deg))
        nCasesWSmultiplier = 2 if self.sweepWS else 1
        nCasesYMmultiplier = 2 if self.sweepYM else 1
        nCasesROmultiplier = len(self.EDmodel)
        
        # Yaw multiplier, setup in the form of repeated wind directions with changing yaw
        nCasesYawmultiplier = int(len(self.inflow_deg)/len(np.unique(self.inflow_deg)))
        
        # Aux nCases vars
        nCases        = int(nWindDir * nCasesWSmultiplier * nCasesYMmultiplier * nCasesROmultiplier * nCasesYawmultiplier)
        nCasesWSfalse = int(nCases/nCasesWSmultiplier)
        nCasesWStrue  = int(nCases - nCasesWSfalse)
        nCasesYMfalse = int(nCases/nCasesYMmultiplier)
        nCasesYMtrue  = int(nCases - nCasesYMfalse)
        if self.verbose>2: print(f'    Cases: nWindDir = {nWindDir}                                   ')
        if self.verbose>2: print(f'    Cases: nCases = {nCases}')
        if self.verbose>2: print(f'    Cases: nCasesWStrue = {nCasesWStrue}')
        if self.verbose>2: print(f'    Cases: nCasesWSfalse = {nCasesWSfalse}')
        if self.verbose>2: print(f'    Cases: nCasesYMtrue = {nCasesYMtrue}')
        if self.verbose>2: print(f'    Cases: nCasesYMfalse = {nCasesYMfalse}')
        
        # Build an array of wind directions, with repeated values to account for wake steering and yaw misalign bools, and ROM options
        windDir = np.repeat(self.inflow_deg, nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier)
        yawInit = np.repeat(self.yaw_init, nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier, axis=0)
        
        # Build arrays of wake steering and yaw misalignment bools (done this way for clarity)
        if self.sweepWS and self.sweepYM:
            wakeSteering = np.tile([False, True,  False, True], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
            misalignment = np.tile([False, False, True,  True], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
        elif self.sweepWS and not self.sweepYM:
            wakeSteering = np.tile([False, True ], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
            misalignment = np.tile([False, False], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
        elif not self.sweepWS and self.sweepYM:
            wakeSteering = np.tile([False, False], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
            misalignment = np.tile([False, True ], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
        elif not self.sweepWS and not self.sweepYM:
            wakeSteering = np.tile([False], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
            misalignment = np.tile([False], nWindDir*nCasesROmultiplier*nCasesYawmultiplier)
        
            
        # Create array of random numbers for yaw misalignment, and set it to zero where no yaw misalign is requested
        yawMisalignedValue = np.random.uniform(size = [nCases,self.nTurbines], low=-8.0, high=8.0)
        yawMisalignedValue[~misalignment,:] = 0
        
        # Count number of simplified models to add that information to the xarray. If their length is 1, it means they weren't requested
        if len(self.ADmodel) == 1:
            nADyn = self.nTurbines
        else:
            nADyn = [ self.ADmodel[i].count('ADyn') for i in range(len(self.ADmodel)) ]
        if len(self.EDmodel) == 1:
           nFED = self.nTurbines
        else:
            nFED = [ self.EDmodel[i].count('FED') for i in range(len(self.EDmodel)) ]
        
        # Come up with an ordered "yaw case" numbering for dir name
        yawCase =  np.arange(nCasesYawmultiplier)+1
        
        # Assemble main case dataset, containing turbine info
        self.nCases = nCases
        self.allCases = xr.Dataset(
            {
                'Tx':     (['case','turbine'], np.repeat(self.wts_rot_ds['x'].values   , nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier*nCasesYawmultiplier, axis=0)),
                'Ty':     (['case','turbine'], np.repeat(self.wts_rot_ds['y'].values   , nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier*nCasesYawmultiplier, axis=0)),
                'Tz':     (['case','turbine'], np.repeat(self.wts_rot_ds['z'].values   , nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier*nCasesYawmultiplier, axis=0)),
                'D':      (['case','turbine'], np.repeat(self.wts_rot_ds['D'].values   , nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier*nCasesYawmultiplier, axis=0)),
                'zhub':   (['case','turbine'], np.repeat(self.wts_rot_ds['zhub'].values, nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier*nCasesYawmultiplier, axis=0)),
                'yawmis': (['case','turbine'], yawMisalignedValue),
                'yaw':    (['case','turbine'], yawInit),
                'yawCase': (['case'], np.repeat(yawCase, nWindDir*nCasesWSmultiplier*nCasesYMmultiplier*nCasesROmultiplier)),
                'ADmodel': (['case','turbine'], np.tile(np.repeat(self.ADmodel, nCasesWSmultiplier*nCasesYMmultiplier*nCasesYawmultiplier, axis=0),(nWindDir,1)) ),
                'EDmodel': (['case','turbine'], np.tile(np.repeat(self.EDmodel, nCasesWSmultiplier*nCasesYMmultiplier*nCasesYawmultiplier, axis=0),(nWindDir,1)) ),
                'nFullAeroDyn':    (['case'], np.repeat(np.tile(nADyn, nWindDir), nCasesWSmultiplier*nCasesYMmultiplier*nCasesYawmultiplier)),
                'nFulllElastoDyn': (['case'], np.repeat(np.tile(nFED,  nWindDir), nCasesWSmultiplier*nCasesYMmultiplier*nCasesYawmultiplier)),
                'wakeSteering': (['case'], wakeSteering),
                'misalignment': (['case'], misalignment),
                'inflow_deg':   (['case'], windDir),
            },
            coords={
                'case':    range(nCases),
                'turbine': range(self.nTurbines),
            },
        )


    def _rotate_wts(self):
        # Calculate the rotated positions of the turbines wrt the reference turbine
        wts_rot={}
        for inflow in self.inflow_deg:
            for i , turb in self.wts.items():
                ref = self.wts[self.refTurb_rot]
  
                xori = self.wts[i]['x']
                x = ref['x'] + (self.wts[i]['x']-ref['x'])*cosd(inflow) - (self.wts[i]['y']-ref['y'])*sind(inflow)
                yori = self.wts[i]['y']
                y = ref['y'] - (self.wts[i]['x']-ref['x'])*sind(-inflow) + (self.wts[i]['y']-ref['y'])*cosd(-inflow)
                z = self.wts[i]['z']
                D = self.wts[i]['D']
                zhub = self.wts[i]['zhub']
  
                wts_rot[inflow,i] = {'x':x, 'y':y, 'z':z,
                                     'D':D, 'zhub':zhub,
                                    }

        self.wts_rot_ds = pd.DataFrame.from_dict(wts_rot, orient='index').to_xarray().rename({'level_0':'inflow_deg','level_1':'turbine'})
  

    def _determine_box_parameters(self):

        if self.dt_high_les is not None:
            # Box paramters given. Only one check is needed since it passed `checkInputs`
            return

        # todo: compute the boxes paraters
        raise NotImplementedError(f'The ability to automatically determine the box paraters is not implemented yet.')



  
    def _setRotorParameters(self):
  
  
        if self.D == 220: # 12 MW turbine
            self.bins = xr.Dataset({'WaveHs':      (['wspd'], [ 1.429, 1.429]), # 1.429 comes from Matt's hydrodyn input file
                               'WaveTp':      (['wspd'], [ 7.073, 7.073]), # 7.073 comes from Matt's hydrodyn input file
                               'RotSpeed':    (['wspd'], [ 4.0, 4.0]), # 4 rpm comes from Matt's ED input file
                               'BlPitch':     (['wspd'], [ 0.0, 0.0]), # 0 deg comes from Matt's ED input file
                               #'WvHiCOffD':   (['wspd'], [0,   0]), # 2nd order wave info. Unused for now
                               #'WvLowCOffS':  (['wspd'], [0,   0]), # 2nd order wave info. Unused for now
                              },  coords={'wspd': [10, 15]} )  # 15 m/s is 'else', since method='nearest' is used on the variable `bins`
            
        elif self.D == 250: # IEA 15 MW
            self.bins = xr.Dataset({'WaveHs':      (['wspd'], [1.172, 1.323, 1.523, 1.764, 2.255]),  # higher values on default input from the repository (4.52)
                               'WaveTp':      (['wspd'], [7.287, 6.963, 7.115, 6.959, 7.067]),  # higher values on default input from the repository (9.45)
                               'RotSpeed':    (['wspd'], [4.995, 6.087, 7.557, 7.557, 7.557]),
                               'BlPitch':     (['wspd'], [0.315, 0,     0.645, 7.6,   13.8 ]),
                               #'WvHiCOffD':   (['wspd'], [0,     0,     0,     0,     0    ]), # 2nd order wave info. Unused for now. 3.04292 from repo; 0.862 from KS
                               #'WvLowCOffS':  (['wspd'], [0,     0,     0,     0,     0    ]), # 2nd order wave info. Unused for now  0.314159 from repo; 0.862 from KS
                              },  coords={'wspd': [6.6, 8.6, 10.6, 12.6, 15]} )  # 15 m/s is 'else', since method='nearest' is used on the variable `bins`
            
        else:
            raise ValueError(f'Unknown turbine with diameter {self.D}. Add values to the `_setRotorParameters` function.')
  


    
    def TS_low_setup(self, writeFiles=True, runOnce=False):
        # Loops on all conditions/seeds creating Low-res TurbSim box  (following python-toolbox/pyFAST/fastfarm/examples/Ex1_TurbSimInputSetup.py)
        boxType='lowres'
        for cond in range(self.nConditions):
            for seed in range(self.nSeeds):
                seedPath = os.path.join(self.path, self.condDirList[cond], f'Seed_{seed}')
                
                # ---------------- TurbSim Low boxes setup ------------------ #
                # Set file to be created
                currentTSLowFile = os.path.join(seedPath, 'Low_stillToBeModified.inp')
                
                # Get properties needed for the creation of the low-res turbsim inp file
                D_       = self.allCases['D'   ].max().values
                HubHt_   = self.allCases['zhub'].max().values
                xlocs_   = self.allCases['Tx'  ].values.flatten() # All turbines are needed for proper
                ylocs_   = self.allCases['Ty'  ].values.flatten() # and consistent extent calculation
                Vhub_    = self.allCond.sel(cond=cond)['vhub'   ].values
                shear_   = self.allCond.sel(cond=cond)['shear'  ].values
                tivalue_ = self.allCond.sel(cond=cond)['TIvalue'].values
                # Coherence parameters
                a = 12;  b=0.12                            # IEC 61400-3 ed4, app C, eq C.16
                Lambda1 = 0.7*HubHt_ if HubHt_<60 else 42  # IEC 61400-3 ed4, sec 6.3.1, eq 5 
        
                # Create and write new Low.inp files creating the proper box with proper resolution
                # By passing low_ext, manual mode for the domain size is activated, and by passing ds_low, manual mode
                # for discretization (and further domain size) is also activated
                currentTS = TSCaseCreation(D_, HubHt_, Vhub_, tivalue_, shear_, x=xlocs_, y=ylocs_, zbot=self.zbot, cmax=self.cmax,
                                           fmax=self.fmax, Cmeander=self.Cmeander, boxType=boxType, low_ext=self.extent_low, ds_low=self.ds_low_les)
                self.TSlowbox = currentTS
                if runOnce:
                    return
                currentTS.writeTSFile(self.turbsimLowfilepath, currentTSLowFile, tmax=self.tmax, verbose=self.verbose)
                
                # Modify some values and save file (some have already been set in the call above)
                Lowinp = FASTInputFile(currentTSLowFile)
                Lowinp['RandSeed1'] = self.seedValues[seed]
                Lowinp['PLExp']     = shear_
                #Lowinp['latitude']  = latitude # Not used when IECKAI model is selected.
                Lowinp['InCDec1']   = Lowinp['InCDec2'] = Lowinp['InCDec3'] = f'"{a} {b/(8.1*Lambda1):.8f}"'
                # The dt was computed for a proper low-res box but here we will want to compare with the high-res
                # and it is convenient to have the same time step. Let's do that change here
                Lowinp['TimeStep']  = 1/(2*self.fmax)
                if writeFiles:
                    Lowinp.write( os.path.join(seedPath, 'Low.inp') )
                
                # Let's remove the original file
                os.remove(os.path.join(seedPath, 'Low_stillToBeModified.inp'))

        self.TSlowBoxFilesCreatedBool = True


    def TS_low_slurm_prepare(self, slurmfilepath):



       # # --------------------------------------------------
       # # ----- Prepare SLURM script for Low-res boxes -----
       # # --------------- ONE SCRIPT PER CASE --------------
       # # --------------------------------------------------
       # 
       # if not os.path.isfile(slurmfilepath):
       #     raise ValueError (f'SLURM script for FAST.Farm {slurmfilepath} does not exist.')
       # self.slurmfilename_ff = os.path.basename(slurmfilepath)


       # for cond in range(self.nConditions):
       #     for case in range(self.nCases):
       #         for seed in range(self.nSeeds):
       #             
       #             fname = f'runFASTFarm_cond{cond}_case{case}_seed{seed}.sh'
       #             status = os.system(f'cp {slurmfilepath} {os.path.join(self.path,fname)}')
       # 
       #             # Change job name (for convenience only)
       #             sed_command = f"sed -i 's|#SBATCH --job-name=runFF|#SBATCH --job-name=c{cond}_c{case}_s{seed}_runFF_{os.path.basename(self.path)}|g' {fname}"
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)
       #             # Change logfile name (for convenience only)
       #             sed_command = f"sed -i 's|#SBATCH --output log.fastfarm_c0_c0_seed0|#SBATCH --output log.fastfarm_c{cond}_c{case}_s{seed}|g' {fname}"
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)
       #             # Change the fastfarm binary to be called
       #             sed_command = f"""sed -i "s|^fastfarmbin.*|fastfarmbin='{self.ffbin}'|g" {fname}"""
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)
       #             # Change the path inside the script to the desired one
       #             sed_command = f"sed -i 's|/projects/shellwind/rthedin/Task2_2regis|{self.path}|g' {fname}"
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)
       #             # Write condition
       #             sed_command = f"""sed -i "s|^cond.*|cond='{self.condDirList[cond]}'|g" {fname}"""
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)
       #             # Write case
       #             sed_command = f"""sed -i "s|^case.*|case='{self.caseDirList[case]}'|g" {fname}"""
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)
       #             # Write seed
       #             sed_command = f"""sed -i "s|^seed.*|seed={seed}|g" {fname}"""
       #             _ = subprocess.call(sed_command, cwd=self.path, shell=True)



        # --------------------------------------------------
        # ----- Prepare SLURM script for Low-res boxes -----
        # --------------------------------------------------

        if not os.path.isfile(slurmfilepath):
            raise ValueError (f'SLURM script for low-res box {slurmfilepath} does not exist.')
        self.slurmfilename_low = os.path.basename(slurmfilepath)

        import subprocess
        status = os.system(f'cp {slurmfilepath} {self.path}/{self.slurmfilename_low}')
        
        # Change job name (for convenience only)
        _ = subprocess.call(f"sed -i 's|#SBATCH --job-name=lowBox|#SBATCH --job-name=lowBox_{os.path.basename(self.path)}|g' {self.slurmfilename_low}", cwd=self.path, shell=True)
        # Change the path inside the script to the desired one
        sed_command = f"sed -i 's|/projects/shellwind/rthedin/Task2_2regis|{self.path}|g' {self.slurmfilename_low}"
        _ = subprocess.call(sed_command, cwd=self.path, shell=True)
        # Change number of nodes values 
        _ = subprocess.call(f"sed -i 's|#SBATCH --nodes=2|#SBATCH --nodes={int(np.ceil(self.nConditions*self.nSeeds/6))}|g' {self.slurmfilename_low}", cwd=self.path, shell=True)
        # Assemble list of conditions and write it
        listtoprint = "' '".join(self.condDirList)
        sed_command = f"""sed -i "s|^condList.*|condList=('{listtoprint}')|g" {self.slurmfilename_low}"""
        _ = subprocess.call(sed_command, cwd=self.path, shell=True)
        # Change the number of seeds
        _ = subprocess.call(f"sed -i 's|nSeeds=6|nSeeds={self.nSeeds}|g' {self.slurmfilename_low}", cwd=self.path, shell=True)


        if self.nSeeds != 6:
            print(f'--- WARNING: The memory-per-cpu on the low-res boxes SLURM script is configured for 6 seeds, not {self.nSeeds}.')

    def TS_low_slurm_submit(self):
        # ---------------------------------
        # ----- Run turbSim Low boxes -----
        # ---------------------------------
        # Submit the script to SLURM
        _ = subprocess.call('sbatch {self.slurmfilename_low}', cwd=self.path, shell=True)


    def TS_low_createSymlinks(self):
        # Create symbolic links for all of the time-series and the Low.bts files too
        
        notepath = os.getcwd()
        os.chdir(self.path)
        for cond in range(self.nConditions):
            for case in range(self.nCases):
                for seed in range(self.nSeeds):
                    try:
                        src = os.path.join('../../../..', self.condDirList[cond], f'Seed_{seed}', 'Low.bts')
                        dst = os.path.join(self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', 'TurbSim', 'Low.bts')
                        os.symlink(src, dst)                
                    except FileExistsError:
                        print(f'    File {dst} already exists. Skipping symlink.')
        os.chdir(notepath)


    def getDomainParameters(self):

        # If the low box setup hasn't been called (e.g. LES run), do it once to get domain extents
        if not self.TSlowBoxFilesCreatedBool:
            self.TS_low_setup(writeFiles=False, runOnce=True)

        # Figure out how many (and which) high boxes actually need to be executed. Remember that wake steering, yaw misalignment, SED/ADsk models,
        # and sweep in yaw do not require extra TurbSim runs
        self.nHighBoxCases = len(np.unique(self.inflow_deg))  # some wind dir might be repeated for sweep on yaws
        
        self.allHighBoxCases = self.allCases.where(~self.allCases['wakeSteering'],drop=True).drop_vars('wakeSteering')\
                                            .where(~self.allCases['misalignment'], drop=True).drop_vars('misalignment')\
                                            .where(self.allCases['nFullAeroDyn']==self.nTurbines, drop=True).drop_vars('ADmodel')\
                                            .where(self.allCases['nFulllElastoDyn']==self.nTurbines, drop=True).drop_vars('EDmodel')\
                                            .where(self.allCases['yawCase']==1, drop=True).drop_vars('yawCase')
        
        if self.nHighBoxCases != len(self.allHighBoxCases.case):
            raise ValueError(f'The number of cases do not match as expected. {self.nHighBoxCases} unique wind directions, but {len(self.allHighBoxCases.case)} unique cases.')
        
        # Determine offsets from turbines coordinate frame to TurbSim coordinate frame
        self.yoffset_turbsOrigin2TSOrigin = -( (self.TSlowbox.ymax - self.TSlowbox.ymin)/2 + self.TSlowbox.ymin )
        self.xoffset_turbsOrigin2TSOrigin = - self.extent_low[0]*self.D
        
        if self.verbose>0:
            print(f"    The y offset between the turbine ref frame and turbsim is {self.yoffset_turbsOrigin2TSOrigin}")
            print(f"    The x offset between the turbine ref frame and turbsim is {self.xoffset_turbsOrigin2TSOrigin}")

        if self.verbose>2:
            print(f'allHighBoxCases is:')
            print(self.allHighBoxCases)


    def TS_high_get_time_series(self):

        # Loop on all conditions/seeds extracting time series from the Low box at turbines location
        boxType='highres'
        for cond in range(self.nConditions):
            for seed in range(self.nSeeds):
                condSeedPath = os.path.join(path, self.condDirList[cond], f'Seed_{seed}')
        
                # Read output .bts for current seed
                bts = TurbSimFile(os.path.join(condSeedPath, 'Low.bts'))
                bts['t']  = np.round(bts['t'],  6) # rounding single precision read as double precision
                bts['dt'] = np.round(bts['dt'], 6)
        
                for case in range(self.nHighBoxCases):
                    # Get actual case number given the high-box that need to be saved
                    case = self.allHighBoxCases.isel(case=case)['case'].values
                    
                    caseSeedPath = os.path.join(self.path, self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', 'TurbSim')
                    
                    for t in range(nTurbines):
                        # Recover turbine properties of the current case
                        HubHt_   = self.allCases.sel(case=case, turbine=t)['zhub'].values
                        xloc_    = self.allCases.sel(case=case, turbine=t)['Tx'  ].values
                        yloc_    = self.allCases.sel(case=case, turbine=t)['Ty'  ].values
        
                        # Turbine location in TurbSim reference frame
                        xt = xloc_ + self.xoffset_turbsOrigin2TSOrigin
                        yt = yloc_ + self.yoffset_turbsOrigin2TSOrigin
        
                        # Get indices of turbine location in TurbSim files
                        jTurb, kTurb = bts.closestPoint(y=yt,z=HubHt_)
                        Hub_series = bts['z'][kTurb]
        
                        # Get indices of the half height position (TurbSim's hub height)
                        jMid, kMid = bts.iMid
        
                        # Get time series at the box center to get mean vhub and create time array. JJ: get at the Turbsim hub height
                        #Vhub = bts['u'][0,:,jTurb,kTurb]
                        Vmid = bts['u'][0,:,jMid,kMid]
                        time = bts.t
        
                        # The time-series need to be shifted depending on the turbine location, so we need to find how many
                        # grid points (time steps) the data have convected. We use the mean streamwise component for that
                        start_time_step = round( (xt/Vmid.mean())/bts.dt ) 
        
                        # Get time-series given rolling
                        uvel = np.roll(bts['u'][0, :, jTurb, kTurb], start_time_step)
                        vvel = np.roll(bts['u'][1, :, jTurb, kTurb], start_time_step)
                        wvel = np.roll(bts['u'][2, :, jTurb, kTurb], start_time_step)
        
                        # Checks
                        assert len(time)==len(uvel)
                        assert len(uvel)==len(vvel)
                        assert len(vvel)==len(wvel)
        
                        # Save timeseries as CondXX/Seed_Z/USRTimeSeries_T*.txt. This file will later be copied to CondXX/CaseYY/Seed_Z
                        timeSeriesOutputFile = os.path.join(caseSeedPath, f'USRTimeSeries_T{t+1}.txt')
        
                        # The reference frame used in the time-series is the inertial frame of the high-res box (local).
                        # Sometimes the point where we want to place the turbine at exists and then we can set y=0. For example, suppose the low-res
                        # grid has y = ..., 980, 1000, 1020, ..., and we want to place a turbine at y=1000. The high-res box has 5m resolution. Then,
                        # the time-series will be pulled from _exactly_ y=1000, and since the low-res grid has a grid point there too, so we can put
                        # y=0 on the time-series input file. However, if we want the turbine at y=998, we must account for the difference since this 
                        # y value is not a grid point of the low-res box. So we compute an offset between the turbine location and the nearest grid
                        # point in the low-res box, and then pass this offset to the time-series file. In this example, the offset is 2 m, thus the
                        # time-series file will have a y of 2 m.
                        yoffset = bts['y'][jTurb] - yt
                        if yoffset != 0:
                            print(f"Seed {seed}, Case {case}: Turbine {t+1} is not at a grid point location. Tubine is at y={yloc_}",\
                                  f"on the turbine reference frame, which is y={yt} on the low-res TurbSim reference frame. The",\
                                  f"nearest grid point in y is {bts['y'][jTurb]} so printing y={yoffset} to the time-series file.")
                        writeTimeSeriesFile(timeSeriesOutputFile, yoffset, Hub_series, uvel, vvel, wvel, time)




    def TS_high_setup(self, writeFiles=True):

        #todo: Check if the low-res boxes were created successfully

        # Create symbolic links for the low-res boxes
        self.TS_low_createSymlinks()

        # Loop on all conditions/cases/seeds setting up the High boxes
        boxType='highres'
        for cond in range(self.nConditions):
            for case in range(self.nHighBoxCases):
                # Get actual case number given the high-box that need to be saved
                case = self.allHighBoxCases.isel(case=case)['case'].values
                if self.verbose>3:
                    print(f'Generating high-res box setup for cond {cond} ({self.condDirList[cond]}), case {case} ({self.caseDirList[case]}).')
                for seed in range(self.nSeeds):
                    seedPath = os.path.join(self.path, self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', 'TurbSim')
        
                    for t in range(self.nTurbines):
        
                        # ---------------- TurbSim High boxes setup ------------------ #
                        currentTSHighFile = os.path.join(seedPath, f'HighT{t+1}_stillToBeModified.inp')
        
                        # Get properties needed for the creation of the high-res turbsim inp file
                        D_       = self.allCases.sel(case=case, turbine=t)['D'   ].values
                        HubHt_   = self.allCases.sel(case=case, turbine=t)['zhub'].values
                        xloc_    = self.allCases.sel(case=case, turbine=t)['Tx'  ].values
                        yloc_    = self.allCases.sel(case=case, turbine=t)['Ty'  ].values
                        Vhub_    = self.allCond.sel(cond=cond)['vhub'   ].values
                        shear_   = self.allCond.sel(cond=cond)['shear'  ].values
                        tivalue_ = self.allCond.sel(cond=cond)['TIvalue'].values
                        
                        # Coherence parameters
                        a = 12;  b=0.12                            # IEC 61400-3 ed4, app C, eq C.16
                        Lambda1 = 0.7*HubHt_ if HubHt_<60 else 42  # IEC 61400-3 ed4, sec 6.3.1, eq 5 
        
                        # Create and write new Low.inp files creating the proper box with proper resolution
                        currentTS = TSCaseCreation(D_, HubHt_, Vhub_, tivalue_, shear_, x=xloc_, y=yloc_, zbot=self.zbot,
                                                   cmax=self.cmax, fmax=self.fmax, Cmeander=self.Cmeander, boxType=boxType, high_ext=self.extent_high)
                        currentTS.writeTSFile(self.turbsimHighfilepath, currentTSHighFile, tmax=self.tmax, turb=t, verbose=self.verbose)
        
                        # Modify some values and save file (some have already been set in the call above)
                        Highinp = FASTInputFile(currentTSHighFile)
                        Highinp['RandSeed1'] = self.seedValues[seed]
                        Highinp['TimeStep']  = 1/(2*self.fmax)
                        Highinp['TurbModel'] = f'"TIMESR"'
                        Highinp['UserFile']  = f'"USRTimeSeries_T{t+1}.txt"'
                        Highinp['RefHt']     = HubHt_
                        Highinp['URef']      = Vhub_
                        Highinp['PLExp']     = shear_
                        #Highinp['latitude']  = latitude # Not used when IECKAI model is selected.
                        Highinp['InCDec1']   = Highinp['InCDec2'] = Highinp['InCDec3'] = f'"{a} {b/(8.1*Lambda1):.8f}"'
                        if writeFiles:
                            Highinp.write( os.path.join(seedPath, f'HighT{t+1}.inp') )
                
                        # Let's remove the original file
                        os.remove(os.path.join(seedPath, f'HighT{t+1}_stillToBeModified.inp'))
        

    def TS_high_slurm_prepare(self, slurmfilepath):
        # ---------------------------------------------------
        # ----- Prepare SLURM script for High-res boxes -----
        # ---------------------------------------------------
        
        if not os.path.isfile(slurmfilepath):
            raise ValueError (f'SLURM script for high-res box {slurmfilepath} does not exist.')
        self.slurmfilename_high = os.path.basename(slurmfilepath)

        ntasks = self.nConditions*self.nHighBoxCases*self.nSeeds*self.nTurbines
        status = os.system(f'cp {slurmfilepath} {self.path}/{self.slurmfilename_high}')
        
        # Change job name (for convenience only)
        _ = subprocess.call(f"sed -i 's|#SBATCH --job-name=highBox|#SBATCH --job-name=highBox_{os.path.basename(self.path)}|g' {self.slurmfilename_high}", cwd=self.path, shell=True)
        # Change the path inside the script to the desired one
        sed_command = f"sed -i 's|/projects/shellwind/rthedin/Task2_2regis|{self.path}|g' {self.slurmfilename_high}"
        _ = subprocess.call(sed_command, cwd=self.path, shell=True)
        # Change number of turbines
        _ = subprocess.call(f"sed -i 's|nTurbines=12|nTurbines={self.nTurbines}|g' {self.slurmfilename_high}", cwd=self.path, shell=True)
        # Change number of seeds
        _ = subprocess.call(f"sed -i 's|nSeeds=6|nSeeds={self.nSeeds}|g' {self.slurmfilename_high}", cwd=self.path, shell=True)
        # Change number of nodes values
        _ = subprocess.call(f"sed -i 's|#SBATCH --nodes=3|#SBATCH --nodes={int(np.ceil(ntasks/36))}|g' {self.slurmfilename_high}", cwd=self.path, shell=True)
        # Assemble list of conditions and write it
        listtoprint = "' '".join(self.condDirList)
        sed_command = f"""sed -i "s|^condList.*|condList=('{listtoprint}')|g" {self.slurmfilename_high}"""
        _ = subprocess.call(sed_command, cwd=self.path, shell=True)
        # Assemble list of cases and write it
        highBoxesCaseDirList = [self.caseDirList[c] for c in self.allHighBoxCases.case.values]
        listtoprint = "' '".join(highBoxesCaseDirList)
        sed_command = f"""sed -i "s|^caseList.*|caseList=('{listtoprint}')|g" {self.slurmfilename_high}"""
        _ = subprocess.call(sed_command, cwd=self.path, shell=True)



    def TS_high_slurm_submit(self):
        # ----------------------------------
        # ----- Run turbSim High boxes -----
        # ----------------------------------
        # Submit the script to SLURM
        _ = subprocess.call('sbatch {self.slurmfilename_high}', cwd=self.path, shell=True)

    
    def TS_high_create_symlink(self):
        # Create symlink of all the high boxes for the cases with wake steering and yaw misalignment. These are the "repeated" boxes
        
        notepath = os.getcwd()
        os.chdir(self.path)
        for cond in range(self.nConditions):
            for t in range(self.nTurbines):
                for seed in range(self.nSeeds):
                    for case in range(self.nCases):
                        # Let's check if the current case is source (has bts) or destination (needs a symlink to bts)
                        varsToDrop = ['wakeSteering','misalignment','yawmis','yaw','yawCase','ADmodel','EDmodel','nFullAeroDyn','nFulllElastoDyn']
                        if case in self.allHighBoxCases['case']:
                            src = os.path.join('../../../..', self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', 'TurbSim', f'HighT{t+1}.bts')
                            xr_src = self.allCases.sel(case=case, drop=True).drop_vars(varsToDrop)
                            continue
                        else:
                            dst = os.path.join(self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', 'TurbSim', f'HighT{t+1}.bts')
                            xr_dst = self.allCases.sel(case=case, drop=True).drop_vars(varsToDrop)
                            
                        # Let's make sure the src and destination are the same case, except wake steering and yaw misalign bools
                        xr.testing.assert_equal(xr_src, xr_dst)
                        
                        try:
                            os.symlink(src, dst)
                        except FileExistsError:
                            if self.verbose>1: print(f'File {dst} already exists. Skipping symlink.')
        os.chdir(notepath)



    def FF_setup(self, outlistFF=None, **kwargs):
        '''

        **kwargs:
            seedsToKeep: int
                For the LES setup. Often 1, but if you want to run multiple times the same thing, pick a different value
        '''

        if outlistFF is None:
            # Output list for FAST.Farm runs. Use 1 at the end for turbines (they will be replicated for all turbines). Combination of sample input file and KS's input
            outlistFF = [
                "RtAxsXT1     , RtAxsYT1     , RtAxsZT1",
                "RtPosXT1     , RtPosYT1     , RtPosZT1",
                #"RtDiamT1",
                "YawErrT1",
                "TIAmbT1",
                'RtVAmbT1',
                'RtVRelT1',
                'RtSkewT1',
                'RtCtAvgT1',
                'W1VAmbX, W1VAmbY, W1VAmbZ',
                "W1VDisX, W1VDisY, W1VDisZ",
                "CtT1N01      , CtT1N02      , CtT1N03      , CtT1N04      , CtT1N05      , CtT1N06      , CtT1N07      , CtT1N08      , CtT1N09      , CtT1N10      , CtT1N11      , CtT1N12      , CtT1N13      , CtT1N14      , CtT1N15      , CtT1N16      , CtT1N17      , CtT1N18      , CtT1N19      ,  CtT1N20",
            #    "WkAxsXT1D1   , WkAxsXT1D2   , WkAxsXT1D3   , WkAxsXT1D4   , WkAxsXT1D5   , WkAxsXT1D6   , WkAxsXT1D7",
            #    "WkAxsYT1D1   , WkAxsYT1D2   , WkAxsYT1D3   , WkAxsYT1D4   , WkAxsYT1D5   , WkAxsYT1D6   , WkAxsYT1D7",
            #    "WkAxsZT1D1   , WkAxsZT1D2   , WkAxsZT1D3   , WkAxsZT1D4   , WkAxsZT1D5   , WkAxsZT1D6   , WkAxsZT1D7",
            #    "WkPosXT1D1   , WkPosXT1D2   , WkPosXT1D3   , WkPosXT1D4   , WkPosXT1D5   , WkPosXT1D6   , WkPosXT1D7",
            #    "WkPosYT1D1   , WkPosYT1D2   , WkPosYT1D3   , WkPosYT1D4   , WkPosYT1D5   , WkPosYT1D6   , WkPosYT1D7",
            #    "WkPosZT1D1   , WkPosZT1D2   , WkPosZT1D3   , WkPosZT1D4   , WkPosZT1D5   , WkPosZT1D6   , WkPosZT1D7",
            #    "WkDfVxT1N01D1, WkDfVxT1N02D1, WkDfVxT1N03D1, WkDfVxT1N04D1, WkDfVxT1N05D1, WkDfVxT1N06D1, WkDfVxT1N07D1, WkDfVxT1N08D1, WkDfVxT1N09D1, WkDfVxT1N10D1, WkDfVxT1N11D1, WkDfVxT1N12D1, WkDfVxT1N13D1, WkDfVxT1N14D1, WkDfVxT1N15D1, WkDfVxT1N16D1, WkDfVxT1N17D1, WkDfVxT1N18D1, WkDfVxT1N19D1, WkDfVxT1N20D1",
            #    "WkDfVxT1N01D2, WkDfVxT1N02D2, WkDfVxT1N03D2, WkDfVxT1N04D2, WkDfVxT1N05D2, WkDfVxT1N06D2, WkDfVxT1N07D2, WkDfVxT1N08D2, WkDfVxT1N09D2, WkDfVxT1N10D2, WkDfVxT1N11D2, WkDfVxT1N12D2, WkDfVxT1N13D2, WkDfVxT1N14D2, WkDfVxT1N15D2, WkDfVxT1N16D2, WkDfVxT1N17D2, WkDfVxT1N18D2, WkDfVxT1N19D2, WkDfVxT1N20D2",
            #    "WkDfVxT1N01D3, WkDfVxT1N02D3, WkDfVxT1N03D3, WkDfVxT1N04D3, WkDfVxT1N05D3, WkDfVxT1N06D3, WkDfVxT1N07D3, WkDfVxT1N08D3, WkDfVxT1N09D3, WkDfVxT1N10D3, WkDfVxT1N11D3, WkDfVxT1N12D3, WkDfVxT1N13D3, WkDfVxT1N14D3, WkDfVxT1N15D3, WkDfVxT1N16D3, WkDfVxT1N17D3, WkDfVxT1N18D3, WkDfVxT1N19D3, WkDfVxT1N20D3",
            #    "WkDfVxT1N01D4, WkDfVxT1N02D4, WkDfVxT1N03D4, WkDfVxT1N04D4, WkDfVxT1N05D4, WkDfVxT1N06D4, WkDfVxT1N07D4, WkDfVxT1N08D4, WkDfVxT1N09D4, WkDfVxT1N10D4, WkDfVxT1N11D4, WkDfVxT1N12D4, WkDfVxT1N13D4, WkDfVxT1N14D4, WkDfVxT1N15D4, WkDfVxT1N16D4, WkDfVxT1N17D4, WkDfVxT1N18D4, WkDfVxT1N19D4, WkDfVxT1N20D4",
            #    "WkDfVxT1N01D5, WkDfVxT1N02D5, WkDfVxT1N03D5, WkDfVxT1N04D5, WkDfVxT1N05D5, WkDfVxT1N06D5, WkDfVxT1N07D5, WkDfVxT1N08D5, WkDfVxT1N09D5, WkDfVxT1N10D5, WkDfVxT1N11D5, WkDfVxT1N12D5, WkDfVxT1N13D5, WkDfVxT1N14D5, WkDfVxT1N15D5, WkDfVxT1N16D5, WkDfVxT1N17D5, WkDfVxT1N18D5, WkDfVxT1N19D5, WkDfVxT1N20D5",
            #    "WkDfVxT1N01D6, WkDfVxT1N02D6, WkDfVxT1N03D6, WkDfVxT1N04D6, WkDfVxT1N05D6, WkDfVxT1N06D6, WkDfVxT1N07D6, WkDfVxT1N08D6, WkDfVxT1N09D6, WkDfVxT1N10D6, WkDfVxT1N11D6, WkDfVxT1N12D6, WkDfVxT1N13D6, WkDfVxT1N14D6, WkDfVxT1N15D6, WkDfVxT1N16D6, WkDfVxT1N17D6, WkDfVxT1N18D6, WkDfVxT1N19D6, WkDfVxT1N20D6",
            #    "WkDfVxT1N01D7, WkDfVxT1N02D7, WkDfVxT1N03D7, WkDfVxT1N04D7, WkDfVxT1N05D7, WkDfVxT1N06D7, WkDfVxT1N07D7, WkDfVxT1N08D7, WkDfVxT1N09D7, WkDfVxT1N10D7, WkDfVxT1N11D7, WkDfVxT1N12D7, WkDfVxT1N13D7, WkDfVxT1N14D7, WkDfVxT1N15D7, WkDfVxT1N16D7, WkDfVxT1N17D7, WkDfVxT1N18D7, WkDfVxT1N19D7, WkDfVxT1N20D7",
            #    "WkDfVrT1N01D1, WkDfVrT1N02D1, WkDfVrT1N03D1, WkDfVrT1N04D1, WkDfVrT1N05D1, WkDfVrT1N06D1, WkDfVrT1N07D1, WkDfVrT1N08D1, WkDfVrT1N09D1, WkDfVrT1N10D1, WkDfVrT1N11D1, WkDfVrT1N12D1, WkDfVrT1N13D1, WkDfVrT1N14D1, WkDfVrT1N15D1, WkDfVrT1N16D1, WkDfVrT1N17D1, WkDfVrT1N18D1, WkDfVrT1N19D1, WkDfVrT1N20D1",
            #    "WkDfVrT1N01D2, WkDfVrT1N02D2, WkDfVrT1N03D2, WkDfVrT1N04D2, WkDfVrT1N05D2, WkDfVrT1N06D2, WkDfVrT1N07D2, WkDfVrT1N08D2, WkDfVrT1N09D2, WkDfVrT1N10D2, WkDfVrT1N11D2, WkDfVrT1N12D2, WkDfVrT1N13D2, WkDfVrT1N14D2, WkDfVrT1N15D2, WkDfVrT1N16D2, WkDfVrT1N17D2, WkDfVrT1N18D2, WkDfVrT1N19D2, WkDfVrT1N20D2",
            #    "WkDfVrT1N01D3, WkDfVrT1N02D3, WkDfVrT1N03D3, WkDfVrT1N04D3, WkDfVrT1N05D3, WkDfVrT1N06D3, WkDfVrT1N07D3, WkDfVrT1N08D3, WkDfVrT1N09D3, WkDfVrT1N10D3, WkDfVrT1N11D3, WkDfVrT1N12D3, WkDfVrT1N13D3, WkDfVrT1N14D3, WkDfVrT1N15D3, WkDfVrT1N16D3, WkDfVrT1N17D3, WkDfVrT1N18D3, WkDfVrT1N19D3, WkDfVrT1N20D3",
            #    "WkDfVrT1N01D4, WkDfVrT1N02D4, WkDfVrT1N03D4, WkDfVrT1N04D4, WkDfVrT1N05D4, WkDfVrT1N06D4, WkDfVrT1N07D4, WkDfVrT1N08D4, WkDfVrT1N09D4, WkDfVrT1N10D4, WkDfVrT1N11D4, WkDfVrT1N12D4, WkDfVrT1N13D4, WkDfVrT1N14D4, WkDfVrT1N15D4, WkDfVrT1N16D4, WkDfVrT1N17D4, WkDfVrT1N18D4, WkDfVrT1N19D4, WkDfVrT1N20D4",
            #    "WkDfVrT1N01D5, WkDfVrT1N02D5, WkDfVrT1N03D5, WkDfVrT1N04D5, WkDfVrT1N05D5, WkDfVrT1N06D5, WkDfVrT1N07D5, WkDfVrT1N08D5, WkDfVrT1N09D5, WkDfVrT1N10D5, WkDfVrT1N11D5, WkDfVrT1N12D5, WkDfVrT1N13D5, WkDfVrT1N14D5, WkDfVrT1N15D5, WkDfVrT1N16D5, WkDfVrT1N17D5, WkDfVrT1N18D5, WkDfVrT1N19D5, WkDfVrT1N20D5",
            #    "WkDfVrT1N01D6, WkDfVrT1N02D6, WkDfVrT1N03D6, WkDfVrT1N04D6, WkDfVrT1N05D6, WkDfVrT1N06D6, WkDfVrT1N07D6, WkDfVrT1N08D6, WkDfVrT1N09D6, WkDfVrT1N10D6, WkDfVrT1N11D6, WkDfVrT1N12D6, WkDfVrT1N13D6, WkDfVrT1N14D6, WkDfVrT1N15D6, WkDfVrT1N16D6, WkDfVrT1N17D6, WkDfVrT1N18D6, WkDfVrT1N19D6, WkDfVrT1N20D6",
            #    "WkDfVrT1N01D7, WkDfVrT1N02D7, WkDfVrT1N03D7, WkDfVrT1N04D7, WkDfVrT1N05D7, WkDfVrT1N06D7, WkDfVrT1N07D7, WkDfVrT1N08D7, WkDfVrT1N09D7, WkDfVrT1N10D7, WkDfVrT1N11D7, WkDfVrT1N12D7, WkDfVrT1N13D7, WkDfVrT1N14D7, WkDfVrT1N15D7, WkDfVrT1N16D7, WkDfVrT1N17D7, WkDfVrT1N18D7, WkDfVrT1N19D7, WkDfVrT1N20D7",
            ]
        self.outlistFF = outlistFF


        # Planes to save in FAST.Farm. We want the planes through the original farm, so let's get the position of the turbines at wdir=0
        alignedTurbs = self.allCases.where(self.allCases['inflow_deg']==0, drop=True).isel(case=0)
        if self.inflowStr == 'TurbSim':
            # Turbine location in TurbSim reference frame
            xWT = alignedTurbs['Tx'].values + self.xoffset_turbsOrigin2TSOrigin
            yWT = alignedTurbs['Ty'].values + self.yoffset_turbsOrigin2TSOrigin
        elif self.inflowStr == 'LES':
            # Turbine location in LES reference frame
            xWT = alignedTurbs['Tx'].values
            yWT = alignedTurbs['Ty'].values
        
        offset=50
        planes_xy = [self.zhub+self.zbot]
        planes_yz = np.unique(xWT+offset)
        planes_xz = np.unique(yWT)
        
        # Number of planes must be at most 9
        self.planes_xy = planes_xy[0:9]
        self.planes_yz = planes_yz[0:9]
        self.planes_xz = planes_xz[0:9]      
        

        if self.inflowStr == 'LES':
            self._FF_setup_LES(**kwargs)

        elif self.inflowStr == 'TurbSim':
            # We need to make sure the TurbSim boxes have been executed. Let's check the last line of the logfile
            highboxlog_path = os.path.join(self.path, self.condDirList[0], self.caseDirList[0], 'Seed_0', 'TurbSim', 'log.hight1.seed0.txt')
            if not os.path.isfile(highboxlog_path):
                raise ValueError(f'All TurbSim boxes need to be completed before this step can be done.')

            with open(highboxlog_path) as f:
                last = None
                for last in (line for line in f if line.rstrip('\n')):  pass

            if 'TurbSim terminated normally' not in last:
                raise ValueError(f'All TurbSim boxes need to be completed before this step can be done.')

            self._FF_setup_TS(**kwargs)


    def _FF_setup_LES(self, seedsToKeep=1):

        self.seedsToKeep = seedsToKeep

        # Clean unnecessary directories and files created by the general setup
        for cond in range(self.nConditions):
            for seed in range(self.nSeeds):
                os.system(f"rm -rf {os.path.join(self.path,self.condDirList[cond],f'Seed_{seed}')}")
    
            for case in range(self.nCases):
                #os.system(f"rm -rf {os.path.join(path,condDirList[cond],caseDirList[case], f'Seed_0','InflowWind.dat')}") # needs to exist
                for seed in range(seedsToKeep,self.nSeeds):
                    os.system(f"rm -rf {os.path.join(self.path,self.condDirList[cond],self.caseDirList[case], f'Seed_{seed}')}")
                        
                        
        
       # Create symlinks for the processed-and-renamed vtk files
        LESboxesDirName = 'LESboxes'
        
        for cond in range(self.nConditions):
            for case in range(self.nCases):
                for seed in range(self.seedsToKeep):
                    # Remove TurbSim dir and create LES boxes dir
                    _ = os.system(f"rm -rf {os.path.join(self.path,self.condDirList[cond],self.caseDirList[case], f'Seed_{seed}', 'TurbSim')}")
                    if not os.path.exists(os.path.join(self.path,self.condDirList[cond],self.caseDirList[case],f'Seed_{seed}',LESboxesDirName)):
                        os.makedirs(os.path.join(self.path,self.condDirList[cond],self.caseDirList[case],f'Seed_{seed}',LESboxesDirName))
        
                    # Low-res box
                    try:
                        src = os.path.join(self.LESpath, 'Low')
                        dst = os.path.join(self.path, self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', LESboxesDirName, 'Low')
                        os.symlink(src, dst)                
                    except FileExistsError:
                        print(f'Directory {dst} already exists. Skipping symlink.')
        
                    # High-res boxes
                    for t in range(self.nTurbines):
                        try:
                            src = os.path.join(self.LESpath, f"HighT{t+1}_inflow{str(self.allCases.sel(case=case).inflow_deg.values).replace('-','m')}deg")
                            dst = os.path.join(self.path,self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}', LESboxesDirName, f'HighT{t+1}')
                            os.symlink(src, dst)                
                        except FileExistsError:
                            print(f'Directory {dst} already exists. Skipping symlink.')
        
                        # Make a link for Amb.t0.vtk pointing to Amb.t1.vtk 
        
        
        
        # Loops on all conditions/cases and cases for FAST.Farm  (following python-toolbox/pyFAST/fastfarm/examples/Ex2_FFarmInputSetup.py)
        for cond in range(self.nConditions):
            for case in range(self.nCases):
                for seed in range(seedsToKeep):
                    seedPath = os.path.join(self.path, self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}')
        
                    # Recover case properties
                    D_       = self.allCases['D'   ].max().values # Getting the maximum in case different turbines are present
                    Vhub_    = self.allCond.sel(cond=cond)['vhub'   ].values
                    # Recover turbine properties (array of length nTurbines)
                    xWT = self.allCases.sel(case=case)['Tx'].values
                    yWT = self.allCases.sel(case=case)['Ty'].values
                    zWT = self.allCases.sel(case=case)['Tz'].values             
        
                    # --------------- FAST.Farm ----------------- #
                    templateFSTF = os.path.join(self.templatePath, self.FFfilename)
                    outputFSTF   = os.path.join(seedPath, 'FFarm_mod.fstf')
        
                    # Write the file (mostly for turbine locations here
                    writeFastFarm(outputFSTF, templateFSTF, xWT, yWT, zWT, FFTS=None, OutListT1=self.outlistFF, noLeadingZero=True)
        
                    # Open saved file and change additional values manually or make sure we have the correct ones
                    ff_file = FASTInputFile(outputFSTF)
        
                    # Open output file and change additional values manually or make sure we have the correct ones
                    ff_file['InflowFile']  = f'"../{self.IWfilename}"'   
                    ff_file['Mod_AmbWind'] = 1  # 1: LES boxes; 2: single TurbSim; 3: multiple TurbSim
                    ff_file['TMax'] = self.tmax
        
                    # LES-related parameters
                    ff_file['DT_Low-VTK']   = self.dt_low_les
                    ff_file['DT_High-VTK']  = self.dt_high_les
                    ff_file['WindFilePath'] = f'''"{os.path.join(seedPath, LESboxesDirName)}"'''
                    #if checkWindFiles:
                    #    ff_file['ChkWndFiles'] = 'TRUE'
        
                    # Super controller
                    ff_file['UseSC'] = False
                    ff_file['SC_FileName'] = '/path/to/SC_DLL.dll'
        
                    # Wake dynamics
                    ff_file['dr'] = self.cmax
                    ff_file['NumRadii']  = int(np.ceil(3*D_/(2*self.cmax) + 1))
                    ff_file['NumPlanes'] = int(np.ceil( 20*D_/(self.dt_low_les*Vhub_*(1-1/6)) ) )
        
                    # Vizualization outputs
                    ff_file['WrDisWind'] = 'False'
                    ff_file['WrDisDT']   = ff_file['DT_Low-VTK']*10    # writeFastFarm sets this to be the same as DT_Low
                    ff_file['NOutDisWindXY'] = len(self.planes_xy)
                    ff_file['OutDisWindZ']   = ' '.join(map(str, self.planes_xy))
                    ff_file['NOutDisWindYZ'] = len(self.planes_yz)
                    ff_file['OutDisWindX']   = ' '.join(map(str, self.planes_yz))
                    ff_file['NOutDisWindXZ'] = len(self.planes_xz)
                    ff_file['OutDisWindY']   = ' '.join(map(str, self.planes_xz))
        
                    # Modify wake outputs
                    ff_file['NOutDist'] = 7
                    ff_file['OutDist']  = ' '.join(map(str,  [1,1.5,2,2.5,3,3.5,4]*D_))
                    # Mofidy wind output
                    ff_file['NWindVel'] = 9
                    ff_file['WindVelX'] = ' '.join(map(str, xWT[:9]))
                    ff_file['WindVelY'] = ' '.join(map(str, yWT[:9]))
                    ff_file['WindVelZ'] = ' '.join(map(str, zWT[:9]))
        
                    ff_file.write(outputFSTF)

        # Update the number of seeds variable for the LES case
        self.nSeeds = self.seedsToKeep



    def _FF_setup_TS(self):

        # Loops on all conditions/cases and cases for FAST.Farm  (following python-toolbox/pyFAST/fastfarm/examples/Ex2_FFarmInputSetup.py)
        for cond in range(self.nConditions):
            for case in range(self.nCases):
                for seed in range(self.nSeeds):
                    seedPath = os.path.join(self.path, self.condDirList[cond], self.caseDirList[case], f'Seed_{seed}')
        
                    # Recover case properties
                    D_       = self.allCases['D'   ].max().values # Getting the maximum in case different turbines are present
                    HubHt_   = self.allCases['zhub'].max().values # Getting the maximum in case different turbines are present
                    Vhub_    = self.allCond.sel(cond=cond)['vhub'   ].values
                    shear_   = self.allCond.sel(cond=cond)['shear'  ].values
                    tivalue_ = self.allCond.sel(cond=cond)['TIvalue'].values
                    # Recover turbine properties (array of length nTurbines)
                    xWT = self.allCases.sel(case=case)['Tx'].values
                    yWT = self.allCases.sel(case=case)['Ty'].values
                    zWT = self.allCases.sel(case=case)['Tz'].values  
        
                    # # Turbine location in TurbSim reference frame
                    xt = xWT + self.xoffset_turbsOrigin2TSOrigin
                    yt = yWT + self.yoffset_turbsOrigin2TSOrigin
                    
        
                    # --------------- FAST.Farm ----------------- #
                    templateFSTF = os.path.join(self.templatePath, self.FFfilename)
                    outputFSTF   = os.path.join(seedPath, 'FFarm_mod.fstf')
        
                    # Open TurbSim outputs for the Low box and one High box (they are all of the same size)
                    lowbts = TurbSimFile(os.path.join(seedPath,'TurbSim', 'Low.bts'))
                    highbts  = TurbSimFile(os.path.join(seedPath,'TurbSim', f'HighT1.bts'))
        
                    # Get dictionary with all the D{X,Y,Z,t}, L{X,Y,Z,t}, N{X,Y,Z,t}, {X,Y,Z}0
                    dt_low_desired = self.Cmeander*D_/(10*Vhub_) # will be made multiple of dT_High inside _getBoxesParamsForFF
                    d = self._getBoxesParamsForFF(lowbts, highbts, dt_low_desired, D_, HubHt_, xWT, yt)
        
                    # Write the file
                    writeFastFarm(outputFSTF, templateFSTF, xWT, yt, zWT, d, OutListT1=self.outlistFF, noLeadingZero=True)
        
        
                    # Open saved file and change additional values manually or make sure we have the correct ones
                    ff_file = FASTInputFile(outputFSTF)
                    ff_file['InflowFile'] = f'"../{self.IWfilename}"'     #!!!!!!!! this path is not filled. should it be?
                    #ff_file['DT']=1.0
                    ff_file['Mod_AmbWind'] = 3  # 1: LES boxes; 2: single TurbSim; 3: multiple TurbSim
                    ff_file['TMax'] = self.tmax
        
                    # Super controller
                    ff_file['UseSC'] = False
                    ff_file['SC_FileName'] = '/path/to/SC_DLL.dll'
        
                    # Wake dynamics # !!!!!!!!!!!!!!!!!!!!!! are these values good? (Emmanuel's sample file had 3, 50, 80. KS had 5, 75, 80)
                    ff_file['dr'] = self.cmax
                    ff_file['NumRadii']  = int(np.ceil(3*D_/(2*self.cmax) + 1))
                    ff_file['NumPlanes'] = int(np.ceil( 20*D_/(dt_low_desired*Vhub_*(1-1/6)) ) )
        
                    # Vizualization outputs
                    ff_file['WrDisWind'] = 'False'
                    ff_file['WrDisDT']   = ff_file['DT_Low']*10    # writeFastFarm sets this to be the same as DT_Low
                    ff_file['NOutDisWindXY'] = len(self.planes_xy)
                    ff_file['OutDisWindZ']   = ' '.join(map(str, self.planes_xy))
                    ff_file['NOutDisWindYZ'] = len(self.planes_yz)
                    ff_file['OutDisWindX']   = ' '.join(map(str, self.planes_yz))
                    ff_file['NOutDisWindXZ'] = len(self.planes_xz)
                    ff_file['OutDisWindY']   = ' '.join(map(str, self.planes_xz))
        
                    # Modify wake outputs
                    ff_file['NOutDist'] = 7
                    ff_file['OutDist']  = ' '.join(map(str,  [1,1.5,2,2.5,3,3.5,4]*D_))
                    # Mofidy wind output # !!!!! JJ why only 9?
                    ff_file['NWindVel'] = 9
                    ff_file['WindVelX'] = ' '.join(map(str, xWT[:9]))
                    ff_file['WindVelY'] = ' '.join(map(str, yWT[:9]))
                    ff_file['WindVelZ'] = ' '.join(map(str, zWT[:9]+self.zhub))
        
                    ff_file.write(outputFSTF)

        return



    def _getMultipleOf(self, val, multipleof):
        valmult = int(val/multipleof)*multipleof
        return round(valmult, 4)


    def _getBoxesParamsForFF(self, lowbts, highbts, dt_low_desired, D, HubHt, xWT, yt):
        # Get mean wind speeds at the half height location (advection speed)
        _, meanU_High =  highbts.midValues() # !!!!!!!!!!!!!!!! JJ: does it make sense to get both? the meanu for low will be a lot higher than vhub,
        _, meanU_Low  =  lowbts.midValues()  # !!!!!!!!!!!!!!!! JJ: and the meanu for high can be a bit higher than vhub
    
        dT_High = np.round(highbts.dt, 4) # 0.3
        # dX_High can sometimes be too high. So get the closest to the cmax, but multiple of what should have been
        dX_High = round(meanU_High*dT_High)
        if self.verbose>1:
            print(f'original dX_High is {dX_High}')
        dX_High = round(self.cmax/dX_High) * dX_High
        if self.verbose>1:
            print(f'after adjusting to closes multiple of cmax, dX_High is {dX_High}')
        dY_High = highbts.y[1] - highbts.y[0]
        dZ_High = highbts.z[1] - highbts.z[0]
    
    
        # ----- Low
        dT_Low = self._getMultipleOf(dt_low_desired, multipleof=dT_High)
        dX_Low = self._getMultipleOf(meanU_Low*dT_Low, multipleof=dX_High) # dy and dz high are 5. 
        dY_Low = lowbts.y[1] - lowbts.y[0]
        dZ_Low = lowbts.z[1] - lowbts.z[0]
    
        LY_Low = lowbts.y[-1]-lowbts.y[0]
        LZ_Low = lowbts.z[-1]-lowbts.z[0]
        LT_Low = np.round(lowbts.t[-1]-lowbts.t[0], 4)
    
        X0_Low = np.floor( (min(xWT) - self.extent_low[0]*D ))# - dX_Low)) #  # JJ!!!!!!! # removing EB's -dX_Low from the X0_Low specification
        X0_Low = self._getMultipleOf(X0_Low, multipleof=dX_Low)
        Y0_Low = np.floor( -LY_Low/2                   )   # Starting on integer value for aesthetics
        Z0_Low = lowbts.z[0]                               # we start at lowest to include tower
    
        XMax_Low = self._getMultipleOf(max(xWT) + self.extent_low[1]*D, multipleof=dX_Low)
        LX_Low = XMax_Low-X0_Low
    
        nX_Low = int(np.ceil(LX_Low/dX_Low)+1)  # plus 1 from the guidance
        nY_Low = len(lowbts.y)    # !!!!!!! JJ: different from what EB has
        nZ_Low = len(lowbts.z)
    
        assert nY_Low == int(np.ceil(LY_Low/dY_Low)+1)
        assert nZ_Low == int(np.ceil(LZ_Low/dZ_Low)+1)
        assert (nY_Low-1)*dY_Low == LY_Low
        assert (nZ_Low-1)*dZ_Low == LZ_Low
    
    
        # ----- High 
        Z0_High   = highbts.z[0]        # we start at lowest to include tower
    
        LX_High = self.extent_high*D
        LY_High = highbts.y[-1]-highbts.y[0]
        LZ_High = highbts.z[-1]-highbts.z[0]
        LT_High = np.round(highbts.t[-1]-highbts.t[0], 4)
        # !!!!!!!!!!!!!!!!! JJ: EB had these below. Why is that?
        #Max_High = HubHt+extent_high*D  # !!!!!!!!!!!!!!!!!! actual hh or midpoint?
        #LY_High =  min(LY_Box, extent_YZ*D      ) # Bounding to not exceed the box dimension
        #LZ_High =  min(LZ_Box, ZMax_High-Z0_High) # Bounding to not exceed the box dimension
    
        nX_High = int(np.ceil(LX_High/dX_High) + 1)  # plus 1 from the guidance
        nY_High = len(highbts.y)
        nZ_High = len(highbts.z)
    
        assert nY_High == int(np.ceil(LY_High/dY_High)+1)
        assert nZ_High == int(np.ceil(LZ_High/dZ_High)+1)
        assert (nY_High-1)*dY_High == LY_High
        assert (nZ_High-1)*dZ_High == LZ_High
    
        # --- High-res location per turbine
        X0_desired = np.asarray(xWT)-LX_High/2  # high-res is centered on turbine location
        Y0_desired = np.asarray(yt)-LY_High/2   # high-res is centered on turbine location
        X0_High    = X0_Low + np.floor((X0_desired-X0_Low)/dX_High)*dX_High
        Y0_High    = Y0_Low + np.floor((Y0_desired-Y0_Low)/dY_High)*dY_High
    
        if self.verbose>2:
            print(f'  Low Box  \t\t  High box   ')
            print(f'dT_Low: {dT_Low}\t\t dT_High: {dT_High}')
            print(f'dX_Low: {dX_Low}\t\t dX_High: {dX_High}')
            print(f'dY_Low: {dY_Low}\t\t dY_High: {dY_High}')
            print(f'dZ_Low: {dZ_Low}\t\t dZ_High: {dZ_High}')
            print(f'LX_Low: {LX_Low}\t\t LX_High: {LX_High}')
            print(f'LY_Low: {LY_Low}\t\t LY_High: {LY_High}')
            print(f'LZ_Low: {LZ_Low}\t\t LZ_High: {LZ_High}')
            print(f'LT_Low: {LT_Low}\t\t LT_High: {LT_High}')
            print(f'nX_Low: {nX_Low}\t\t nX_High: {nX_High}')
            print(f'nY_Low: {nY_Low}\t\t nY_High: {nY_High}')
            print(f'nZ_Low: {nZ_Low}\t\t nZ_High: {nZ_High}')
            print(f'X0_Low: {X0_Low}\t\t X0_High: {X0_High}')
            print(f'Y0_Low: {Y0_Low}  \t Y0_High: {Y0_High}')
            print(f'Z0_Low: {Z0_Low}\t\t Z0_High: {Z0_High}')
    
    
        # Fill dictionary with all values
        d = dict()
        d['DT_Low']      = np.around(dT_Low ,4)
        d['DT_High']     = np.around(dT_High,4)
        d['NX_Low']      = nX_Low
        d['NY_Low']      = nY_Low
        d['NZ_Low']      = nZ_Low
        d['X0_Low']      = np.around(X0_Low,4)
        d['Y0_Low']      = np.around(Y0_Low,4)
        d['Z0_Low']      = np.around(Z0_Low,4)
        d['dX_Low']      = np.around(dX_Low,4)
        d['dY_Low']      = np.around(dY_Low,4)
        d['dZ_Low']      = np.around(dZ_Low,4)
        d['NX_High']     = nX_High
        d['NY_High']     = nY_High
        d['NZ_High']     = nZ_High
        # --- High extent info for turbine outputs
        d['dX_High']     = np.around(dX_High,4)
        d['dY_High']     = np.around(dY_High,4)
        d['dZ_High']     = np.around(dZ_High,4)
        d['X0_High']     = np.around(X0_High,4)
        d['Y0_High']     = np.around(Y0_High,4)
        d['Z0_High']     = np.around(Z0_High,4)
        # --- Misc
        d['U_mean_Low']  = meanU_Low
        d['U_mean_High'] = meanU_High
    
    
        # --- Sanity check: check that the high res is at "almost" an integer location
        X_rel = (np.array(d['X0_High'])-d['X0_Low'])/d['dX_High']
        Y_rel = (np.array(d['Y0_High'])-d['Y0_Low'])/d['dY_High']
        dX = X_rel - np.round(X_rel) # Should be close to zero
        dY = Y_rel - np.round(Y_rel) # Should be close to zero
    
        if any(abs(dX)>1e-3):
            print('Deltas:',dX)
            raise Exception('Some X0_High are not on an integer multiple of the high-res grid')
        if any(abs(dY)>1e-3):
            print('Deltas:',dY)
            raise Exception('Some Y0_High are not on an integer multiple of the high-res grid')
            
        return d





    def FF_slurm_prepare(self, slurmfilepath):
        # ----------------------------------------------
        # ----- Prepare SLURM script for FAST.Farm -----
        # ------------- ONE SCRIPT PER CASE ------------
        # ----------------------------------------------
        
        if not os.path.isfile(slurmfilepath):
            raise ValueError (f'SLURM script for FAST.Farm {slurmfilepath} does not exist.')
        self.slurmfilename_ff = os.path.basename(slurmfilepath)


        for cond in range(self.nConditions):
            for case in range(self.nCases):
                for seed in range(self.nSeeds):
                    
                    fname = f'runFASTFarm_cond{cond}_case{case}_seed{seed}.sh'
                    status = os.system(f'cp {slurmfilepath} {os.path.join(self.path,fname)}')
        
                    # Change job name (for convenience only)
                    sed_command = f"sed -i 's|#SBATCH --job-name=runFF|#SBATCH --job-name=c{cond}_c{case}_s{seed}_runFF_{os.path.basename(self.path)}|g' {fname}"
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)
                    # Change logfile name (for convenience only)
                    sed_command = f"sed -i 's|#SBATCH --output log.fastfarm_c0_c0_seed0|#SBATCH --output log.fastfarm_c{cond}_c{case}_s{seed}|g' {fname}"
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)
                    # Change the fastfarm binary to be called
                    sed_command = f"""sed -i "s|^fastfarmbin.*|fastfarmbin='{self.ffbin}'|g" {fname}"""
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)
                    # Change the path inside the script to the desired one
                    sed_command = f"sed -i 's|/projects/shellwind/rthedin/Task2_2regis|{self.path}|g' {fname}"
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)
                    # Write condition
                    sed_command = f"""sed -i "s|^cond.*|cond='{self.condDirList[cond]}'|g" {fname}"""
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)
                    # Write case
                    sed_command = f"""sed -i "s|^case.*|case='{self.caseDirList[case]}'|g" {fname}"""
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)
                    # Write seed
                    sed_command = f"""sed -i "s|^seed.*|seed={seed}|g" {fname}"""
                    _ = subprocess.call(sed_command, cwd=self.path, shell=True)



    def FF_slurm_submit(self):

        # ----------------------------------
        # ---------- Run FAST.Farm ---------
        # ------- ONE SCRIPT PER CASE ------
        # ----------------------------------
        import time
            
        for cond in range(self.nConditions):
            for case in range(self.nCases):
                for seed in range(self.nSeeds):
                    
                    # Submit the script to SLURM
                    fname = f'runFASTFarm_cond{cond}_case{case}_seed{seed}.sh'
                    sub_command = f"sbatch {fname}"
                    print(f'Calling: {sub_command}')
                    subprocess.call(sub_command, cwd=self.path, shell=True)
                    time.sleep(4) # Sometimes the same job gets submitted twice. This gets around it.
            




















































