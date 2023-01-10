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


    def __init__(self, wts, cmax, fmax, Cmeander, tmax, zbot, vhub, shear, TIvalue, path, dt_high_les, ds_high_les, extent_high, dt_low_les, ds_low_les, extent_low, ffbin, ADmodel=None, EDmodel=None, yaw=None, nSeeds=6, LESpath=None, sweepWakeSteering=False, sweepYawMisalignment=False, seedValues=None):
        '''
        
        ffbin: str
            Full path of the FAST.Farm binary to be executed
        '''

        self.wts         = wts
        self.cmax        = cmax
        self.fmax        = fmax
        self.Cmeander    = Cmeander
        self.tmax        = tmax
        self.zbot        = zbot
        self.vhub        = vhub
        self.shear       = shear
        self.TIvalue     = TIvalue
        self.path        = path
        self.dt_high_les = dt_high_les
        self.ds_high_les = ds_high_les
        self.dt_low_les  = dt_low_les
        self.ds_low_les  = ds_low_les
        self.extent_low  = extent_low
        self.extent_high = extent_high
        self.ffbin       = ffbin
        self.ADmodel     = ADmodel
        self.EDmodel     = EDmodel
        self.yaw         = yaw
        self.nSeeds      = nSeeds
        self.LESpath     = LESpath
        self.sweepWS     = sweepWakeSteering
        self.sweepYM     = sweepYawMisalignment
        self.seedValues  = seedValues
        

        self._setAndCheckInputs()


        self.createAuxArrays()


  def createAuxArrays(self):



  def _setAndCheckInputs(self):

      if not  isinstance(self.wts,dict):
          raise ValueError (f'`wts` needs to be a dictionary with the following entries for each turbine: x, y, z, D, zhub.')
      if not isinstance(self.wts[0]['x'],(float,int)):
          raise ValueError (f'The `x` value for the turbines should be an integer or float')
      # todo: check if all entries of the dict are the same regarding D and zhub
      self.D = self.wts[0]['D']


      # Auxiliary variables and some checks
      if self.seedValues is None:
          self.seedValues = [2318573, 122299, 123456, 389432, -432443, 9849898]


      self.nTurbines = len(self.wts)


      self._setRotorParameters()

      # Ensure quantities are list
      self.vhub    = [vhub]    if isinstance(vhub,(float,int))    else vhub
      self.shear   = [shear]   if isinstance(shear,(float,int))   else shear
      self.TIvalue = [TIvalue] if isinstance(TIvalue,(float,int)) else TIvalue

      if self.ADmodel is None:
          self.ADmodel = np.tile(['ADyn'],(1,self.nTurbines))

      if self.EDmodel is None:
          self.EDmodel = np.tile(['FED'],(1,self.nTurbines))

      if self.yaw is None:
          self.yaw = np.zeros((1,self.nTurbines)))

      if not os.path.isfile(self.ffbin):
          raise ValueError (f'The FAST.Farm binary given does not appear to exist')



      if self.LESpath is None:
          print('Setting up FAST.Farm based on TurbSim inflow')
          self.inflowStr = 'TurbSim'
      else:
          print('Seeing up FAST.Farm based on LES inflow')
          if not os.path.isdir(self.LESpath):
              raise ValueError (f'The path {self.LESpath} does not exist')
          self.inflowStr = 'LES'



      # Create case path is doesn't exist
      if not os.path.exists(self.path):
          os.makedirs(self.path)
          
      # Perform some checks
      if len(self.seedValues) != self.nSeeds:
          raise ValueError(f'Number of seeds is {self.nSeeds} but {len(self.seedValues)} seed values were given. Adjust the seedValues array accordingly')
          
      if not (np.array(self.extent_low)>=0).all():
          raise ValueError(f'The array for low-res box extents should be given with positive values')

      if self.dt_low_les < self.dt_high_les:
          raise ValueError(f'The temporal resolution dT_High should not be greater than dT_Low on the LES side')

      if self.ds_low_les < ds_high_les:
          raise ValueError(f'The grid resolution dS_High should not be greater than dS_Low on the LES side')

      assert isinstance(self.extent_high, (float,int))
      if self.extent_high<=0:
          raise ValueError(f'The extent of high boxes should be positive')
          
      #assert dt_low_desired%(1/(2*fmax)) < 1e-10


      if np.shape(self.ADmodel) != np.shape(self.EDmodel):
          raise ValueError('Every case should have the aerodynamic and elastic model selected. The number of cases (lines) in `ADmodel` and `EDmodel` should be the same')
          
      if self.nTurbines != np.shape(self.ADmodel)[1]:
          raise ValueError(f'The number of turbines in wts ({len(wts)}) should match the number of turbines in the ADmodel and EDmodel arrays ({np.shape(ADmodel)[1]})')
          
      if len(inflow_deg) != len(yaw_init):
          raise ValueError(f'One row for each inflow angle should be given in yaw_init. Currently {len(inflow_deg)} inflow angle(s) and {len(yaw_init)} yaw entrie(s)')




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









