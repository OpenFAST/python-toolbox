# -*- coding: utf-8 -*-
"""
  - The x-y extents of the box are large enough to accomodate all turbines
  - The z extent is large enough to include the rotor, start from the user specified `zbot`,
      and accomodates for the meandering in the vertical direction .
  - The dy and dz resolution is set of the maximum chord of the turbine
  - The dt resolution is set based on the maximum frequency expected to be relevant for dynamics

@author: kshaler
"""
import os, glob, struct
import numpy as np

class TSCaseCreation:
    
    def __init__(self, D, HubHt, Vhub, TI, PLexp, x, y, z=None, zbot=1.0, cmax=5.0, fmax=5.0, Cmeander=1.9):
        """
        Instantiate the object. 
        
        Parameters
        ----------
        D       :   float,
                   rotor diameter (m)
        HubHt   :   float,
                   turbine hub height (m)
        Vhub	:   float,
                   mean wind speed at hub height (m/s)
        TI	:   float,
                   turbulence intensity at hub height
        PLexp	:   float,
                   power law exponent for shear (-)
        x, y, z   :   float,
               x-, y-, and z-location of turbine, respectively
                  if z is None, z is set to 0
        cmax    :   float,
                   maximum blade chord (m). If not specified, set to NREL 5MW value.
        fmax    :   float,
                   maximum excitation frequency (Hz). If not specified set to NREL 5MW tower value.
        """
        # Turbine parameters
        self.Turb(D, HubHt, cmax, fmax)
        # Turbine location
        self.turbLocs(x,y,z)
        # Discretization
        self.discretization(Vhub, TI, PLexp)
        # Setup domain size
        self.domainSize(zbot=zbot, Cmeander=Cmeander)

    def Turb(self, D, HubHt, cmax=5.0, fmax=5.0):
        """
        Define turbine parameters
        
        Parameters
        __________
        D       :   float,
                   rotor diameter (m)
        HubHt   :   float,
                   turbine hub height (m)
        tpath   :   string,
                   path to base turbine location (.fst)
        cmax    :   float,
                   maximum blade chord (m). If not specified, set to NREL 5MW value.
        fmax    :   float,
                   maximum excitation frequency (Hz). If not specified set to NREL 5MW tower value.
        """
        
        self.D     = D
        self.RefHt = HubHt
        self.cmax  = cmax
        self.fmax  = fmax
        
    def turbLocs(self,x,y,z=None):
        """
        Specify turbine locations
        
        Parameters
        ----------
        x, y, z   :   float,
               x-, y-, and z-location of turbine, respectively
        """
        self.x     = np.asarray(x)
        self.y     = np.asarray(y)
        if z is None:
            self.z     = np.asarray(y)*0
        else:
            self.z     = np.asarray(z)

    def discretization(self,Vhub,TI,Shear):
        
        self.URef = Vhub
        self.TI = TI
        self.PLexp = Shear
        
        # Derived properties
        self.dt = 1.0/(2.0*self.fmax)
        self.dy = self.cmax
        self.dz = self.cmax

    def domainSize(self, zbot, Cmeander=1.9):
    
        ymin = min(self.y)-2.23313*Cmeander*self.D/2
        ymax = max(self.y)+2.23313*Cmeander*self.D/2
        
        width_des = ymax-ymin
        height_des = self.RefHt+self.D/2+2.23313*Cmeander*self.D/2
        
        self.ny = round(width_des/self.dy)+1
        self.nz = round(height_des/self.dz)+1
        
        self.Width  = self.dy*(self.ny-1)
        self.Height = self.dz*(self.nz-1)
        
        Dgrid=min(self.Height,self.Width)
        self.HubHt_for_TS = zbot-0.5*Dgrid+self.Height


    def plotSetup(self, fig=None, ax=None):
        """ Plot a figure showing the turbine locations and the extent of the turbulence box"""
        if fig is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(6,5))
            ax  = fig.add_subplot(111,aspect="equal")
        xmin = min(self.x)-self.D
        xmax = max(self.x)+self.D

        ymid = np.mean(self.y)
        ymin = -self.Width/2 + ymid
        ymax =  self.Width/2 + ymid


        # high-res boxes
        for wt in range(len(self.x)):
            ax.plot(self.x[wt],self.y[wt],'x',ms=8,mew=2,label="WT{0}".format(wt+1))

        # low-res box
        #         ax.plot([xmin,xmax,xmax,xmin,xmin],
        #                 [ymin,ymin,ymax,ymax,ymin],'--k',lw=2,label='Low')
        ax.axhline(ymin, ls='--', c='k', lw=2, label='Low')
        ax.axhline(ymax, ls='--', c='k', lw=2)

        ax.legend(bbox_to_anchor=(1.05,1.015),frameon=False)
        ax.set_xlabel("x-location [m]")
        ax.set_ylabel("y-location [m]")
        fig.tight_layout
        return fig, ax

    def writeTSFile(self, fileIn, fileOut, NewFile=True, tpath=None, tmax=50):
        """ Write a TurbSim primary input file, 
        See WriteTSFile below.
        """
        WriteTSFile(fileIn, fileOut, self, NewFile=NewFile, tpath=tpath, tmax=tmax)


        
def WriteTSFile(fileIn, fileOut, params, NewFile=True, tpath=None, tmax=50):
    """ Write a TurbSim primary input file, 

        tpath:     string,
                   path to base turbine location (.fst)
                   only used if NewFile is False

    """

    if NewFile == True:
        print('Writing a new {0} file from scratch'.format(fileOut))
        # --- Writing FFarm input file from scratch
        with open(fileOut, 'w') as f:
            f.write('--------TurbSim v2.00.* Input File------------------------\n')
            f.write('for Certification Test #1 (Kaimal Spectrum, formatted FF files).\n')
            f.write('---------Runtime Options-----------------------------------\n')
            f.write('False\tEcho\t\t- Echo input data to <RootName>.ech (flag)\n')
            f.write('123456\tRandSeed1\t\t- First random seed  (-2147483648 to 2147483647)\n')
            f.write('RanLux\tRandSeed2\t\t- Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n')
            f.write('False\tWrBHHTP\t\t- Output hub-height turbulence parameters in binary form?  (Generates RootName.bin)\n')
            f.write('False\tWrFHHTP\t\t- Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat)\n')
            f.write('False\tWrADHH\t\t- Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh)\n')
            f.write('True\tWrADFF\t\t- Output full-field time-series data in TurbSim/AeroDyn form? (Generates RootName.bts)\n')
            f.write('False\tWrBLFF\t\t- Output full-field time-series data in BLADED/AeroDyn form?  (Generates RootName.wnd)\n')
            f.write('False\tWrADTWR\t\t- Output tower time-series data? (Generates RootName.twr)\n')
            f.write('False\tWrFMTFF\t\t- Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, RootName.v, RootName.w)\n')
            f.write('False\tWrACT\t\t- Output coherent turbulence time steps in AeroDyn form? (Generates RootName.cts)\n')
            f.write('True\tClockwise\t\t- Clockwise rotation looking downwind? (used only for full-field binary files - not necessary for AeroDyn)\n')
            f.write('0\tScaleIEC\t\t- Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale uniformly; 2=use individual scales]\n')
            f.write('\n')
            f.write('--------Turbine/Model Specifications-----------------------\n')
            f.write('{:.0f}\tNumGrid_Z\t\t- Vertical grid-point matrix dimension\n'.format(params.nz))
            f.write('{:.0f}\tNumGrid_Y\t\t- Horizontal grid-point matrix dimension\n'.format(params.ny))
            f.write('{:.6f}\tTimeStep\t\t- Time step [seconds]\n'.format(params.dt))
            f.write('{:.4f}\tAnalysisTime\t\t- Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n'.format(tmax))
            f.write('"ALL"\tUsableTime\t\t- Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds unless UsableTime is "ALL")\n')
            f.write('{:.3f}\tHubHt\t\t- Hub height [m] (should be > 0.5*GridHeight)\n'.format(params.HubHt_for_TS))
            f.write('{:.3f}\tGridHeight\t\t- Grid height [m]\n'.format(params.Height))
            f.write('{:.3f}\tGridWidth\t\t- Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n'.format(params.Width))
            f.write('0\tVFlowAng\t\t- Vertical mean flow (uptilt) angle [degrees]\n')
            f.write('0\tHFlowAng\t\t- Horizontal mean flow (skew) angle [degrees]\n')
            f.write('\n')
            f.write('--------Meteorological Boundary Conditions-------------------\n')
            f.write('"IECKAI"\tTurbModel\t\t- Turbulence model ("IECKAI","IECVKM","GP_LLJ","NWTCUP","SMOOTH","WF_UPW","WF_07D","WF_14D","TIDAL","API","IECKAI","TIMESR", or "NONE")\n')
            f.write('"unused"\tUserFile\t\t- Name of the file that contains inputs for user-defined spectra or time series inputs (used only for "IECKAI" and "TIMESR" models)\n')
            f.write('1\tIECstandard\t\t- Number of IEC 61400-x standard (x=1,2, or 3 with optional 61400-1 edition number (i.e. "1-Ed2") )\n')
            f.write('"{:.3f}\t"\tIECturbc\t\t- IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP model, not used for other models)\n'.format(params.TI))
            f.write('"NTM"\tIEC_WindType\t\t- IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n')
            f.write('"default"\tETMc\t\t- IEC Extreme Turbulence Model "c" parameter [m/s]\n')
            f.write('"PL"\tWindProfileType\t\t- Velocity profile type ("LOG";"PL"=power law;"JET";"H2L"=Log law for TIDAL model;"API";"PL";"TS";"IEC"=PL on rotor disk, LOG elsewhere; or "default")\n')
            f.write('"PowerLaw_6ms02.dat"\tProfileFile\t\t- Name of the file that contains input profiles for WindProfileType="PL" and/or TurbModel="USRVKM" [-]\n')
            f.write('{:.3f}\tRefHt\t\t- Height of the reference velocity (URef) [m]\n'.format(params.RefHt))
            f.write('{:.3f}\tURef\t\t- Mean (total) velocity at the reference height [m/s] (or "default" for JET velocity profile) [must be 1-hr mean for API model; otherwise is the mean over AnalysisTime seconds]\n'.format(params.URef))
            f.write('350\tZJetMax\t\t- Jet height [m] (used only for JET velocity profile, valid 70-490 m)\n')
            f.write('"{:.3f}"\tPLExp\t\t- Power law exponent [-] (or "default")\n'.format(params.PLexp))
            f.write('"default"\tZ0\t\t- Surface roughness length [m] (or "default")\n')
            f.write('\n')
            f.write('--------Non-IEC Meteorological Boundary Conditions------------\n')
            f.write('"default"\tLatitude\t\t- Site latitude [degrees] (or "default")\n')
            f.write('0.05\tRICH_NO\t\t- Gradient Richardson number [-]\n')
            f.write('"default"\tUStar\t\t- Friction or shear velocity [m/s] (or "default")\n')
            f.write('"default"\tZI\t\t- Mixing layer depth [m] (or "default")\n')
            f.write('"default"\tPC_UW\t\t- Hub mean u\'w\' Reynolds stress [m^2/s^2] (or "default" or "none")\n')
            f.write('"default"\tPC_UV\t\t- Hub mean u\'v\' Reynolds stress [m^2/s^2] (or "default" or "none")\n')
            f.write('"default"\tPC_VW\t\t- Hub mean v\'w\' Reynolds stress [m^2/s^2] (or "default" or "none")\n')
            f.write('\n')
            f.write('--------Spatial Coherence Parameters----------------------------\n')
            f.write('"IEC"\tSCMod1\t\t- u-component coherence model ("GENERAL","IEC","API","NONE", or "default")\n')
            f.write('"IEC"\tSCMod2\t\t- v-component coherence model ("GENERAL","IEC","NONE", or "default")\n')
            f.write('"IEC"\tSCMod3\t\t- w-component coherence model ("GENERAL","IEC","NONE", or "default")\n')
            f.write('"12.0 0.000659"\tInCDec1\t- u-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n')
            f.write('"12.0 0.000659"\tInCDec2\t- v-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n')
            f.write('"12.0 0.000659"\tInCDec3\t- w-component coherence parameters for general or IEC models [-, m^-1] (e.g. "10.0  0.3e-3" in quotes) (or "default")\n')
            f.write('"0.0"\tCohExp\t\t- Coherence exponent for general model [-] (or "default")\n')
            f.write('\n')
            f.write('--------Coherent Turbulence Scaling Parameters-------------------\n')
            f.write('".\EventData"\tCTEventPath\t\t- Name of the path where event data files are located\n')
            f.write('"random"\tCTEventFile\t\t- Type of event files ("LES", "DNS", or "RANDOM")\n')
            f.write('true\tRandomize\t\t- Randomize the disturbance scale and locations? (true/false)\n')
            f.write('1\tDistScl\t\t- Disturbance scale [-] (ratio of event dataset height to rotor disk). (Ignored when Randomize = true.)\n')
            f.write('0.5\tCTLy\t\t- Fractional location of tower centerline from right [-] (looking downwind) to left side of the dataset. (Ignored when Randomize = true.)\n')
            f.write('0.5\tCTLz\t\t- Fractional location of hub height from the bottom of the dataset. [-] (Ignored when Randomize = true.)\n')
            f.write('30\tCTStartTime\t\t- Minimum start time for coherent structures in RootName.cts [seconds]\n')
            f.write('\n')
            f.write('====================================================\n')
            f.write('! NOTE: Do not add or remove any lines in this file!\n')
            f.write('====================================================\n')

    else:
        print('Modifying {0} to be {1}'.format(fileIn,fileOut))

        NewPars = [int(params.nz), int(params.ny), int(params.dt), format(params.HubHt,'.2f'), format(params.Height,'.2f'), format(params.Width,'.2f'), format(params.TI,'.2f'), format(params.RefHt,'.2f'), format(params.URef,'.2f'), int(params.PLexp)]
        ModVars = ['NumGrid_Z','NumGrid_Y','TimeStep','HubHt','GridHeight','GridWidth','IECturb','RefHt','URef','PLExp']
        wt=0
        with open(fileOut, 'w+') as new_file:
            with open(fileIn) as old_file:
                for line in old_file.readlines():
                    newline = line
                    for index,tmpVar in enumerate(ModVars):
                        if tmpVar in line:
                            newline = str(NewPars[index])+'\t!!Orig is:  '+line
                    if '.fst' in line:
                        newline =str('{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{}_WT{:d}.fst"\t{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}\n'.format(params.x[wt],params.y[wt],params.z[wt],tpath,wt+1,params.X0_High[wt],params.Y0_High[wt],params.Z0_High,params.dX_High,params.dY_High,params.dZ_High))
                        wt+=1
                    new_file.write(newline)
