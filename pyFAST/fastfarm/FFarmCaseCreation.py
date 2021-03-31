# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 2020

@author: kshaler
"""
import os, glob, struct
import numpy as np

class FFarmCaseCreation:
    
    def __init__(self,nTurbs,prefix='FFarm'):
        """
        Instantiate the object. 
        
        Parameters
        ----------
        prefix  :   string,
                   prefix for FAST.Farm simulations
        nTurbs:   integer,
                   number of turbines in the simulations
        """
        self.prefix = prefix
        self.nTurbs = nTurbs
        
    def Turb(self, D, HubHt, tpath, cmax=5, fmax=5):
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
        
        self.D = D
        self.HubHt = HubHt
        self.cmax = cmax
        self.fmax = fmax
        self.tpath = tpath
        
    def turbLocs(self,x,y,z):
        """
        Specify turbine locations
        
        Parameters
        ----------
        x, y, z   :   float,
               x-, y-, and z-location of turbine, respectively
        """
        self.x = x
        self.y = y
        self.z = z
        
    def discretization(self,Vhub,Cmeander=1.9):
        
        self.dX_High_desired = self.cmax
        self.dX_Low_desired = Cmeander*self.D*Vhub/150.0

    def highResDomain(self,ifdata,meanU,dx_des,extent_X,extent_Y):
        
        self.dT_High = ifdata.dT
        
        X0_des = {wt:{} for wt in range(self.nTurbs)}
        Y0_des = {wt:{} for wt in range(self.nTurbs)}
        self.X0_High = {wt:{} for wt in range(self.nTurbs)}
        self.Y0_High = {wt:{} for wt in range(self.nTurbs)}
        
        # high-box extent in x and y [D]
        self.Xdist = extent_X*self.D
        self.Ydist = extent_Y*self.D
        
        X0_rel = self.Xdist/2.0
        Y0_rel = self.Ydist/2.0

        for wt in range(self.nTurbs):
            X0_des[wt] = self.x[wt]-X0_rel#*self.D
            Y0_des[wt] = self.y[wt]-Y0_rel#*self.D)/2.0
        
        self.dY_High = ifdata.dY
        self.dZ_High = ifdata.dZ
        self.Z0_Low = ifdata.zBot
        self.Z0_High = ifdata.zBot
        
        self.Width = ifdata.dY*(ifdata.nY-1)
        height = ifdata.dZ*(ifdata.nZ-1)
        
        effSimLength = ifdata.nSeconds+self.Width/meanU
        
        length = effSimLength*meanU
        nx = round(effSimLength/ifdata.dT)
        dx_TS = length/(nx-1)
        
        self.dX_High = round(dx_des/dx_TS)*dx_TS
        
        self.Zdist_High = height
        
        self.nX_High = round(self.Xdist/self.dX_High)+1
        self.nY_High = round(self.Ydist/self.dY_High)+1
        self.nZ_High = round(self.Zdist_High/self.dZ_High)+1

        for wt in range(self.nTurbs):
            self.X0_High[wt] = round(X0_des[wt]/self.dX_High)*self.dX_High
            self.Y0_High[wt] = round(Y0_des[wt]/self.dY_High)*self.dY_High
            
    def lowResDomain(self,dx_des,Vhub,Cmeander=1.9):
        
        dt_des = Cmeander*self.D/(10.0*Vhub)
        self.dT_Low = round(dt_des/self.dT_High)*self.dT_High

        dy_des = dx_des
        dz_des = dx_des

        self.X0_Low = min(self.x)-3*self.D
        self.Y0_Low = -self.Width/2

        self.dX_Low = round(dx_des/self.dX_High)*self.dX_High
        self.dY_Low = round(dy_des/self.dY_High)*self.dY_High
        self.dZ_Low = round(dz_des/self.dZ_High)*self.dZ_High

        xdist = max(self.x)+8.0*self.D-self.X0_Low

        self.nX_Low = round(xdist/self.dX_Low)+1
        self.nY_Low = round(self.Width/self.dY_Low)+1
        self.nZ_Low = round(self.Zdist_High/self.dZ_Low)+1
        
def WriteFFarmFile(fileIn,fileOut,params,NewFile=True):
    """ Write a FAST.Farm primary input file, 

    """

    if NewFile == True:
        print('Writing a new {0} file from scratch'.format(fileOut))
        # --- Writing FFarm input file from scratch
        with open(fileOut, 'w') as f:
            f.write('FAST.Farm v1.00.* INPUT FILE\n')
            f.write('Sample FAST.Farm input file\n')
            f.write('--- SIMULATION CONTROL ---\n')
            f.write('False\tEcho               Echo input data to <RootName>.ech? (flag)\n')
            f.write('FATAL\tAbortLevel         Error level when simulation should abort (string) {"WARNING", "SEVERE", "FATAL"}\n')
            f.write('2000.0\tTMax               Total run time (s) [>=0.0]\n')
            f.write('False\tUseSC              Use a super controller? (flag)\n')
            f.write('2\tMod_AmbWind        Ambient wind model (-) (switch) {1: high-fidelity precursor in VTK format, 2: InflowWind module}\n')
            f.write('--- SUPER CONTROLLER --- [used only for UseSC=True]\n')
            f.write('"SC_DLL.dll"       SC_FileName        Name/location of the dynamic library {.dll [Windows] or .so [Linux]} containing the Super Controller algorithms (quoated string)\n')
            f.write('--- AMBIENT WIND ---\n')
            f.write('{:.1f}\tDT                 Time step for low -resolution wind data input files; will be used as the global FAST.Farm time step (s) [>0.0]\n'.format(params.dT_Low))
            f.write('{:.6f}\tDT_High            Time step for high-resolution wind data input files (s) [>0.0]\n'.format(params.dT_High))
            f.write('"/projects/isda/kshaler/WESC/freestream/08ms"          WindFilePath       Path name to wind data files from precursor (string)\n')
            f.write('False              ChkWndFiles        Check all the ambient wind files for data consistency? (flag)\n')
            f.write('--- AMBIENT WIND: INFLOWWIND MODULE --- [used only for Mod_AmbWind=2]\n')
            f.write('{:.3f}\t\tDT                 Time step for low -resolution wind data interpolation; will be used as the global FAST.Farm time step (s) [>0.0]\n'.format(params.dT_Low))
            f.write('{:.6f}\t\tDT_High            Time step for high-resolution wind data interpolation (s) [>0.0]\n'.format(params.dT_High))
            f.write('{:.0f}\t\tNX_Low             Number  of low -resolution spatial nodes in X direction for wind data interpolation (-) [>=2]\n'.format(params.nX_Low))
            f.write('{:.0f}\t\tNY_Low             Number  of low -resolution spatial nodes in Y direction for wind data interpolation (-) [>=2]\n'.format(params.nY_Low))
            f.write('{:.0f}\t\tNZ_Low             Number  of low -resolution spatial nodes in Z direction for wind data interpolation (-) [>=2]\n'.format(params.nZ_Low))
            f.write('{:.3f}\tX0_Low             Origin  of low -resolution spatial nodes in X direction for wind data interpolation (m)\n'.format(params.X0_Low))
            f.write('{:.3f}\tY0_Low             Origin  of low -resolution spatial nodes in Y direction for wind data interpolation (m)\n'.format(params.Y0_Low))
            f.write('{:.3f}\t\tZ0_Low             Origin  of low -resolution spatial nodes in Z direction for wind data interpolation (m)\n'.format(params.Z0_Low))
            f.write('{:.3f}\t\tdX_Low             Spacing of low -resolution spatial nodes in X direction for wind data interpolation (m) [>0.0]\n'.format(params.dX_Low))
            f.write('{:.3f}\t\tdY_Low             Spacing of low -resolution spatial nodes in Y direction for wind data interpolation (m) [>0.0]\n'.format(params.dY_Low))
            f.write('{:.3f}\t\tdZ_Low             Spacing of low -resolution spatial nodes in Z direction for wind data interpolation (m) [>0.0]\n'.format(params.dZ_Low))
            f.write('{:.0f}\t\tNX_High            Number  of high-resolution spatial nodes in X direction for wind data interpolation (-) [>=2]\n'.format(params.nX_High))
            f.write('{:.0f}\t\tNY_High            Number  of high-resolution spatial nodes in Y direction for wind data interpolation (-) [>=2]\n'.format(params.nY_High))
            f.write('{:.0f}\t\tNZ_High            Number  of high-resolution spatial nodes in Z direction for wind data interpolation (-) [>=2]\n'.format(params.nZ_High))
            f.write('"InflowWind.dat"   InflowFile         Name of file containing InflowWind module input parameters (quoted string)\n')
            f.write('--- WIND TURBINES ---\n')
            f.write('{:.0f}                  NumTurbines        Number of wind turbines (-) [>=1]                          [last 6 columns used only for Mod_AmbWind=2]\n'.format(params.nTurbs))                
            f.write('WT_X   WT_Y   WT_Z   WT_FASTInFile                                                               X0_High  Y0_High  Z0_High  dX_High  dY_High  dZ_High\n')
            f.write('(m)    (m)    (m)    (string)                                                                    (m)      (m)      (m)      (m)      (m)      (m)\n')
            for wt in range(params.nTurbs):
                f.write('{:.3f}\t{:.3f}\t{:.3f}\t{}_WT{:d}.fst"\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(params.x[wt],params.y[wt],params.z[wt],params.tpath,wt+1,params.X0_High[wt],params.Y0_High[wt],params.Z0_High,params.dX_High,params.dY_High,params.dZ_High))
            f.write('--- WAKE DYNAMICS ---\n')
            f.write('3.0\t\tdr                 Radial increment of radial finite-difference grid (m) [>0.0]\n')
            f.write('50\t\tNumRadii           Number of radii in the radial finite-difference grid (-) [>=2]\n')
            f.write('44\t\tNumPlanes          Number of wake planes (-) [>=2]\n')
            f.write('DEFAULT\t\tf_c                Cut-off (corner) frequency of the low-pass time-filter for the wake advection, deflection, and meandering model (Hz) [>0.0] or DEFAULT [DEFAULT=0.0007]\n')
            f.write('DEFAULT\t\tC_HWkDfl_O         Calibrated parameter in the correction for wake deflection defining the horizontal offset at the rotor                                               (m    ) or DEFAULT [DEFAULT= 0.0  ]\n')
            f.write('DEFAULT\t\tC_HWkDfl_OY        Calibrated parameter in the correction for wake deflection defining the horizontal offset at the rotor scaled with                         yaw error (m/deg) or DEFAULT [DEFAULT= 0.3  ]\n')
            f.write('DEFAULT\t\tC_HWkDfl_x         Calibrated parameter in the correction for wake deflection defining the horizontal offset              scaled with downstream distance               (-    ) or DEFAULT [DEFAULT= 0.0  ]\n')
            f.write('DEFAULT\t\tC_HWkDfl_xY        Calibrated parameter in the correction for wake deflection defining the horizontal offset              scaled with downstream distance and yaw error (1/deg) or DEFAULT [DEFAULT=-0.004]\n')
            f.write('DEFAULT\t\tC_NearWake         Calibrated parameter for the near-wake correction (-) [>1.0] or DEFAULT [DEFAULT=1.8]\n')
            f.write('DEFAULT\t\tk_vAmb             Calibrated parameter for the influence of ambient turbulence in the eddy viscosity (-) [>=0.0] or DEFAULT [DEFAULT=0.05 ]\n')
            f.write('DEFAULT\t\tk_vShr             Calibrated parameter for the influence of the shear layer    in the eddy viscosity (-) [>=0.0] or DEFAULT [DEFAULT=0.016]\n')
            f.write('DEFAULT\t\tC_vAmb_DMin        Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the minimum and exponential regions (-) [>=0.0          ] or DEFAULT [DEFAULT= 0.0 ]\n')
            f.write('DEFAULT\t\tC_vAmb_DMax        Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the transitional diameter fraction between the exponential and maximum regions (-) [> C_vAmb_DMin  ] or DEFAULT [DEFAULT= 1.0 ]\n')
            f.write('DEFAULT\t\tC_vAmb_FMin        Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the value in the minimum region                                                (-) [>=0.0 and <=1.0] or DEFAULT [DEFAULT= 1.0 ]\n')
            f.write('DEFAULT\t\tC_vAmb_Exp         Calibrated parameter in the eddy viscosity filter function for ambient turbulence defining the exponent in the exponential region                                         (-) [> 0.0          ] or DEFAULT [DEFAULT= 0.01]\n')
            f.write('DEFAULT\t\tC_vShr_DMin        Calibrated parameter in the eddy viscosity filter function for the shear layer    defining the transitional diameter fraction between the minimum and exponential regions (-) [>=0.0          ] or DEFAULT [DEFAULT= 3.0 ]\n')
            f.write('DEFAULT\t\tC_vShr_DMax        Calibrated parameter in the eddy viscosity filter function for the shear layer    defining the transitional diameter fraction between the exponential and maximum regions (-) [> C_vShr_DMin  ] or DEFAULT [DEFAULT=25.0 ]\n')
            f.write('DEFAULT\t\tC_vShr_FMin        Calibrated parameter in the eddy viscosity filter function for the shear layer    defining the value in the minimum region                                                (-) [>=0.0 and <=1.0] or DEFAULT [DEFAULT= 0.2 ]\n')
            f.write('DEFAULT\t\tC_vShr_Exp         Calibrated parameter in the eddy viscosity filter function for the shear layer    defining the exponent in the exponential region                                         (-) [> 0.0          ] or DEFAULT [DEFAULT= 0.1 ]\n')
            f.write('DEFAULT\t\tMod_WakeDiam       Wake diameter calculation model (-) (switch) {1: rotor diameter, 2: velocity based, 3: mass-flux based, 4: momentum-flux based} or DEFAULT [DEFAULT=1]\n')
            f.write('DEFAULT\t\tC_WakeDiam         Calibrated parameter for wake diameter calculation (-) [>0.0 and <0.99] or DEFAULT [DEFAULT=0.95] [unused for Mod_WakeDiam=1]\n')
            f.write('DEFAULT\t\tMod_Meander        Spatial filter model for wake meandering (-) (switch) {1: uniform, 2: truncated jinc, 3: windowed jinc} or DEFAULT [DEFAULT=3]\n')
            f.write('DEFAULT\t\tC_Meander          Calibrated parameter for wake meandering (-) [>=1.0] or DEFAULT [DEFAULT=1.9]\n')
            f.write('--- VISUALIZATION ---\n')
            f.write('False\t\tWrDisWind          Write low- and high-resolution disturbed wind data to <RootName>.Low.Dis.t<n>.vtk etc.? (flag)\n')
            f.write('1\t\tNOutDisWindXY      Number of XY planes for output of disturbed wind data across the low-resolution domain to <RootName>.Low.DisXY<n_out>.t<n>.vtk (-) [0 to 9]\n')
            f.write('80.0\t\tOutDisWindZ        Z coordinates of XY planes for output of disturbed wind data across the low-resolution domain (m) [1 to NOutDisWindXY] [unused for NOutDisWindXY=0]\n')
            f.write('0\t\tNOutDisWindYZ      Number of YZ planes for output of disturbed wind data across the low-resolution domain to <RootName>/Low.DisYZ<n_out>.t<n>.vtk (-) [0 to 9]\n')
            f.write('748.0, 1252.0, 1378.0, 1504.0, 1630.0, 1756.0, 1882.0, 2008.0   OutDisWindX        X coordinates of YZ planes for output of disturbed wind data across the low-resolution domain (m) [1 to NOutDisWindYZ] [unused for NOutDisWindYZ=0]\n')
            f.write('0\t\tNOutDisWindXZ      Number of XZ planes for output of disturbed wind data across the low-resolution domain to <RootName>/Low.DisXZ<n_out>.t<n>.vtk (-) [0 to 9]\n')
            f.write('0.0\tO\tutDisWindY        Y coordinates of XZ planes for output of disturbed wind data across the low-resolution domain (m) [1 to NOutDisWindXZ] [unused for NOutDisWindXZ=0]\n')
            f.write('2.0\t\tWrDisDT            Time step for disturbed wind visualization output (s) [>0.0] or DEFAULT [DEFAULT=DT] [unused for WrDisWind=False and NOutDisWindXY=NOutDisWindYZ=NOutDisWindXZ=0]\n')
            f.write('--- OUTPUT ---\n')
            f.write('True\t\tSumPrint           Print summary data to <RootName>.sum? (flag)\n')
            f.write('99999.9\t\tChkptTime          Amount of time between creating checkpoint files for potential restart (s) [>0.0]\n')
            f.write('200.0\t\tTStart             Time to begin tabular output (s) [>=0.0]\n')
            f.write('1\t\tOutFileFmt         Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}\n')
            f.write('True\t\tTabDelim           Use tab delimiters in text tabular output file? (flag) {uses spaces if False}\n')
            f.write('"ES10.3E2"\tOutFmt             Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)\n')
            f.write('20\t\tNOutRadii          Number of radial nodes         for wake output for an individual rotor (-) [0 to 20]\n')
            f.write('0, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 16, 17, 18, 19, 21, 24, 28, 33, 39  OutRadii           List of radial nodes         for wake output for an individual rotor (-) [1 to NOutRadii] [unused for NOutRadii=0]\n')
            f.write('7\t\tNOutDist           Number of downstream distances for wake output for an individual rotor (-) [0 to 9 ]\n')
            f.write('252.0, 378.0, 504.0, 630.0, 756.0, 882.0, 1008.0          OutDist            List of downstream distances for wake output for an individual rotor (m) [1 to NOutDist ] [unused for NOutDist =0]\n')
            f.write('1\t\tNWindVel           Number of points for wind output (-) [0 to 9]\n')
            f.write('0.0\t\tWindVelX           List of coordinates in the X direction for wind output (m) [1 to NWindVel] [unused for NWindVel=0]\n')
            f.write('0.0\t\tWindVelY           List of coordinates in the Y direction for wind output (m) [1 to NWindVel] [unused for NWindVel=0]\n')
            f.write('80.0\t\tWindVelZ           List of coordinates in the Z direction for wind output (m) [1 to NWindVel] [unused for NWindVel=0]\n')
            f.write('\tOutList            The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels (quoted string)\n')
            f.write('"RtAxsXT1     , RtAxsYT1     , RtAxsZT1"\n')
            f.write('"RtPosXT1     , RtPosYT1     , RtPosZT1"\n')
            f.write('"YawErrT1"\n')
            f.write('"TIAmbT1"\n')
            f.write('"CtT1N01      , CtT1N02      , CtT1N03      , CtT1N04      , CtT1N05      , CtT1N06      , CtT1N07      , CtT1N08      , CtT1N09      , CtT1N10      , CtT1N11      , CtT1N12      , CtT1N13      , CtT1N14      , CtT1N15      , CtT1N16      , CtT1N17      , CtT1N18      , CtT1N19      , CtT1N20"\n')
            f.write('"WkAxsXT1D1   , WkAxsXT1D2   , WkAxsXT1D3   , WkAxsXT1D4   , WkAxsXT1D5   , WkAxsXT1D6   , WkAxsXT1D7"\n')
            f.write('"WkAxsYT1D1   , WkAxsYT1D2   , WkAxsYT1D3   , WkAxsYT1D4   , WkAxsYT1D5   , WkAxsYT1D6   , WkAxsYT1D7"\n')
            f.write('"WkAxsZT1D1   , WkAxsZT1D2   , WkAxsZT1D3   , WkAxsZT1D4   , WkAxsZT1D5   , WkAxsZT1D6   , WkAxsZT1D7"\n')
            f.write('"WkPosXT1D1   , WkPosXT1D2   , WkPosXT1D3   , WkPosXT1D4   , WkPosXT1D5   , WkPosXT1D6   , WkPosXT1D7"\n')
            f.write('"WkPosYT1D1   , WkPosYT1D2   , WkPosYT1D3   , WkPosYT1D4   , WkPosYT1D5   , WkPosYT1D6   , WkPosYT1D7"\n')
            f.write('"WkPosZT1D1   , WkPosZT1D2   , WkPosZT1D3   , WkPosZT1D4   , WkPosZT1D5   , WkPosZT1D6   , WkPosZT1D7"\n')
            f.write('"WkDfVxT1N01D1, WkDfVxT1N02D1, WkDfVxT1N03D1, WkDfVxT1N04D1, WkDfVxT1N05D1, WkDfVxT1N06D1, WkDfVxT1N07D1, WkDfVxT1N08D1, WkDfVxT1N09D1, WkDfVxT1N10D1, WkDfVxT1N11D1, WkDfVxT1N12D1, WkDfVxT1N13D1, WkDfVxT1N14D1, WkDfVxT1N15D1, WkDfVxT1N16D1, WkDfVxT1N17D1, WkDfVxT1N18D1, WkDfVxT1N19D1, WkDfVxT1N20D1"\n')
            f.write('"WkDfVxT1N01D2, WkDfVxT1N02D2, WkDfVxT1N03D2, WkDfVxT1N04D2, WkDfVxT1N05D2, WkDfVxT1N06D2, WkDfVxT1N07D2, WkDfVxT1N08D2, WkDfVxT1N09D2, WkDfVxT1N10D2, WkDfVxT1N11D2, WkDfVxT1N12D2, WkDfVxT1N13D2, WkDfVxT1N14D2, WkDfVxT1N15D2, WkDfVxT1N16D2, WkDfVxT1N17D2, WkDfVxT1N18D2, WkDfVxT1N19D2, WkDfVxT1N20D2"\n')
            f.write('"WkDfVxT1N01D3, WkDfVxT1N02D3, WkDfVxT1N03D3, WkDfVxT1N04D3, WkDfVxT1N05D3, WkDfVxT1N06D3, WkDfVxT1N07D3, WkDfVxT1N08D3, WkDfVxT1N09D3, WkDfVxT1N10D3, WkDfVxT1N11D3, WkDfVxT1N12D3, WkDfVxT1N13D3, WkDfVxT1N14D3, WkDfVxT1N15D3, WkDfVxT1N16D3, WkDfVxT1N17D3, WkDfVxT1N18D3, WkDfVxT1N19D3, WkDfVxT1N20D3"\n')
            f.write('"WkDfVxT1N01D4, WkDfVxT1N02D4, WkDfVxT1N03D4, WkDfVxT1N04D4, WkDfVxT1N05D4, WkDfVxT1N06D4, WkDfVxT1N07D4, WkDfVxT1N08D4, WkDfVxT1N09D4, WkDfVxT1N10D4, WkDfVxT1N11D4, WkDfVxT1N12D4, WkDfVxT1N13D4, WkDfVxT1N14D4, WkDfVxT1N15D4, WkDfVxT1N16D4, WkDfVxT1N17D4, WkDfVxT1N18D4, WkDfVxT1N19D4, WkDfVxT1N20D4"\n')
            f.write('"WkDfVxT1N01D5, WkDfVxT1N02D5, WkDfVxT1N03D5, WkDfVxT1N04D5, WkDfVxT1N05D5, WkDfVxT1N06D5, WkDfVxT1N07D5, WkDfVxT1N08D5, WkDfVxT1N09D5, WkDfVxT1N10D5, WkDfVxT1N11D5, WkDfVxT1N12D5, WkDfVxT1N13D5, WkDfVxT1N14D5, WkDfVxT1N15D5, WkDfVxT1N16D5, WkDfVxT1N17D5, WkDfVxT1N18D5, WkDfVxT1N19D5, WkDfVxT1N20D5"\n')
            f.write('"WkDfVxT1N01D6, WkDfVxT1N02D6, WkDfVxT1N03D6, WkDfVxT1N04D6, WkDfVxT1N05D6, WkDfVxT1N06D6, WkDfVxT1N07D6, WkDfVxT1N08D6, WkDfVxT1N09D6, WkDfVxT1N10D6, WkDfVxT1N11D6, WkDfVxT1N12D6, WkDfVxT1N13D6, WkDfVxT1N14D6, WkDfVxT1N15D6, WkDfVxT1N16D6, WkDfVxT1N17D6, WkDfVxT1N18D6, WkDfVxT1N19D6, WkDfVxT1N20D6"\n')
            f.write('"WkDfVxT1N01D7, WkDfVxT1N02D7, WkDfVxT1N03D7, WkDfVxT1N04D7, WkDfVxT1N05D7, WkDfVxT1N06D7, WkDfVxT1N07D7, WkDfVxT1N08D7, WkDfVxT1N09D7, WkDfVxT1N10D7, WkDfVxT1N11D7, WkDfVxT1N12D7, WkDfVxT1N13D7, WkDfVxT1N14D7, WkDfVxT1N15D7, WkDfVxT1N16D7, WkDfVxT1N17D7, WkDfVxT1N18D7, WkDfVxT1N19D7, WkDfVxT1N20D7"\n')
            f.write('"WkDfVrT1N01D1, WkDfVrT1N02D1, WkDfVrT1N03D1, WkDfVrT1N04D1, WkDfVrT1N05D1, WkDfVrT1N06D1, WkDfVrT1N07D1, WkDfVrT1N08D1, WkDfVrT1N09D1, WkDfVrT1N10D1, WkDfVrT1N11D1, WkDfVrT1N12D1, WkDfVrT1N13D1, WkDfVrT1N14D1, WkDfVrT1N15D1, WkDfVrT1N16D1, WkDfVrT1N17D1, WkDfVrT1N18D1, WkDfVrT1N19D1, WkDfVrT1N20D1"\n')
            f.write('"WkDfVrT1N01D2, WkDfVrT1N02D2, WkDfVrT1N03D2, WkDfVrT1N04D2, WkDfVrT1N05D2, WkDfVrT1N06D2, WkDfVrT1N07D2, WkDfVrT1N08D2, WkDfVrT1N09D2, WkDfVrT1N10D2, WkDfVrT1N11D2, WkDfVrT1N12D2, WkDfVrT1N13D2, WkDfVrT1N14D2, WkDfVrT1N15D2, WkDfVrT1N16D2, WkDfVrT1N17D2, WkDfVrT1N18D2, WkDfVrT1N19D2, WkDfVrT1N20D2"\n')
            f.write('"WkDfVrT1N01D3, WkDfVrT1N02D3, WkDfVrT1N03D3, WkDfVrT1N04D3, WkDfVrT1N05D3, WkDfVrT1N06D3, WkDfVrT1N07D3, WkDfVrT1N08D3, WkDfVrT1N09D3, WkDfVrT1N10D3, WkDfVrT1N11D3, WkDfVrT1N12D3, WkDfVrT1N13D3, WkDfVrT1N14D3, WkDfVrT1N15D3, WkDfVrT1N16D3, WkDfVrT1N17D3, WkDfVrT1N18D3, WkDfVrT1N19D3, WkDfVrT1N20D3"\n')
            f.write('"WkDfVrT1N01D4, WkDfVrT1N02D4, WkDfVrT1N03D4, WkDfVrT1N04D4, WkDfVrT1N05D4, WkDfVrT1N06D4, WkDfVrT1N07D4, WkDfVrT1N08D4, WkDfVrT1N09D4, WkDfVrT1N10D4, WkDfVrT1N11D4, WkDfVrT1N12D4, WkDfVrT1N13D4, WkDfVrT1N14D4, WkDfVrT1N15D4, WkDfVrT1N16D4, WkDfVrT1N17D4, WkDfVrT1N18D4, WkDfVrT1N19D4, WkDfVrT1N20D4"\n')
            f.write('"WkDfVrT1N01D5, WkDfVrT1N02D5, WkDfVrT1N03D5, WkDfVrT1N04D5, WkDfVrT1N05D5, WkDfVrT1N06D5, WkDfVrT1N07D5, WkDfVrT1N08D5, WkDfVrT1N09D5, WkDfVrT1N10D5, WkDfVrT1N11D5, WkDfVrT1N12D5, WkDfVrT1N13D5, WkDfVrT1N14D5, WkDfVrT1N15D5, WkDfVrT1N16D5, WkDfVrT1N17D5, WkDfVrT1N18D5, WkDfVrT1N19D5, WkDfVrT1N20D5"\n')
            f.write('"WkDfVrT1N01D6, WkDfVrT1N02D6, WkDfVrT1N03D6, WkDfVrT1N04D6, WkDfVrT1N05D6, WkDfVrT1N06D6, WkDfVrT1N07D6, WkDfVrT1N08D6, WkDfVrT1N09D6, WkDfVrT1N10D6, WkDfVrT1N11D6, WkDfVrT1N12D6, WkDfVrT1N13D6, WkDfVrT1N14D6, WkDfVrT1N15D6, WkDfVrT1N16D6, WkDfVrT1N17D6, WkDfVrT1N18D6, WkDfVrT1N19D6, WkDfVrT1N20D6"\n')
            f.write('"WkDfVrT1N01D7, WkDfVrT1N02D7, WkDfVrT1N03D7, WkDfVrT1N04D7, WkDfVrT1N05D7, WkDfVrT1N06D7, WkDfVrT1N07D7, WkDfVrT1N08D7, WkDfVrT1N09D7, WkDfVrT1N10D7, WkDfVrT1N11D7, WkDfVrT1N12D7, WkDfVrT1N13D7, WkDfVrT1N14D7, WkDfVrT1N15D7, WkDfVrT1N16D7, WkDfVrT1N17D7, WkDfVrT1N18D7, WkDfVrT1N19D7, WkDfVrT1N20D7"\n')
            f.write('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n')

    else:
        print('Modifying {0} to be {1}'.format(fileIn,fileOut))

        NewPars = [format(params.dT_Low,'.1f'), format(params.dT_High,'.6f'), int(params.nX_Low), int(params.nY_Low), int(params.nZ_Low), format(params.X0_Low,'.2f'), format(params.Y0_Low,'.2f'), format(params.Z0_Low,'.2f'), format(params.dX_Low,'.2f'), format(params.dY_Low,'.2f'), format(params.dZ_Low,'.2f'), int(params.nX_High), int(params.nY_High), int(params.nZ_High)]
        ModVars = ['DT','DT_High','NX_Low','NY_Low','NZ_Low','X0_Low','Y0_Low','Z0_Low','dX_Low','dY_Low','dZ_Low','NX_High','NY_High','NZ_High']
        wt=0
        with open(fileOut, 'w+') as new_file:
            with open(fileIn) as old_file:
                for line in old_file.readlines():
                    newline = line
                    for index,tmpVar in enumerate(ModVars):
                        if tmpVar in line:
                            newline = str(NewPars[index])+'\t!!Orig is:  '+line
                    if '.fst' in line:
                        newline =str('{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{}_WT{:d}.fst"\t{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}\t\t{:.3f}\n'.format(params.x[wt],params.y[wt],params.z[wt],params.tpath,wt+1,params.X0_High[wt],params.Y0_High[wt],params.Z0_High,params.dX_High,params.dY_High,params.dZ_High))
                        wt+=1
                    new_file.write(newline)