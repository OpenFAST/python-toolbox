# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:37:53 2023

@author: mahfouz
"""

import os, glob, struct,math
import numpy as np
import subprocess
import pdb

class MBCreation:
    
    def __init__(self, D, HubHt, Vhub, PLexp,  y, zbot=1.0, cmax=5.0, low_ext=None,fmax=1, 
                 Lm=None, gamma=None, alphaepsilon=None, seedvalue=None, time=None):
        """
        Instantiate the object. 
        
        Parameters
        ----------
        D        :   float,
                    rotor diameter (m)
        HubHt    :   float,
                    turbine hub height (m)
        Vhub	 :   float,
                    mean wind speed at hub height (m/s)
        PLexp    :   float,
                    power law exponent for shear (-)
        y        :   float,
                     y- location of turbine, respectively
        cmax     :   float,
                    maximum blade chord (m). If not specified, set to NREL 5MW value.
        low_ext  :  list of floats [xmin, xmax, ymin, ymax, zabovehub]
                    extents for the low-res box. All values should be positive  If not specified, resorts to
                    computations by the manual
        """

        # Perform some checks on the input
        if low_ext is not None and len(low_ext) != 5:
            raise ValueError('low_ext not defined properly. It should be [xmin, xmax, ymin, ymax, zabovehub]')

        # Set parameters for convenience
        self.low_ext            = low_ext
        self.Lm                 = Lm
        self.gamma              = gamma
        self.alphaepsilon       = alphaepsilon      
        self.seedvalue          = seedvalue
        self.time               = time
        self.Vhub               = Vhub
        self.cmax               = cmax
        self.HubHt              = HubHt
        self.D                  = D
        self.zbot               = zbot
        self.PLexp              = PLexp
        self.y                  = y
        # pdb.set_trace()
        
        
        self.ymin = min(self.y) - self.low_ext[2]*self.D
        self.ymax = max(self.y) + self.low_ext[3]*self.D
        self.dt_targeted         = 1.0/(2.0*fmax)
        
        
        
    def writeMannbatch(self, Mannpath, windfilename, outputpath):
        i=1  
        dT_Box=5
        box_width = self.ymax - self.ymin
        box_height= self.HubHt + 0.5* self.D + self.low_ext[4]*self.D
        # pdb.set_trace()
        while dT_Box > self.dt_targeted:
            dx=self.cmax-i
            dy=dx
            dz=dx
            ny=box_width/dy
            nz=box_height/dz
            # pdb.set_trace()
            ny=2**(math.ceil(math.log(ny, 2)))
            nz=2**(math.ceil(math.log(nz, 2)))
    
    
            RefHt=0.5*dz*(nz-1)+self.zbot 
            Uref=self.Vhub/((self.HubHt/RefHt)**self.PLexp)
            dT_Box      = dx/Uref
            nx=2**(math.ceil(math.log(Uref*self.time/dx, 2)))
            i=i+1
        # pdb.set_trace()
        if not os.path.exists(outputpath):
            os.makedirs(outputpath)
        filename=windfilename
        
        batch_file_content = "@echo off\n"        
        batch_file_content += "start \"\"  cmd /C " + Mannpath +" "+ os.path.join(outputpath, filename)  + " " + str(self.alphaepsilon) + " "+ str(self.Lm) + " "+ str(self.gamma) + ' '+ str(self.seedvalue )+ ' '+ str(nx)+ ' '+ str(ny)+ ' '+ str(nz) + ' '+ str(dx)+ ' '+ str(dy)+ ' '+ str(dz) +' false ' + "\n"
        batch_file_content += "echo Script execution completed\n"

        with open(os.path.join(outputpath, filename+'.bat'), "w") as batch_file:
            batch_file.write(batch_file_content)  
            
        subprocess.run([os.path.join(outputpath, filename+'.bat')],stdout=subprocess.PIPE)
