# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 08:54:48 2017

@author: pdoubraw

Edited by Kelsey Shaler
"""
import os, glob, struct
import numpy as np
from scipy import fftpack 
import matplotlib.pyplot as plt

class stochasticTurbulence:
    
    def __init__(self,D,prefix='prefix'):
        """
        Instantiate the object. 
        
        Parameters
        ----------
        prefix : string,
            can be used to read existing files or generate new ones.
        """
        self.prefix = prefix         
        self.D = D
        self.R = D/2.0
        
    def readBTS(self,pathToBinary, zHub):
        """
        Read TurbSim Binary FF.
        
        Parameters
        ----------
        pathToBinaries : string,
            where to find binary outputs        
        """      
        fPath = os.path.join(pathToBinary,"{0}.bts".format(self.prefix))
        
        filePath = glob.glob(fPath)        
        
        if len(filePath)==0:
            raise Exception("Could not find file at {0}.".format(fPath))        

        print('Opening file {0}...'.format(filePath[0]))
        self.filePath = filePath
        self.fileDir  = pathToBinary
        components = ['u','v','w']

        with open(filePath[0], mode='rb') as file:            
            fileContent = file.read()

            self.nZ, self.nY, self.nTower, self.nTimeSteps = \
                                    struct.unpack('i'*4,fileContent[2:18])
            self.dZ, self.dY, self.dT, self.uHub, dummy, self.zBot = \
                                    struct.unpack('f'*6,fileContent[18:42])
            self.nSeconds = self.nTimeSteps * self.dT

            self.zHub = zHub ## KS -- adding b/c HubHeight specified in TurbSim file is NOT the real HH
            vSlope = {} ; vIntercept = {} ; byteStart = 42    
            for component in components:
                vSlope[component], vIntercept[component] = struct.unpack('f'*2,
                          fileContent[byteStart:byteStart+8])
                byteStart += 8
            
            nChar = struct.unpack('i',fileContent[byteStart:byteStart+4])
            byteStart += 4
            
            vNumber = ""
            for i in range(int(nChar[0])):        
                vNumber += str(chr(struct.unpack('B',fileContent[byteStart:byteStart+1])[0]))
                byteStart += 1
            self.info = vNumber
               
            data = np.zeros((3,self.nY,self.nZ,self.nTimeSteps),float)
            for it in range(self.nTimeSteps):
                for iz in range(self.nZ):
                    for iy in range(self.nY):
                        for iComponent,component in enumerate(components):
                            data[iComponent,iy,iz,it] = struct.unpack('h',fileContent[byteStart:byteStart+2])[0]                        
                            byteStart += 2
                        
            for iComponent,component in enumerate(components):                    
                data[iComponent,:,:,:] = (data[iComponent,:,:,:] - vIntercept[component])/vSlope[component]
                   
            data = np.einsum('ljki->lijk',data)
            self.u = data[0,:,:,:]
            self.v = data[1,:,:,:]
            self.w = data[2,:,:,:]
            self.U = np.sqrt(self.u**2+self.v**2)
                
            if self.nTower>0:
                dataTower = np.zeros((self.nTower,self.nZ,self.nTimeSteps),float)
                for iz in range(self.nTower):
                    for iComponent,component in enumerate(components):
                        dataTower[iComponent,iz,it] = struct.unpack('h',fileContent[byteStart:byteStart+2])[0]                        
                        byteStart += 2
                for iComponent,component in enumerate(components):
                    dataTower[iComponent,:,:] = (dataTower[iComponent,:,:] - vIntercept[component])/vSlope[component]
            
        self.generateGrid()
        
    def generateGrid(self):
        """
        Generates mesh attributes.
        """
        self.y = np.array([ 0 + i*self.dY for i in range(self.nY) ])
        self.z = np.array([ self.zBot + i*self.dZ for i in range(self.nZ) ])
        [self.Y,self.Z] = np.meshgrid(self.y,self.z)
        self.kHub = self.z2k(self.zHub)
        self.yHub = np.mean(self.y)
        self.jHub = int(self.nY/2)
        self.tiHub = self.TI(j=self.jHub,k=self.kHub)
        self.getRotorPoints()
        
    def y2j(self,y):
        """
        Computes j grid index for a given y.
        """
        return np.argmin(np.abs(self.y-y))            
    
    def z2k(self,z):
        """
        Computes k grid index for a given z.
        """
        return np.argmin(np.abs(self.z-z))
  
    def TI(self,y=None,z=None,j=None,k=None):       
        """
        If no argument is given, compute TI over entire grid and return array of size (nY,nZ). Else, compute TI at the specified point.
        
        Parameters
        ----------
        y : float,
            cross-stream position [m]
        z : float,
            vertical position AGL [m]
        j : int,
            grid index along cross-stream
        k : int,
            grid index along vertical
        """
        if ((y==None) & (j==None)):
            return np.std(self.U,axis=0) / np.mean(self.U,axis=0)
        if ((y==None) & (j!=None)):
            return (np.std(self.U[:,j,k])/np.mean(self.U[:,j,k]))
        if ((y!=None) & (j==None)):
            uSeries = self.U[:,self.y2j(y),self.z2k(z)]
            return np.std(uSeries)/np.mean(uSeries)
    
    def visualize(self,component='U',time=0):
        """
        Quick peak at the data for a given component, at a specific time.
        """
        data    = getattr(self,component)[time,:,:]
        plt.figure() ; 
        plt.imshow(data) ; 
        plt.colorbar()   
        plt.show()
     
    def spectrum(self,component='u',y=None,z=None):
        """
        Calculate spectrum of a specific component, given time series at ~ hub.
        
        Parameters
        ----------
        component : string,
            which component to use
        y : float,
            y coordinate [m] of specific location
        z : float,
            z coordinate [m] of specific location

        """
        if y==None:
            k = self.kHub
            j = self.jHub
        data    = getattr(self,component)      
        data    = data[:,j,k]
        N       = data.size
        freqs   = fftpack.fftfreq(N,self.dT)[1:N/2]
        psd     = (np.abs(fftpack.fft(data,N)[1:N/2]))**2
        return [freqs, psd]        

    def getRotorPoints(self):
        """
        In the square y-z slice, return which points are at the edge of the rotor in the horizontal and vertical directions.
        
        Returns
        -------
        jLeft : int,
            index for grid point that matches the left side of the rotor (when looking towards upstream)
        jRight : int,
            index for grid point that matches the right side of the rotor (when looking towards upstream)
        kBot : int,
            index for grid point that matches the bottom of the rotor
        kTop : int,
            index for grid point that matches the top of the rotor
        """
        self.zBotRotor      = self.zHub - self.R
        self.zTopRotor      = self.zHub + self.R
        self.yLeftRotor     = self.yHub - self.R
        self.yRightRotor    = self.yHub + self.R        
        self.jLeftRotor  = self.y2j(self.yLeftRotor)
        self.jRightRotor = self.y2j(self.yRightRotor)
        self.kBotRotor   = self.z2k(self.zBotRotor)
        self.kTopRotor   = self.z2k(self.zTopRotor)
        
