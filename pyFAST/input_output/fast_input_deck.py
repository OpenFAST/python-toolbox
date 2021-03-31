import os
import numpy as np
import re
import pandas as pd

from .fast_input_file import FASTInputFile

__all__  = ['FASTInputDeck']
# --------------------------------------------------------------------------------}
# --- Full FAST input deck
# --------------------------------------------------------------------------------{
class FASTInputDeck(dict):
    """Container for input files that make up a FAST input deck"""

    def __init__(self, fullFstPath, readlist=['all'], verbose=False):
        """Read FAST master file and read inputs for FAST modules

        INPUTS:
          - fullFstPath: 
          - readlist: list of module files to be read, or ['all'], modules are identified as follows:
                ['Fst','ED','AD','BD','BDbld','EDtwr','EDbld','ADbld','AF','IW','HD','SrvD','SD','MD']

        """
        self.filename = fullFstPath
        self.verbose  = verbose
        self.readlist = readlist
        if not type(self.readlist) is list:
            self.readlist=[readlist]
        if 'all' in self.readlist:
            self.readlist = ['Fst','ED','AD','BD','BDbld','EDtwr','EDbld','ADbld','AF','IW','HD','SrvD','SD','MD']
        else:
            self.readlist = ['Fst']+self.readlist

        self.inputfiles = {}

        # --- Harmonization with AeroElasticSE
        self.FAST_ver       = 'OPENFAST'
        self.path2dll       = None   # Path to dll file

        self.fst_vt={}
        self.fst_vt['description']       = ''            
        self.fst_vt['Fst']               = None
        self.fst_vt['ElastoDyn']         = None
        self.fst_vt['ElastoDynBlade']    = None
        self.fst_vt['ElastoDynTower']    = None
        self.fst_vt['InflowWind']        = None
        self.fst_vt['AeroDyn14']         = None
        self.fst_vt['AeroDyn15']         = None
        self.fst_vt['AeroDynBlade']      = None 
        self.fst_vt['AeroDynTower']      = None
        self.fst_vt['AeroDynPolar']      = None
        self.fst_vt['ServoDyn']          = None
        self.fst_vt['DISCON_in']         = None
        self.fst_vt['HydroDyn']          = None
        self.fst_vt['MoorDyn']           = None
        self.fst_vt['SubDyn']            = None
        self.fst_vt['MAP']               = None
        self.fst_vt['BeamDyn']           = None
        self.fst_vt['BeamDynBlade']      = None # Small change of interface
        self.fst_vt['af_data']           = [] # Small change of interface


        # Read all inputs files
        self.read()


    @property
    def FAST_InputFile(self):
        return os.path.basename(self.filename)   # FAST input file (ext=.fst)
    @property
    def FAST_directory(self):
        return os.path.dirname(self.filename)    # Path to fst directory files

    def read(self, filename=None):
        if filename is not None:
            self.filename = filename

        # Read OpenFAST files
        self.fst_vt['Fst'] = self._read(self.FAST_InputFile, 'Fst')
        if self.fst_vt['Fst'] is None:
            raise Exception('Error reading main file {}'.format(self.filename))
        keys = self.fst_vt['Fst'].keys()


        if 'NumTurbines' in keys:
            self.version='AD_driver'
        elif 'InterpOrder' in self.fst_vt['Fst'].keys():
            self.version='OF2'
        else:
            self.version='F7'


        if self.version=='AD_driver':
            # ---- AD Driver
            # InflowWind
            if self.fst_vt['Fst']['CompInflow']>0:
                self.fst_vt['InflowWind'] = self._read(self.fst_vt['Fst']['InflowFile'],'IW')
            self.fst_vt['AeroDyn15'] = self._read(self.fst_vt['Fst']['AeroFile'],'AD')
            if self.fst_vt['AeroDyn15'] is not None:
                # Blades
                bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']), self.fst_vt['AeroDyn15']['ADBlFile(1)'])
                self.fst_vt['AeroDynBlade'] = self._read(bld_file,'ADbld')

                #self.fst_vt['AeroDynBlade'] = []
                #for i in range(3):
                #    bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']), self.fst_vt[key]['ADBlFile({})'.format(i+1)])
                #    self.fst_vt['AeroDynBlade'].append(self._read(bld_file,'ADbld'))
                # Polars
                self.fst_vt['af_data']=[] # TODO add to "AeroDyn"
                for afi, af_filename in enumerate(self.fst_vt['AeroDyn15']['AFNames']):
                    af_filename = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']),af_filename)
                    polar = self._read(af_filename, 'AF')
                    self.fst_vt['af_data'].append(polar)
            # --- Backward compatibility
            self.AD  = self.fst_vt['AeroDyn15']
            self.ADversion='AD15'

        elif self.version=='OF2':
            # ---- Regular OpenFAST file
            # ElastoDyn
            if 'EDFile' in self.fst_vt['Fst'].keys():
                self.fst_vt['ElastoDyn'] = self._read(self.fst_vt['Fst']['EDFile'],'ED')
                if self.fst_vt['ElastoDyn'] is not None:
                    twr_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['EDFile']), self.fst_vt['ElastoDyn']['TwrFile'])
                    try:
                        bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['EDFile']), self.fst_vt['ElastoDyn']['BldFile(1)'])
                    except:
                        bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['EDFile']), self.fst_vt['ElastoDyn']['BldFile1'])
                    self.fst_vt['ElastoDynTower'] = self._read(twr_file,'EDtwr')
                    self.fst_vt['ElastoDynBlade'] = self._read(bld_file,'EDbld')

            # InflowWind
            if self.fst_vt['Fst']['CompInflow']>0:
                self.fst_vt['InflowWind'] = self._read(self.fst_vt['Fst']['InflowFile'],'IW')

            # AeroDyn
            if self.fst_vt['Fst']['CompAero']>0:
                key = 'AeroDyn14' if self.fst_vt['Fst']['CompAero']==1 else 'AeroDyn15'
                self.fst_vt[key] = self._read(self.fst_vt['Fst']['AeroFile'],'AD')

                if self.fst_vt[key] is not None:
                    # Blades
                    bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']), self.fst_vt[key]['ADBlFile(1)'])
                    self.fst_vt['AeroDynBlade'] = self._read(bld_file,'ADbld')
                    #self.fst_vt['AeroDynBlade'] = []
                    #for i in range(3):
                    #    bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']), self.fst_vt[key]['ADBlFile({})'.format(i+1)])
                    #    self.fst_vt['AeroDynBlade'].append(self._read(bld_file,'ADbld'))
                    # Polars
                    self.fst_vt['af_data']=[] # TODO add to "AeroDyn"
                    for afi, af_filename in enumerate(self.fst_vt['AeroDyn15']['AFNames']):
                        af_filename = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']),af_filename)
                        polar = self._read(af_filename, 'AF')
                        self.fst_vt['af_data'].append(polar)

            # ServoDyn
            if self.fst_vt['Fst']['CompServo']>0:
                self.fst_vt['ServoDyn'] = self._read(self.fst_vt['Fst']['ServoFile'],'SrvD')
                # TODO Discon

            # HydroDyn
            if self.fst_vt['Fst']['CompHydro']== 1:
                self.fst_vt['HydroDyn'] = self._read(self.fst_vt['Fst']['HydroFile'],'HD')

            # SubDyn
            if self.fst_vt['Fst']['CompSub'] == 1:
                self.fst_vt['SubDyn'] = self._read(self.fst_vt['Fst']['SubFile'],'HD')

            # Mooring
            if self.fst_vt['Fst']['CompMooring']==1:
                self.fst_vt['MAP'] = self._read(self.fst_vt['Fst']['MooringFile'],'MD')
            if self.fst_vt['Fst']['CompMooring']==2:
                self.fst_vt['MoorDyn'] = self._read(self.fst_vt['Fst']['MooringFile'],'MD')

            # BeamDyn
            if self.fst_vt['Fst']['CompElast'] == 2:
                self.fst_vt['BeamDyn'] = self._read(self.fst_vt['Fst']['BDBldFile(1)'],'BD')
                if self.fst_vt['BeamDyn'] is not None:
                    # Blades
                    bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['BDBldFile(1)']), self.fst_vt['BeamDyn']['BldFile'])
                    self.fst_vt['BeamDynBlade']= self._read(bld_file,'BDbld')
            # --- Backward compatibility
            if self.fst_vt['Fst']['CompAero']==1:
                self.AD  = self.fst_vt['AeroDyn14']
                self.ADversion='AD14'
            else:
                self.AD  = self.fst_vt['AeroDyn15']
                self.ADversion='AD15'

        elif self.version=='F7':
            raise NotImplementedError('')

        # --- Backward compatibility
        self.fst = self.fst_vt['Fst']
        self.ED  = self.fst_vt['ElastoDyn']
        if self.AD is not None:
            self.AD.Bld1 = self.fst_vt['AeroDynBlade']
            self.AD.AF  = self.fst_vt['af_data']
        self.IW  = self.fst_vt['InflowWind']
        self.BD  = self.fst_vt['BeamDyn']

    @ property
    def unusedNames(self):
        return ['unused','nan','na','none']

    def _read(self, relfilepath, shortkey):
        """ read any openfast input """
        relfilepath =relfilepath.replace('"','')
        basename = os.path.basename(relfilepath)

        # Only read what the user requested to be read
        if shortkey not in self.readlist:
            if self.verbose:
                print('>>> Skipping ',shortkey)
            return None

        # Skip "unused" and "NA"
        if basename.lower() in self.unusedNames:
            if self.verbose:
                print('>>> Unused ',shortkey)
            return None

        # Attempt reading
        fullpath =os.path.join(self.FAST_directory, relfilepath)
        try:
            data = FASTInputFile(fullpath)
            if self.verbose:
                print('>>> Read: ',fullpath)
            self.inputfiles[shortkey] = fullpath
            return data
        except FileNotFoundError:
            print('[WARN] File not found '+fullpath)
            return None

    def __repr__(self):
        s='<weio.FastInputDeck object>'+'\n'
        s+='filename   : '+self.filename+'\n'
        s+='version    : '+self.version+'\n'
        s+='AD version : '+self.ADversion+'\n'
        s+='fst_vt     : dict{'+','.join([k for k,v in self.fst_vt.items() if v is not None])+'}'
        s+='\n'
        return s

if __name__ == "__main__":
    fst=FASTInputDeck('NREL5MW.fst')
    print(fst)
