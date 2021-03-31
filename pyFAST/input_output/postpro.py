# --- For cmd.py
from __future__ import division, print_function
import os
import pandas as pd
import numpy as np
import re

# --- fast libraries
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.fast_output_file import FASTOutputFile
from pyFAST.input_output.fast_input_deck import FASTInputDeck

# --------------------------------------------------------------------------------}
# --- Tools for IO 
# --------------------------------------------------------------------------------{
def ED_BldStations(ED):
    """ Returns ElastoDyn Blade Station positions, useful to know where the outputs are.
    INPUTS:
       - ED: either:
           - a filename of a ElastoDyn input file
           - an instance of FileCl, as returned by reading the file, ED = weio.read(ED_filename)

    OUTUPTS:
        - bld_fract: fraction of the blade length were stations are defined
        - r_nodes: spanwise position from the rotor apex of the Blade stations
    """
    if not isinstance(ED,FASTInputFile):
        ED = FASTInputFile(ED)

    nBldNodes = ED['BldNodes']
    bld_fract    = np.arange(1./nBldNodes/2., 1, 1./nBldNodes)
    r_nodes      = bld_fract*(ED['TipRad']-ED['HubRad']) + ED['HubRad']
    return bld_fract, r_nodes

def ED_TwrStations(ED):
    """ Returns ElastoDyn Tower Station positions, useful to know where the outputs are.
    INPUTS:
       - ED: either:
           - a filename of a ElastoDyn input file
           - an instance of FileCl, as returned by reading the file, ED = weio.read(ED_filename)

    OUTPUTS:
        - r_fract: fraction of the towet length were stations are defined
        - h_nodes: height from the *ground* of the stations  (not from the Tower base)
    """
    if not isinstance(ED,FASTInputFile):
        ED = FASTInputFile(ED)

    nTwrNodes = ED['TwrNodes']
    twr_fract    = np.arange(1./nTwrNodes/2., 1, 1./nTwrNodes)
    h_nodes      = twr_fract*(ED['TowerHt']-ED['TowerBsHt']) + ED['TowerBsHt']
    return twr_fract, h_nodes



def ED_BldGag(ED):
    """ Returns the radial position of ElastoDyn blade gages 
    INPUTS:
       - ED: either:
           - a filename of a ElastoDyn input file
           - an instance of FileCl, as returned by reading the file, ED = weio.read(ED_filename)
    OUTPUTS:
       - r_gag: The radial positions of the gages, given from the rotor apex
    """
    if not isinstance(ED,FASTInputFile):
        ED = FASTInputFile(ED)
    _,r_nodes= ED_BldStations(ED)
    
    #     if ED.hasNodal:
    #         return r_nodes, None
    nOuts = ED['NBlGages']
    if nOuts<=0:
        return np.array([]), np.array([])
    if type(ED['BldGagNd']) is list:
        Inodes = np.asarray(ED['BldGagNd'])
    else:
        Inodes = np.array([ED['BldGagNd']])
    r_gag = r_nodes[ Inodes[:nOuts] -1]
    return r_gag, Inodes

def ED_TwrGag(ED):
    """ Returns the heights of ElastoDyn blade gages 
    INPUTS:
       - ED: either:
           - a filename of a ElastoDyn input file
           - an instance of FileCl, as returned by reading the file, ED = weio.read(ED_filename)
    OUTPUTS:
       - h_gag: The heights of the gages, given from the ground height (tower base + TowerBsHt)
    """
    if not isinstance(ED,FASTInputFile):
        ED = FASTInputFile(ED)
    _,h_nodes= ED_TwrStations(ED)
    nOuts = ED['NTwGages']
    if nOuts<=0:
        return np.array([])
    if type(ED['TwrGagNd']) is list:
        Inodes = np.asarray(ED['TwrGagNd'])
    else:
        Inodes = np.array([ED['TwrGagNd']])
    h_gag = h_nodes[ Inodes[:nOuts] -1]
    return h_gag


def AD14_BldGag(AD):
    """ Returns the radial position of AeroDyn 14 blade gages (based on "print" in column 6)
    INPUTS:
       - AD: either:
           - a filename of a AeroDyn input file
           - an instance of FileCl, as returned by reading the file, AD = weio.read(AD_filename)
    OUTPUTS:
       - r_gag: The radial positions of the gages, given from the blade root
    """
    if not isinstance(AD,FASTInputFile):
        AD = FASTInputFile(AD)

    Nodes=AD['BldAeroNodes']  
    if Nodes.shape[1]==6:
       doPrint= np.array([ n.lower().find('p')==0  for n in Nodes[:,5]])
    else:
       doPrint=np.array([ True  for n in Nodes[:,0]])

    r_gag = Nodes[doPrint,0].astype(float)
    IR    = np.arange(1,len(Nodes)+1)[doPrint]
    return r_gag, IR

def AD_BldGag(AD,AD_bld,chordOut=False):
    """ Returns the radial position of AeroDyn blade gages 
    INPUTS:
       - AD: either:
           - a filename of a AeroDyn input file
           - an instance of FileCl, as returned by reading the file, AD = weio.read(AD_filename)
       - AD_bld: either:
           - a filename of a AeroDyn Blade input file
           - an instance of FileCl, as returned by reading the file, AD_bld = weio.read(AD_bld_filename)
    OUTPUTS:
       - r_gag: The radial positions of the gages, given from the blade root
    """
    if not isinstance(AD,FASTInputFile):
        AD = FASTInputFile(AD)
    if not isinstance(AD_bld,FASTInputFile):
        AD_bld = FASTInputFile(AD_bld)
    #print(AD_bld.keys())

    nOuts=AD['NBlOuts']
    if nOuts<=0:
        if chordOut:
            return np.array([]), np.array([])
        else:
            return np.array([])
    INodes = np.array(AD['BlOutNd'][:nOuts])
    r_gag = AD_bld['BldAeroNodes'][INodes-1,0]
    if chordOut:
        chord_gag = AD_bld['BldAeroNodes'][INodes-1,5]
        return r_gag,chord_gag
    else:
        return r_gag

def BD_BldGag(BD):
    """ Returns the radial position of BeamDyn blade gages 
    INPUTS:
       - BD: either:
           - a filename of a BeamDyn input file
           - an instance of FileCl, as returned by reading the file, BD = weio.read(BD_filename)
    OUTPUTS:
       - r_gag: The radial positions of the gages, given from the rotor apex
    """
    if not isinstance(BD,FASTInputFile):
        BD = FASTInputFile(BD)

    M       = BD['MemberGeom']
    r_nodes = M[:,2] # NOTE: we select the z axis here, and we don't take curvilenear coord
    nOuts = BD['NNodeOuts']
    if nOuts<=0:
        nOuts=0
    if type(BD['OutNd']) is list:
        Inodes = np.asarray(BD['OutNd'])
    else:
        Inodes = np.array([BD['OutNd']])
    r_gag = r_nodes[ Inodes[:nOuts] -1]
    return r_gag, Inodes, r_nodes

# 
# 
# 1, 7, 14, 21, 30, 36, 43, 52, 58 BldGagNd List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]

# --------------------------------------------------------------------------------}
# --- Helper functions for radial data  
# --------------------------------------------------------------------------------{
def _HarmonizeSpanwiseData(Name, Columns, vr, R, IR=None) :
    """ helper function to use with spanwiseAD and spanwiseED """
    # --- Data present
    data     = [c for _,c in Columns if c is not None]
    ColNames = [n for n,_ in Columns if n is not None]
    Lengths  = [len(d) for d in data]
    if len(data)<=0:
        print('[WARN] No spanwise data for '+Name)
        return None, None, None

    # --- Harmonize data so that they all have the same length
    nrMax = np.max(Lengths)
    ids=np.arange(nrMax)
    if vr is None:
        bFakeVr=True
        vr_bar = ids/(nrMax-1)
    else:
        vr_bar=vr/R
        bFakeVr=False
        if (nrMax)<len(vr_bar):
            vr_bar=vr_bar[1:nrMax]
        elif (nrMax)>len(vr_bar):
            raise Exception('Inconsitent length between radial stations and max index present in output chanels')

    for i in np.arange(len(data)):
        d=data[i]
        if len(d)<nrMax:
            Values = np.zeros((nrMax,1))
            Values[:] = np.nan
            Values[1:len(d)] = d
            data[i] = Values

    # --- Combine data and remove 
    dataStack = np.column_stack([d for d in data])
    ValidRow = np.logical_not([np.isnan(dataStack).all(axis=1)])
    dataStack = dataStack[ValidRow[0],:]
    ids       = ids      [ValidRow[0]]
    vr_bar    = vr_bar   [ValidRow[0]]

    # --- Create a dataframe
    dfRad = pd.DataFrame(data= dataStack, columns = ColNames)

    if bFakeVr:
        dfRad.insert(0, 'i/n_[-]', vr_bar)
    else:
        dfRad.insert(0, 'r/R_[-]', vr_bar)
        if R is not None:
            r = vr_bar*R
    if IR is not None:
        dfRad['Node_[#]']=IR[:nrMax]
    dfRad['i_[#]']=ids+1
    if not bFakeVr:
        dfRad['r_[m]'] = r

    return dfRad,  nrMax, ValidRow

def insert_radial_columns(df, vr=None, R=None, IR=None):
    """
    Add some columns to the radial data
    """
    if df is None:
        return df
    if df.shape[1]==0:
        return None
    nrMax=len(df)
    ids=np.arange(nrMax)
    if vr is None or R is None:
        # Radial position unknown
        vr_bar = ids/(nrMax-1)
        df.insert(0, 'i/n_[-]', vr_bar)
    else:
        vr_bar=vr/R
        if (nrMax)<=len(vr_bar):
            vr_bar=vr_bar[:nrMax]
        elif (nrMax)>len(vr_bar):
            print(vr_bar)
            raise Exception('Inconsitent length between radial stations ({:d}) and max index present in output chanels ({:d})'.format(len(vr_bar),nrMax))
        df.insert(0, 'r/R_[-]', vr_bar)

    if IR is not None:
        df['Node_[#]']=IR[:nrMax]
    df['i_[#]']=ids+1
    if vr is not None:
        df['r_[m]'] = vr[:nrMax]
    return df

def find_matching_columns(Cols, PatternMap):
    ColsInfo=[]
    nrMax=0
    for colpattern,colmap in PatternMap.items():
        # Extracting columns matching pattern
        cols, sIdx = find_matching_pattern(Cols, colpattern)
        if len(cols)>0:
            # Sorting by ID
            cols  = np.asarray(cols)
            Idx   = np.array([int(s) for s in sIdx])
            Isort = np.argsort(Idx)
            Idx   = Idx[Isort]
            cols  = cols[Isort]
            col={'name':colmap,'Idx':Idx,'cols':cols}
            nrMax=max(nrMax,np.max(Idx))
            ColsInfo.append(col)
    return ColsInfo,nrMax

def extract_spanwise_data(ColsInfo, nrMax, df=None,ts=None):
    """ 
    Extract spanwise data based on some column info
    ColsInfo: see find_matching_columns
    """
    nCols = len(ColsInfo)
    if nCols==0:
        return None
    if ts is not None:
        Values = np.zeros((nrMax,nCols))
        Values[:] = np.nan
    elif df is not None:
        raise NotImplementedError()

    ColNames =[c['name'] for c in ColsInfo]

    for ic,c in enumerate(ColsInfo):
        Idx, cols, colname = c['Idx'], c['cols'], c['name']
        for idx,col in zip(Idx,cols):
            Values[idx-1,ic]=ts[col]
        nMissing = np.sum(np.isnan(Values[:,ic]))
        if len(cols)<nrMax:
            #print(Values)
            print('[WARN] Not all values found for {}, missing {}/{}'.format(colname,nMissing,nrMax))
        if len(cols)>nrMax:
            print('[WARN] More values found for {}, found {}/{}'.format(colname,len(cols),nrMax))
    df = pd.DataFrame(data=Values, columns=ColNames)
    df = df.reindex(sorted(df.columns), axis=1)
    return df

def spanwiseColBD(Cols):
    """ Return column info, available columns and indices that contain BD spanwise data"""
    BDSpanMap=dict()
    for sB in ['B1','B2','B3']:
        BDSpanMap['^'+sB+r'N(\d)TDxr_\[m\]']=sB+'TDxr_[m]'
        BDSpanMap['^'+sB+r'N(\d)TDyr_\[m\]']=sB+'TDyr_[m]'
        BDSpanMap['^'+sB+r'N(\d)TDzr_\[m\]']=sB+'TDzr_[m]'
    return find_matching_columns(Cols, BDSpanMap)

def spanwiseColED(Cols):
    """ Return column info, available columns and indices that contain ED spanwise data"""
    EDSpanMap=dict()
    # All Outs
    for sB in ['B1','B2','B3']:
        EDSpanMap['^[A]*'+sB+r'N(\d*)ALx_\[m/s^2\]' ] = sB+'ALx_[m/s^2]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)ALy_\[m/s^2\]' ] = sB+'ALy_[m/s^2]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)ALz_\[m/s^2\]' ] = sB+'ALz_[m/s^2]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)TDx_\[m\]'     ] = sB+'TDx_[m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)TDy_\[m\]'     ] = sB+'TDy_[m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)TDz_\[m\]'     ] = sB+'TDz_[m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)RDx_\[deg\]'   ] = sB+'RDx_[deg]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)RDy_\[deg\]'   ] = sB+'RDy_[deg]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)RDz_\[deg\]'   ] = sB+'RDz_[deg]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)MLx_\[kN-m\]'  ] = sB+'MLx_[kN-m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)MLy_\[kN-m\]'  ] = sB+'MLy_[kN-m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)MLz_\[kN-m\]'  ] = sB+'MLz_[kN-m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)FLx_\[kN\]'    ] = sB+'FLx_[kN]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)FLy_\[kN\]'    ] = sB+'FLy_[kN]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)FLz_\[kN\]'    ] = sB+'FLz_[kN]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)FLxNT_\[kN\]'  ] = sB+'FLxNT_[kN]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)FLyNT_\[kN\]'  ] = sB+'FLyNT_[kN]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)FlyNT_\[kN\]'  ] = sB+'FLyNT_[kN]'   # <<< Unfortunate
        EDSpanMap['^[A]*'+sB+r'N(\d*)MLxNT_\[kN-m\]'] = sB+'MLxNT_[kN-m]'
        EDSpanMap['^[A]*'+sB+r'N(\d*)MLyNT_\[kN-m\]'] = sB+'MLyNT_[kN-m]'
    # Old
    for sB in ['b1','b2','b3']:
        SB=sB.upper()
        EDSpanMap[r'^Spn(\d)ALx'+sB+r'_\[m/s^2\]']=SB+'ALx_[m/s^2]'
        EDSpanMap[r'^Spn(\d)ALy'+sB+r'_\[m/s^2\]']=SB+'ALy_[m/s^2]'
        EDSpanMap[r'^Spn(\d)ALz'+sB+r'_\[m/s^2\]']=SB+'ALz_[m/s^2]'
        EDSpanMap[r'^Spn(\d)TDx'+sB+r'_\[m\]'    ]=SB+'TDx_[m]'
        EDSpanMap[r'^Spn(\d)TDy'+sB+r'_\[m\]'    ]=SB+'TDy_[m]'
        EDSpanMap[r'^Spn(\d)TDz'+sB+r'_\[m\]'    ]=SB+'TDz_[m]'
        EDSpanMap[r'^Spn(\d)RDx'+sB+r'_\[deg\]'  ]=SB+'RDx_[deg]'
        EDSpanMap[r'^Spn(\d)RDy'+sB+r'_\[deg\]'  ]=SB+'RDy_[deg]'
        EDSpanMap[r'^Spn(\d)RDz'+sB+r'_\[deg\]'  ]=SB+'RDz_[deg]'
        EDSpanMap[r'^Spn(\d)FLx'+sB+r'_\[kN\]'   ]=SB+'FLx_[kN]'
        EDSpanMap[r'^Spn(\d)FLy'+sB+r'_\[kN\]'   ]=SB+'FLy_[kN]'
        EDSpanMap[r'^Spn(\d)FLz'+sB+r'_\[kN\]'   ]=SB+'FLz_[kN]'
        EDSpanMap[r'^Spn(\d)MLy'+sB+r'_\[kN-m\]' ]=SB+'MLx_[kN-m]'
        EDSpanMap[r'^Spn(\d)MLx'+sB+r'_\[kN-m\]' ]=SB+'MLy_[kN-m]'  
        EDSpanMap[r'^Spn(\d)MLz'+sB+r'_\[kN-m\]' ]=SB+'MLz_[kN-m]'
    return find_matching_columns(Cols, EDSpanMap)

def spanwiseColAD(Cols):
    """ Return column info, available columns and indices that contain AD spanwise data"""
    ADSpanMap=dict()
    for sB in ['B1','B2','B3']:
        ADSpanMap['^[A]*'+sB+r'N(\d*)Alpha_\[deg\]']=sB+'Alpha_[deg]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)AOA_\[deg\]'  ]=sB+'Alpha_[deg]' # DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)AxInd_\[-\]'  ]=sB+'AxInd_[-]'  
        ADSpanMap['^[A]*'+sB+r'N(\d*)TnInd_\[-\]'  ]=sB+'TnInd_[-]'  
        ADSpanMap['^[A]*'+sB+r'N(\d*)AIn_\[deg\]'  ]=sB+'AxInd_[-]'   # DBGOuts NOTE BUG Unit
        ADSpanMap['^[A]*'+sB+r'N(\d*)ApI_\[deg\]'  ]=sB+'TnInd_[-]'   # DBGOuts NOTE BUG Unit
        ADSpanMap['^[A]*'+sB+r'N(\d*)AIn_\[-\]'    ]=sB+'AxInd_[-]'   # DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)ApI_\[-\]'    ]=sB+'TnInd_[-]'   # DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)Uin_\[m/s\]'  ]=sB+'Uin_[m/s]'     # DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)Uit_\[m/s\]'  ]=sB+'Uit_[m/s]'     # DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)Uir_\[m/s\]'  ]=sB+'Uir_[m/s]'     # DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)Cl_\[-\]'     ]=sB+'Cl_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Cd_\[-\]'     ]=sB+'Cd_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Cm_\[-\]'     ]=sB+'Cm_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Cx_\[-\]'     ]=sB+'Cx_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Cy_\[-\]'     ]=sB+'Cy_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Cn_\[-\]'     ]=sB+'Cn_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Ct_\[-\]'     ]=sB+'Ct_[-]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Re_\[-\]'     ]=sB+'Re_[-]' 
        ADSpanMap['^[A]*'+sB+r'N(\d*)Vrel_\[m/s\]' ]=sB+'Vrel_[m/s]' 
        ADSpanMap['^[A]*'+sB+r'N(\d*)Theta_\[deg\]']=sB+'Theta_[deg]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Phi_\[deg\]'  ]=sB+'Phi_[deg]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Twst_\[deg\]' ]=sB+'Twst_[deg]' #DBGOuts
        ADSpanMap['^[A]*'+sB+r'N(\d*)Curve_\[deg\]']=sB+'Curve_[deg]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Vindx_\[m/s\]']=sB+'Vindx_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Vindy_\[m/s\]']=sB+'Vindy_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Fx_\[N/m\]'   ]=sB+'Fx_[N/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Fy_\[N/m\]'   ]=sB+'Fy_[N/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Fl_\[N/m\]'   ]=sB+'Fl_[N/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Fd_\[N/m\]'   ]=sB+'Fd_[N/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Fn_\[N/m\]'   ]=sB+'Fn_[N/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Ft_\[N/m\]'   ]=sB+'Ft_[N/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)VUndx_\[m/s\]']=sB+'VUndx_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)VUndy_\[m/s\]']=sB+'VUndy_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)VUndz_\[m/s\]']=sB+'VUndz_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)VDisx_\[m/s\]']=sB+'VDisx_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)VDisy_\[m/s\]']=sB+'VDisy_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)VDisz_\[m/s\]']=sB+'VDisz_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)STVx_\[m/s\]' ]=sB+'STVx_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)STVy_\[m/s\]' ]=sB+'STVy_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)STVz_\[m/s\]' ]=sB+'STVz_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Vx_\[m/s\]'   ]=sB+'Vx_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Vy_\[m/s\]'   ]=sB+'Vy_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)Vz_\[m/s\]'   ]=sB+'Vz_[m/s]'
        ADSpanMap['^[A]*'+sB+r'N(\d*)DynP_\[Pa\]'  ]=sB+'DynP_[Pa]' 
        ADSpanMap['^[A]*'+sB+r'N(\d*)M_\[-\]'      ]=sB+'M_[-]' 
        ADSpanMap['^[A]*'+sB+r'N(\d*)Mm_\[N-m/m\]' ]=sB+'Mm_[N-m/m]'   
        ADSpanMap['^[A]*'+sB+r'N(\d*)Gam_\['       ]=sB+'Gam_[m^2/s]' #DBGOuts
    # --- AD 14
    ADSpanMap[r'^Alpha(\d*)_\[deg\]'  ]='Alpha_[deg]'  
    ADSpanMap[r'^DynPres(\d*)_\[Pa\]' ]='DynPres_[Pa]' 
    ADSpanMap[r'^CLift(\d*)_\[-\]'    ]='CLift_[-]'    
    ADSpanMap[r'^CDrag(\d*)_\[-\]'    ]='CDrag_[-]'    
    ADSpanMap[r'^CNorm(\d*)_\[-\]'    ]='CNorm_[-]'    
    ADSpanMap[r'^CTang(\d*)_\[-\]'    ]='CTang_[-]'    
    ADSpanMap[r'^CMomt(\d*)_\[-\]'    ]='CMomt_[-]'    
    ADSpanMap[r'^Pitch(\d*)_\[deg\]'  ]='Pitch_[deg]'  
    ADSpanMap[r'^AxInd(\d*)_\[-\]'    ]='AxInd_[-]'    
    ADSpanMap[r'^TanInd(\d*)_\[-\]'   ]='TanInd_[-]'   
    ADSpanMap[r'^ForcN(\d*)_\[N\]'    ]='ForcN_[N]'    
    ADSpanMap[r'^ForcT(\d*)_\[N\]'    ]='ForcT_[N]'    
    ADSpanMap[r'^Pmomt(\d*)_\[N-m\]'  ]='Pmomt_[N-N]'  
    ADSpanMap[r'^ReNum(\d*)_\[x10^6\]']='ReNum_[x10^6]'
    ADSpanMap[r'^Gamma(\d*)_\[m^2/s\]']='Gamma_[m^2/s]'

    return find_matching_columns(Cols, ADSpanMap)

def insert_extra_columns_AD(dfRad, tsAvg, vr=None, rho=None, R=None, nB=None, chord=None):
    # --- Compute additional values (AD15 only)
    if dfRad is None:
        return None
    if dfRad.shape[1]==0:
        return dfRad
    if chord is not None:
        if vr is not None:
            chord =chord[0:len(dfRad)]
    for sB in ['B1','B2','B3']:
        try:
            vr_bar=vr/R
            Fx = dfRad[sB+'Fx_[N/m]']
            U0 = tsAvg['Wind1VelX_[m/s]']
            Ct=nB*Fx/(0.5 * rho * 2 * U0**2 * np.pi * vr)
            Ct[vr<0.01*R] = 0
            dfRad[sB+'Ctloc_[-]'] = Ct
            CT=2*np.trapz(vr_bar*Ct,vr_bar)
            dfRad[sB+'CtAvg_[-]']= CT*np.ones(vr.shape)
        except:
            pass
        try:
            dfRad[sB+'Gamma_[m^2/s]'] = 1/2 * chord*  dfRad[sB+'Vrel_[m/s]'] * dfRad[sB+'Cl_[-]'] 
        except:
            pass
        try: 
            if not sB+'Vindx_[m/s]' in dfRad.columns:
                dfRad[sB+'Vindx_[m/s]']= -dfRad[sB+'AxInd_[-]'].values * dfRad[sB+'Vx_[m/s]'].values 
                dfRad[sB+'Vindy_[m/s]']=  dfRad[sB+'TnInd_[-]'].values * dfRad[sB+'Vy_[m/s]'].values 
        except:
            pass
    return dfRad



def spanwisePostPro(FST_In=None,avgMethod='constantwindow',avgParam=5,out_ext='.outb',df=None):
    """
    Postprocess FAST radial data

    INPUTS:
        - FST_IN: Fast .fst input file
        - avgMethod='periods', avgParam=2:  average over 2 last periods, Needs Azimuth sensors!!!
        - avgMethod='constantwindow', avgParam=5:  average over 5s of simulation
        - postprofile: outputfile to write radial data
    """
    # --- Opens Fast output  and performs averaging
    if df is None:
        df = FASTOutputFile(FST_In.replace('.fst',out_ext).replace('.dvr',out_ext)).toDataFrame()
        returnDF=True
    else:
        returnDF=False
    # NOTE: spanwise script doest not support duplicate columns
    df = df.loc[:,~df.columns.duplicated()]
    dfAvg = averageDF(df,avgMethod=avgMethod ,avgParam=avgParam) # NOTE: average 5 last seconds

    # --- Extract info (e.g. radial positions) from Fast input file
    # We don't have a .fst input file, so we'll rely on some default values for "r"
    rho         = 1.225
    chord       = None
    # --- Extract radial positions of output channels
    r_AD, r_ED, r_BD, IR_AD, IR_ED, IR_BD, R, r_hub, fst = FASTRadialOutputs(FST_In, OutputCols=df.columns.values)
    if R is None: 
        R=1
    try:
        chord  = fst.AD.Bld1['BldAeroNodes'][:,5] # Full span
    except:
        pass
    try:
        rho = fst.AD['Rho']
    except:
        try:
            rho = fst.AD['AirDens']
        except:
            pass
    #print('r_AD:', r_AD)
    #print('r_ED:', r_ED)
    #print('r_BD:', r_BD)
    #print('I_AD:', IR_AD)
    #print('I_ED:', IR_ED)
    #print('I_BD:', IR_BD)
    # --- Extract radial data and export to csv if needed
    dfRad_AD    = None
    dfRad_ED    = None
    dfRad_BD    = None
    Cols=dfAvg.columns.values
    # --- AD
    ColsInfoAD, nrMaxAD = spanwiseColAD(Cols)
    dfRad_AD            = extract_spanwise_data(ColsInfoAD, nrMaxAD, df=None, ts=dfAvg.iloc[0])
    dfRad_AD            = insert_extra_columns_AD(dfRad_AD, dfAvg.iloc[0], vr=r_AD, rho=rho, R=R, nB=3, chord=chord)
    dfRad_AD            = insert_radial_columns(dfRad_AD, r_AD, R=R, IR=IR_AD)
    # --- ED
    ColsInfoED, nrMaxED = spanwiseColED(Cols)
    dfRad_ED            = extract_spanwise_data(ColsInfoED, nrMaxED, df=None, ts=dfAvg.iloc[0])
    dfRad_ED            = insert_radial_columns(dfRad_ED, r_ED, R=R, IR=IR_ED)
    # --- BD
    ColsInfoBD, nrMaxBD = spanwiseColBD(Cols)
    dfRad_BD            = extract_spanwise_data(ColsInfoBD, nrMaxBD, df=None, ts=dfAvg.iloc[0])
    dfRad_BD            = insert_radial_columns(dfRad_BD, r_BD, R=R, IR=IR_BD)
    if returnDF:
        return dfRad_ED , dfRad_AD, dfRad_BD, df
    else:
        return dfRad_ED , dfRad_AD, dfRad_BD



def spanwisePostProRows(df, FST_In=None):
    """ 
    Returns a 3D matrix: n x nSpan x nColumn where df is of size n x nColumn

    NOTE: this is really not optimal. Spanwise columns should be extracted only once..
    """
    # --- Extract info (e.g. radial positions) from Fast input file
    # We don't have a .fst input file, so we'll rely on some default values for "r"
    rho         = 1.225
    chord       = None
    # --- Extract radial positions of output channels
    r_AD, r_ED, r_BD, IR_AD, IR_ED, IR_BD, R, r_hub, fst = FASTRadialOutputs(FST_In, OutputCols=df.columns.values)
    #print('r_AD:', r_AD)
    #print('r_ED:', r_ED)
    #print('r_BD:', r_BD)
    if R is None: 
        R=1
    try:
        chord  = fst.AD.Bld1['BldAeroNodes'][:,5] # Full span
    except:
        pass
    try:
        rho = fst.AD['Rho']
    except:
        try:
            rho = fst.AD['AirDens']
        except:
            pass
    # --- Extract radial data for each azimuthal average
    M_AD=None
    M_ED=None
    M_BD=None
    Col_AD=None
    Col_ED=None
    Col_BD=None
    v = df.index.values

    # --- Getting Column info
    Cols=df.columns.values
    if r_AD is not None:
        ColsInfoAD, nrMaxAD = spanwiseColAD(Cols)
    if r_ED is not None:
        ColsInfoED, nrMaxED = spanwiseColED(Cols)
    if r_BD is not None:
        ColsInfoBD, nrMaxBD = spanwiseColBD(Cols)
    for i,val in enumerate(v):
        if r_AD is not None:
            dfRad_AD = extract_spanwise_data(ColsInfoAD, nrMaxAD, df=None, ts=df.iloc[i])
            dfRad_AD = insert_extra_columns_AD(dfRad_AD, df.iloc[i], vr=r_AD, rho=rho, R=R, nB=3, chord=chord)
            dfRad_AD = insert_radial_columns(dfRad_AD, r_AD, R=R, IR=IR_AD)
            if i==0:
                M_AD = np.zeros((len(v), len(dfRad_AD), len(dfRad_AD.columns)))
                Col_AD=dfRad_AD.columns.values
            M_AD[i, :, : ] = dfRad_AD.values
        if r_ED is not None and len(r_ED)>0:
            dfRad_ED = extract_spanwise_data(ColsInfoED, nrMaxED, df=None, ts=df.iloc[i])
            dfRad_ED = insert_radial_columns(dfRad_ED, r_ED, R=R, IR=IR_ED)
            if i==0:
                M_ED = np.zeros((len(v), len(dfRad_ED), len(dfRad_ED.columns)))
                Col_ED=dfRad_ED.columns.values
            M_ED[i, :, : ] = dfRad_ED.values
        if r_BD is not None and len(r_BD)>0:
            dfRad_BD = extract_spanwise_data(ColsInfoBD, nrMaxBD, df=None, ts=df.iloc[i])
            dfRad_BD = insert_radial_columns(dfRad_BD, r_BD, R=R, IR=IR_BD)
            if i==0:
                M_BD = np.zeros((len(v), len(dfRad_BD), len(dfRad_BD.columns)))
                Col_BD=dfRad_BD.columns.values
            M_BD[i, :, : ] = dfRad_BD.values
    return M_AD, Col_AD, M_ED, Col_ED, M_BD, Col_BD


def FASTRadialOutputs(FST_In, OutputCols=None):
    """ Returns radial positions where FAST has outputs
    INPUTS:
       FST_In: fast input file (.fst)
    OUTPUTS:
       r_AD: radial positions of FAST Outputs from the rotor center
    """
    R           = None
    r_hub =0
    r_AD        = None 
    r_ED        = None
    r_BD        = None
    IR_ED       = None
    IR_AD       = None
    IR_BD       = None
    fst=None
    if FST_In is not None:
        fst = FASTInputDeck(FST_In, readlist=['AD','ADbld','ED','BD'])
        # NOTE: all this below should be in FASTInputDeck
        if fst.version == 'F7':
            # --- FAST7
            if  not hasattr(fst,'AD'):
                raise Exception('The AeroDyn file couldn''t be found or read, from main file: '+FST_In)
            r_AD,IR_AD = AD14_BldGag(fst.AD)
            R   = fst.fst['TipRad']
            try:
                rho = fst.AD['Rho']
            except:
                rho = fst.AD['AirDens']
        else:
            # --- OpenFAST 2
            R = None

            # --- ElastoDyn
            if 'NumTurbines' in fst.fst.keys():
                # AeroDyn driver...
                r_hub       = fst.fst['BldHubRad_bl(1_1)']

            elif  not hasattr(fst,'ED'):
                print('[WARN] The Elastodyn file couldn''t be found or read, from main file: '+FST_In)
                #raise Exception('The Elastodyn file couldn''t be found or read, from main file: '+FST_In)
            else:
                R           = fst.ED['TipRad']
                r_hub       = fst.ED['HubRad']
                if fst.ED.hasNodal:
                    _, r_ED = ED_BldStations(fst.ED)
                    IR_ED =None
                else:
                    r_ED, IR_ED = ED_BldGag(fst.ED)

            # --- BeamDyn
            if  fst.BD is not None:
                r_BD, IR_BD, r_BD_All = BD_BldGag(fst.BD)
                r_BD= r_BD+r_hub
                if R is None:
                    R = r_BD_All[-1] # just in case ED file missing

            # --- AeroDyn
            if  fst.AD is None:
                print('[WARN] The AeroDyn file couldn''t be found or read, from main file: '+FST_In)
                #raise Exception('The AeroDyn file couldn''t be found or read, from main file: '+FST_In)
            else:
                if fst.ADversion == 'AD15':
                    if  fst.AD.Bld1 is None:
                        raise Exception('The AeroDyn blade file couldn''t be found or read, from main file: '+FST_In)
                    
                    if 'B1N001Cl_[-]' in OutputCols or np.any(np.char.find(list(OutputCols),'AB1N')==0):
                        # This was compiled with all outs
                        r_AD   = fst.AD.Bld1['BldAeroNodes'][:,0] # Full span
                        r_AD   += r_hub
                        IR_AD  = None
                    else:
                        r_AD,_ = AD_BldGag(fst.AD,fst.AD.Bld1, chordOut = True) # Only at Gages locations
                        r_AD   += r_hub

                    if R is None:
                        # ElastoDyn was not read, we use R from AD
                        R = fst.AD.Bld1['BldAeroNodes'][-1,0]

                elif fst.ADversion == 'AD14':
                    r_AD,IR_AD = AD14_BldGag(fst.AD)

                else:
                    raise Exception('AeroDyn version unknown')
    return r_AD, r_ED, r_BD, IR_AD, IR_ED, IR_BD, R, r_hub, fst



def addToOutlist(OutList, Signals):
    if not isinstance(Signals,list):
        raise Exception('Signals must be a list')
    for s in Signals:
        ss=s.split()[0].strip().strip('"').strip('\'')
        AlreadyIn = any([o.find(ss)==1 for o in OutList ])
        if not AlreadyIn:
            OutList.append(s)
    return OutList



# --------------------------------------------------------------------------------}
# --- Generic df 
# --------------------------------------------------------------------------------{
def remap_df(df, ColMap, bColKeepNewOnly=False, inPlace=False):
    """ Add/rename columns of a dataframe, potentially perform operations between columns

    Example:

        ColumnMap={
          'WS_[m/s]'         : '{Wind1VelX_[m/s]}'             , # create a new column from existing one
          'RtTSR_[-]'        : '{RtTSR_[-]} * 2  +  {RtAeroCt_[-]}'    , # change value of column
          'RotSpeed_[rad/s]' : '{RotSpeed_[rpm]} * 2*np.pi/60 ', # new column [rpm] -> [rad/s]
        }
        # Read
        df = weio.read('FASTOutBin.outb').toDataFrame()
        # Change columns based on formulae, potentially adding new columns
        df = fastlib.remap_df(df, ColumnMap, inplace=True)

    """
    if not inPlace:
        df=df.copy()
    ColMapMiss=[]
    ColNew=[]
    RenameMap=dict()
    for k0,v in ColMap.items():
        k=k0.strip()
        v=v.strip()
        if v.find('{')>=0:
            search_results = re.finditer(r'\{.*?\}', v)
            expr=v
            # For more advanced operations, we use an eval
            bFail=False
            for item in search_results:
                col=item.group(0)[1:-1]
                if col not in df.columns:
                    ColMapMiss.append(col)
                    bFail=True
                expr=expr.replace(item.group(0),'df[\''+col+'\']')
            #print(k0, '=', expr)
            if not bFail:
                df[k]=eval(expr)
                ColNew.append(k)
            else:
                print('[WARN] Column not present in dataframe, cannot evaluate: ',expr)
        else:
            #print(k0,'=',v)
            if v not in df.columns:
                ColMapMiss.append(v)
                print('[WARN] Column not present in dataframe: ',v)
            else:
                RenameMap[k]=v

    # Applying renaming only now so that expressions may be applied in any order
    for k,v in RenameMap.items():
        k=k.strip()
        iCol = list(df.columns).index(v)
        df.columns.values[iCol]=k
        ColNew.append(k)
    df.columns = df.columns.values # Hack to ensure columns are updated

    if len(ColMapMiss)>0:
        print('[FAIL] The following columns were not found in the dataframe:',ColMapMiss)
        #print('Available columns are:',df.columns.values)

    if bColKeepNewOnly:
        ColNew = [c for c,_ in ColMap.items() if c in ColNew]# Making sure we respec order from user
        ColKeepSafe = [c for c in ColNew if c in df.columns.values]
        ColKeepMiss = [c for c in ColNew if c not in df.columns.values]
        if len(ColKeepMiss)>0:
            print('[WARN] Signals missing and omitted for ColKeep:\n       '+'\n       '.join(ColKeepMiss))
        df=df[ColKeepSafe]
    return df


# --------------------------------------------------------------------------------}
# --- Tools for PostProcessing one or several simulations
# --------------------------------------------------------------------------------{
def _zero_crossings(y,x=None,direction=None):
    """
      Find zero-crossing points in a discrete vector, using linear interpolation.
      direction: 'up' or 'down', to select only up-crossings or down-crossings
      Returns: 
          x values xzc such that y(yzc)==0
          indexes izc, such that the zero is between y[izc] (excluded) and y[izc+1] (included)
      if direction is not provided, also returns:
              sign, equal to 1 for up crossing
    """
    y=np.asarray(y)
    if x is None:
        x=np.arange(len(y))

    if np.any((x[1:] - x[0:-1]) <= 0.0):
        raise Exception('x values need to be in ascending order')

    # Indices before zero-crossing
    iBef = np.where(y[1:]*y[0:-1] < 0.0)[0]
    
    # Find the zero crossing by linear interpolation
    xzc = x[iBef] - y[iBef] * (x[iBef+1] - x[iBef]) / (y[iBef+1] - y[iBef])
    
    # Selecting points that are exactly 0 and where neighbor change sign
    iZero = np.where(y == 0.0)[0]
    iZero = iZero[np.where((iZero > 0) & (iZero < x.size-1))]
    iZero = iZero[np.where(y[iZero-1]*y[iZero+1] < 0.0)]

    # Concatenate 
    xzc  = np.concatenate((xzc, x[iZero]))
    iBef = np.concatenate((iBef, iZero))

    # Sort
    iSort = np.argsort(xzc)
    xzc, iBef = xzc[iSort], iBef[iSort]

    # Return up-crossing, down crossing or both
    sign = np.sign(y[iBef+1]-y[iBef])
    if direction == 'up':
        I= np.where(sign==1)[0]
        return xzc[I],iBef[I]
    elif direction == 'down':
        I= np.where(sign==-1)[0]
        return xzc[I],iBef[I]
    elif direction is not None:
        raise Exception('Direction should be either `up` or `down`')
    return xzc, iBef, sign

def find_matching_pattern(List, pattern):
    """ Return elements of a list of strings that match a pattern
        and return the first matching group
    """
    reg_pattern=re.compile(pattern)
    MatchedElements=[]
    MatchedStrings=[]
    for l in List:
        match=reg_pattern.search(l)
        if match:
            MatchedElements.append(l)
            if len(match.groups(1))>0:
                MatchedStrings.append(match.groups(1)[0])
            else:
                MatchedStrings.append('')
    return MatchedElements, MatchedStrings

        

def extractSpanTSReg(ts, col_pattern, colname, IR=None):
    """ Helper function to extract spanwise results, like B1N1Cl B1N2Cl etc. 

    Example
        col_pattern: 'B1N(\d*)Cl_\[-\]'
        colname    : 'B1Cl_[-]'
    """
    # Extracting columns matching pattern
    cols, sIdx = find_matching_pattern(ts.keys(), col_pattern)
    if len(cols) ==0:
        return (None,None)

    # Sorting by ID
    cols = np.asarray(cols)
    Idx  = np.array([int(s) for s in sIdx])
    Isort = np.argsort(Idx)
    Idx  = Idx[Isort]
    cols = cols[Isort]

    nrMax =  np.max(Idx)
    Values = np.zeros((nrMax,1))
    Values[:] = np.nan
#     if IR is None:
#         cols   = [col_pattern.format(ir+1) for ir in range(nr)]
#     else:
#         cols   = [col_pattern.format(ir) for ir in IR]
    for idx,col in zip(Idx,cols):
        Values[idx-1]=ts[col]
    nMissing = np.sum(np.isnan(Values))
    if nMissing==nrMax:
        return (None,None)
    if len(cols)<nrMax:
        #print(Values)
        print('[WARN] Not all values found for {}, missing {}/{}'.format(colname,nMissing,nrMax))
    if len(cols)>nrMax:
        print('[WARN] More values found for {}, found {}/{}'.format(colname,len(cols),nrMax))
    return (colname,Values)

def extractSpanTS(ts, nr, col_pattern, colname, IR=None):
    """ Helper function to extract spanwise results, like B1N1Cl B1N2Cl etc. 

    Example
        col_pattern: 'B1N{:d}Cl_[-]'
        colname    : 'B1Cl_[-]'
    """
    Values=np.zeros((nr,1))
    if IR is None:
        cols   = [col_pattern.format(ir+1) for ir in range(nr)]
    else:
        cols   = [col_pattern.format(ir) for ir in IR]
    colsExist  = [c for c in cols if c in ts.keys() ]
    if len(colsExist)==0:
        return (None,None)

    Values = [ts[c] if c in ts.keys() else np.nan for c in cols  ]
    nMissing = np.sum(np.isnan(Values))
    #Values = ts[cols].T
    #nCoun=len(Values)
    if nMissing==nr:
        return (None,None)
    if len(colsExist)<nr:
        print(Values)
        print('[WARN] Not all values found for {}, missing {}/{}'.format(colname,nMissing,nr))
    if len(colsExist)>nr:
        print('[WARN] More values found for {}, found {}/{}'.format(colname,len(cols),nr))
    return (colname,Values)

def radialInterpTS(df, r, varName, r_ref, blade=1, bldFmt='AB{:d}', ndFmt='N{:03d}', method='interp'):
    """ 
    Interpolate a time series at a given radial position for a given variable (varName)
    INPUTS:
     - df     : a dataframe (typically with OpenFAST time series)
     - r      : radial positions of node where data is to be interpolated
     - varName: variable name (and unit) to be interpolated. 
                The dataframe column will be assumed to be "BldFmt"+"ndFmt"+varName
     - r_ref  : radial position of nodal data present in the dataframe
     - bldFmt : format for blade number, e.g. 'B{:d}' or 'AB{:d}'
     - ndFmt  : format for node number, e.g.  'N{:d}' or 'N{:03d}'
    OUTPUT:
      - interpolated time series
    """
    # --- Sanity checks
    r_ref = np.asarray(r_ref)
    if not np.all(r_ref[:-1] <= r_ref[1:]):
        raise Exception('This function only works for ascending radial values')

    # No extrapolation
    if r<np.min(r_ref) or r>np.max(r_ref):
        raise Exception('Extrapolation not supported')

    # Exactly on first or last nodes
    if r==r_ref[0]:
        col=bldFmt.format(blade) + ndFmt.format(1) + varName
        if col in df.columns.values:
            return df[col]
        else:
            raise Exception('Column {} not found in dataframe'.format(col))
    elif r==r_ref[-1]:
        col=bldFmt.format(blade) + ndFmt.format(len(r_ref)+1) + varName
        if col in df.columns.values:
            return df[col]
        else:
            raise Exception('Column {} not found in dataframe'.format(col))

    if method=='interp':
        # Interpolation
        iBef = np.where(r_ref<r)[0][-1]
        iAft = iBef+1
        #print(r_ref[iBef], r,  r_ref[iAft], '          ',iBef+1, iAft+1)
        fact= np.interp(r, r_ref[iBef:iAft+1], [0,1])
        col=bldFmt.format(blade) + ndFmt.format(iBef+1) + varName
        if col in df.columns.values:
            bef=df[bldFmt.format(blade) + ndFmt.format(iBef+1) + varName]
            aft=df[bldFmt.format(blade) + ndFmt.format(iAft+1) + varName]
        else:
            raise Exception('Column {} not found in dataframe'.format(col))
        return bef*(1-fact) + aft*fact
    else: 
        raise NotImplementedError()



def bin_mean_DF(df, xbins, colBin ):
    """ 
    Perform bin averaging of a dataframe
    """
    if colBin not in df.columns.values:
        raise Exception('The column `{}` does not appear to be in the dataframe'.format(colBin))
    xmid      = (xbins[:-1]+xbins[1:])/2
    df['Bin'] = pd.cut(df[colBin], bins=xbins, labels=xmid ) # Adding a column that has bin attribute
    df2       = df.groupby('Bin').mean()                     # Average by bin
    # also counting
    df['Counts'] = 1
    dfCount=df[['Counts','Bin']].groupby('Bin').sum()
    df2['Counts'] = dfCount['Counts']
    # Just in case some bins are missing (will be nan)
    df2       = df2.reindex(xmid)
    return df2

def azimuthal_average_DF(df, psiBin=None, colPsi='Azimuth_[deg]', tStart=None, colTime='Time_[s]'):
    """ 
    Average a dataframe based on azimuthal value
    Returns a dataframe with same amount of columns as input, and azimuthal values as index
    """
    if psiBin is None: 
        psiBin = np.arange(0,360+1,10)

    if tStart is not None:
        if colTime not in df.columns.values:
            raise Exception('The column `{}` does not appear to be in the dataframe'.format(colTime))
        df=df[ df[colTime]>tStart].copy()

    dfPsi= bin_mean_DF(df, psiBin, colPsi)
    if np.any(dfPsi['Counts']<1):
        print('[WARN] some bins have no data! Increase the bin size.')

    return dfPsi


def averageDF(df,avgMethod='periods',avgParam=None,ColMap=None,ColKeep=None,ColSort=None,stats=['mean']):
    """
    See average PostPro for documentation, same interface, just does it for one dataframe
    """
    def renameCol(x):
        for k,v in ColMap.items():
            if x==v:
                return k
        return x
    # Before doing the colomn map we store the time
    time = df['Time_[s]'].values
    timenoNA = time[~np.isnan(time)]
    # Column mapping
    if ColMap is not None:
        ColMapMiss = [v for _,v in ColMap.items() if v not in df.columns.values]
        if len(ColMapMiss)>0:
            print('[WARN] Signals missing and omitted for ColMap:\n       '+'\n       '.join(ColMapMiss))
        df.rename(columns=renameCol,inplace=True)
    ## Defining a window for stats (start time and end time)
    if avgMethod.lower()=='constantwindow':
        tEnd = timenoNA[-1]
        if avgParam is None:
            tStart=timenoNA[0]
        else:
            tStart =tEnd-avgParam
    elif avgMethod.lower()=='periods':
        # --- Using azimuth to find periods
        if 'Azimuth_[deg]' not in df.columns:
            raise Exception('The sensor `Azimuth_[deg]` does not appear to be in the output file. You cannot use the averaging method by `periods`, use `constantwindow` instead.')
        # NOTE: potentially we could average over each period and then average
        psi=df['Azimuth_[deg]'].values
        _,iBef = _zero_crossings(psi-psi[-10],direction='up')
        if len(iBef)==0:
            _,iBef = _zero_crossings(psi-180,direction='up')
        if len(iBef)==0:
            print('[WARN] Not able to find a zero crossing!')
            tEnd = time[-1]
            iBef=[0]
        else:
            tEnd = time[iBef[-1]]

        if avgParam is None:
            tStart=time[iBef[0]]
        else:
            avgParam=int(avgParam) 
            if len(iBef)-1<avgParam:
                print('[WARN] Not enough periods found ({}) compared to number requested to average ({})!'.format(len(iBef)-1,avgParam))
                avgParam=len(iBef)-1
            if avgParam==0:
                tStart = time[0]
                tEnd   = time[-1]
            else:
                tStart=time[iBef[-1-avgParam]]
    elif avgMethod.lower()=='periods_omega':
        # --- Using average omega to find periods
        if 'RotSpeed_[rpm]' not in df.columns:
            raise Exception('The sensor `RotSpeed_[rpm]` does not appear to be in the output file. You cannot use the averaging method by `periods_omega`, use `periods` or `constantwindow` instead.')
        Omega=df['RotSpeed_[rpm]'].mean()/60*2*np.pi
        Period = 2*np.pi/Omega 
        if avgParam is None:
            nRotations=np.floor(tEnd/Period)
        else:
            nRotations=avgParam
        tStart =tEnd-Period*nRotations
    else:
        raise Exception('Unknown averaging method {}'.format(avgMethod))
    # Narrowind number of columns here (azimuth needed above)
    if ColKeep is not None:
        ColKeepSafe = [c for c in ColKeep if c in df.columns.values]
        ColKeepMiss = [c for c in ColKeep if c not in df.columns.values]
        if len(ColKeepMiss)>0:
            print('[WARN] Signals missing and omitted for ColKeep:\n       '+'\n       '.join(ColKeepMiss))
        df=df[ColKeepSafe]
    if tStart<time[0]:
        print('[WARN] Simulation time ({}) too short compared to required averaging window ({})!'.format(tEnd-time[0],tStart-tEnd))
    IWindow    = np.where((time>=tStart) & (time<=tEnd) & (~np.isnan(time)))[0]
    iEnd   = IWindow[-1]
    iStart = IWindow[0]
    ## Absolute and relative differences at window extremities
    DeltaValuesAbs=(df.iloc[iEnd]-df.iloc[iStart]).abs()
#         DeltaValuesRel=(df.iloc[iEnd]-df.iloc[iStart]).abs()/df.iloc[iEnd]
    DeltaValuesRel=(df.iloc[IWindow].max()-df.iloc[IWindow].min())/df.iloc[IWindow].mean()
    #EndValues=df.iloc[iEnd]
    #if avgMethod.lower()=='periods_omega':
    #    if DeltaValuesRel['RotSpeed_[rpm]']*100>5:
    #        print('[WARN] Rotational speed vary more than 5% in averaging window ({}%) for simulation: {}'.format(DeltaValuesRel['RotSpeed_[rpm]']*100,f))
    ## Stats values during window
    # MeanValues = df[IWindow].mean()
    # StdValues  = df[IWindow].std()
    if 'mean' in stats:
        MeanValues = pd.DataFrame(df.iloc[IWindow].mean()).transpose()
    else:
        raise NotImplementedError()
    return MeanValues



def averagePostPro(outFiles,avgMethod='periods',avgParam=None,ColMap=None,ColKeep=None,ColSort=None,stats=['mean']):
    """ Opens a list of FAST output files, perform average of its signals and return a panda dataframe
    For now, the scripts only computes the mean within a time window which may be a constant or a time that is a function of the rotational speed (see `avgMethod`).
    The script only computes the mean for now. Other stats will be added

    `ColMap` :  dictionary where the key is the new column name, and v the old column name.
                Default: None, output is not sorted
                NOTE: the mapping is done before sorting and `ColKeep` is applied
                ColMap = {'WS':Wind1VelX_[m/s], 'RPM': 'RotSpeed_[rpm]'}
    `ColKeep` : List of strings corresponding to the signals to analyse. 
                Default: None, all columns are analysed
                Example: ColKeep=['RotSpeed_[rpm]','BldPitch1_[deg]','RtAeroCp_[-]']
                     or: ColKeep=list(ColMap.keys())
    `avgMethod` : string defining the method used to determine the extent of the averaging window:
                - 'periods': use a number of periods(`avgParam`), determined by the azimuth. 
                - 'periods_omega': use a number of periods(`avgParam`), determined by the mean RPM
                - 'constantwindow': the averaging window is constant (defined by `avgParam`).
    `avgParam`: based on `avgMethod` it is either
                - for 'periods_*': the number of revolutions for the window. 
                   Default: None, as many period as possible are used
                - for 'constantwindow': the number of seconds for the window
                   Default: None, full simulation length is used
    """
    result=None
    invalidFiles =[]
    # Loop trough files and populate result
    for i,f in enumerate(outFiles):
        try:
            df=FASTOutputFile(f).toDataFrame()
        except:
            invalidFiles.append(f)
            continue
        postpro=averageDF(df, avgMethod=avgMethod, avgParam=avgParam, ColMap=ColMap, ColKeep=ColKeep,ColSort=ColSort,stats=stats)
        MeanValues=postpro # todo
        if result is None:
            # We create a dataframe here, now that we know the colums
            columns = MeanValues.columns
            result = pd.DataFrame(np.nan, index=np.arange(len(outFiles)), columns=columns)
        result.iloc[i,:] = MeanValues.copy().values

    if ColSort is not None:
        # Sorting 
        result.sort_values([ColSort],inplace=True,ascending=True)
        result.reset_index(drop=True,inplace=True) 

    if len(invalidFiles)==len(outFiles):
        raise Exception('None of the files can be read (or exist)!')
    elif len(invalidFiles)>0:
        print('[WARN] There were {} missing/invalid files: {}'.format(len(invalidFiles),invalidFiles))


    return result 


if __name__ == '__main__':
    main()
