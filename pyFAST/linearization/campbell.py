""" 
Generic functions to help setup a Campbell diagram

def postproMBC(xlsFile=None, csvBase=None, sortedSuffix=None, csvModesID=None, xlssheet=None):
  Generate Cambell diagram data from an xls file, or a set of csv files
 

"""


import os, glob
import pandas as pd
import numpy as np
try:
    import pyFAST.linearization.mbc.mbc3 as mbc
except ImportError:
    import weis.control.mbc.mbc3 as mbc


def postproCampbell(out_or_fstfiles, BladeLen=None, TowerLen=None, verbose=True, nFreqOut=15, WS_legacy=None, removeTwrAzimuth=False):
    """ 
    Postprocess linearization files to extract Campbell diagram (linearization at different Operating points)
    - Run MBC
    - Postprocess to put into "CampbellData" matlab form
    - Perform mode identification (work in progress)
    - Export to disk
    """ 
    if len(out_or_fstfiles)==0:
        raise Exception('postproCampbell requires a list of at least one .fst or .out file')

    if len(out_or_fstfiles)==1 and os.path.splitext(out_or_fstfiles[0])[1]=='.pkl':
        import pickle
        CD = pickle.load(open(out_or_fstfiles[0],'rb'))
    else:
        # --- Run MBC for all operating points
        MBC = run_pyMBC(out_or_fstfiles, verbose, removeTwrAzimuth=removeTwrAzimuth)

        # --- Attemps to extract Blade Length and TowerLen from first file...
        filebase, ext = os.path.splitext(out_or_fstfiles[0])
        if BladeLen is None:
            fstFile = filebase+'.fst'
            if os.path.exists(fstFile):
                # try to read BladeLen and TowerLen from fst file
                # TODO: can be done with pyFAST.io or AeroElasticSE
                # The interface is very similar
                from pyFAST.input_output.fast_input_deck import FASTInputDeck
                fst = FASTInputDeck(fstFile, 'ED')
                ED = fst.fst_vt['ElastoDyn']
                if ED is None:
                    raise Exception('Unable to infer BladeLen and TowerLen, ElastoDyn file not found.')
                BladeLen = ED['TipRad'] - ED['HubRad']
                TowerLen = ED['TowerHt']
            else:
                raise Exception('Provide BladeLen and TowerLen, or existing fst/ED file')

        # --- Transform data into "CampbellData" (similar to matlab)
        CD = [None]*len(MBC)
        for i_lin, mbc_ in enumerate(MBC):
            if mbc_ is not None:
                CD[i_lin] = mbc.campbell_diagram_data(mbc_, BladeLen, TowerLen)

        CD = [cd for cd in CD if cd is not None]
        if len(CD)==0:
            raise Exception('No linearization file found')

    # --- Identify modes
    modeID_table,modesDesc=mbc.IdentifyModes(CD)

    # --- Write files to disk
    # Write csv file for manual identification step..
    baseName    = os.path.join(os.path.dirname(out_or_fstfiles[0]), 'Campbell')
    modeID_file = campbellData2CSV(baseName, CD, modeID_table, modesDesc)
    # Write summary txt file to help manual identification step..
    txtFileName = baseName+'_Summary.txt'
    campbellData2TXT(CD, nFreqOut=nFreqOut, txtFileName=txtFileName)

    # --- Return nice dataframes (assuming the identification is correct)
    # TODO, for now we reread the files...
    OP, Freq, Damp, UnMapped, ModeData = postproMBC(csvModesIDFile=modeID_file, verbose=False, WS_legacy=WS_legacy)
    
    return OP, Freq, Damp, UnMapped, ModeData, modeID_file

def run_pyMBC(out_or_fstfiles, verbose=True, removeTwrAzimuth=False):
    """
    Run MBC transform on set of openfast linear outputs

    INPUTS:
      - out_or_fstfiles
    """
    import re
    def glob_re(pattern, strings):
        return list(filter(re.compile(pattern).match, strings))

    if verbose:
        print('run_pyMBC:')
    MBC = [None]*len(out_or_fstfiles)
    for i_lin, fstfile in enumerate(out_or_fstfiles):
        filebase, ext = os.path.splitext(fstfile)
        # NOTE: the code below is problematic for module lin files ED.1.lin 
        # So we do a re search to filter these out
        # First use glob
        lin_file_fmt    = '{}.*.lin'.format(filebase)
        lin_files       = glob.glob(lin_file_fmt)
        #print(lin_files)
        # Then use re for stricter search
        lin_file_fmt_re    = r'.*\.[0-9]+\.lin'
        lin_files = glob_re(lin_file_fmt_re, lin_files)
        #print(lin_files)

        # --- run MBC3 and campbell post_pro on lin files 
        if len(lin_files)>0:
            if verbose:
                print('       Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
            MBC[i_lin], matData, FAST_linData = mbc.fx_mbc3(lin_files, verbose=False, removeTwrAzimuth=removeTwrAzimuth)
        else:
            if verbose:
                print('[WARN] Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))

    return MBC



def campbellData2TXT(CD, nFreqOut=15, txtFileName=None, skipHighDamp=True, skipNonEDBD=True):
    """ Write frequencies, damping, and mode contents for each operating points to a string
    Write to file if filename provided
    """
    if not isinstance(CD, list):
        CD=[CD]
    txt=''
    for iOP,cd in enumerate(CD):
        WS  = cd['WindSpeed']
        RPM = cd['RotSpeed_rpm']
        nFreqOut_loc = np.min([len(cd['Modes']),nFreqOut])
        txt+='------------------------------------------------------------------------\n'
        txt+='--- OP {:d} - WS {:.1f} - RPM {:.2f} \n'.format(iOP+1, WS, RPM)
        txt+='------------------------------------------------------------------------\n'
        skippedDamp=[]
        skippedEDBD=[]
        for im in np.arange(nFreqOut_loc):
            m    = cd['Modes'][im]
            Desc = cd['ShortModeDescr'][im]
            zeta = m['DampingRatio']
            hasED = Desc.find('ED')>=0
            hasBD = Desc.find('BD')>=0
            hasAD = Desc.find('AD')>=0
            if skipHighDamp and (zeta>0.96 or abs(zeta)<1e-5):
                skippedDamp.append(im)
            elif skipNonEDBD and (not (hasBD or hasED)):
                skippedEDBD.append(im)
            else:
                txt+='{:02d} ; {:8.3f} ; {:7.4f} ; {:s}\n'.format(im+1,m['NaturalFreq_Hz'],m['DampingRatio'],Desc)
        if len(skippedEDBD)>0:
            txt+='---- Skipped (No ED/BD)\n'
            for im in skippedEDBD:
                m    = cd['Modes'][im]
                Desc = cd['ShortModeDescr'][im]
                zeta = m['DampingRatio']
                txt+='{:02d} ; {:8.3f} ; {:7.4f} ; {:s}\n'.format(im+1,m['NaturalFreq_Hz'],m['DampingRatio'],Desc)
        if len(skippedDamp)>0:
            txt+='---- Skipped (High Damping)\n'
            for im in skippedDamp:
                m    = cd['Modes'][im]
                Desc = cd['ShortModeDescr'][im]
                zeta = m['DampingRatio']
                txt+='{:02d} ; {:8.3f} ; {:7.4f} ; {:s}\n'.format(im+1,m['NaturalFreq_Hz'],m['DampingRatio'],Desc)
    if txtFileName is not None:
        with open(txtFileName, 'w') as f:
            f.write(txt)
    return txt

def campbellData2CSV(baseName, CD, modeID_table, modesDesc):
    # Write summary of Campbell data and modes identification to several CSV files
    # The files generated will be:
    #    - [BaseName '_ModesID.csv ']
    #    - [BaseName '_OP.csv      ']
    #    - [BaseName '_PointsXX.csv'] for each operating point where XX is the index.
    #
    # INPUTS:
    #   - BaseName    :  basename that will be used to create the different CSV Files
    #   - ModesData: structure as returned by IdentifyModes

    nModes = modeID_table.shape[0]
    nOP    = modeID_table.shape[1]
    filenames=[]

    # --- Write ModeID using Matlab format
    filename='{:s}_ModesID.csv'.format(baseName)
    filenames.append(filename)
    with open(filename, 'w') as f:
        f.write('Mode Number Table,' +','.join(['']*nOP) +'\n')
        if np.isnan(CD[0]['WindSpeed']):
            f.write('Rotor Speed (rpm),' +','.join([str(cd['RotSpeed_rpm']) for cd in CD]) +'\n')
        else:
            f.write('Wind Speed (mps),' +','.join([str(cd['WindSpeed']) for cd in CD]) +'\n')
        for im, v in enumerate(modesDesc):
            f.write(v[0]+',' +','.join([str(ID) for ID in modeID_table[im,:]]) +'\n')

    # --- Write OP using Matlab format
    filename='{:s}_OP.csv'.format(baseName)
    filenames.append(filename)
    with open(filename, 'w') as f:
        f.write('Operating Points,'  +','.join(['']*nOP) +'\n')
        f.write('Wind Speed (mps),'  +','.join([str(cd['WindSpeed']) for cd in CD]) +'\n')
        f.write('Rotor Speed (rpm),' +','.join([str(cd['RotSpeed_rpm']) for cd in CD]) +'\n')

    # --- Write Modes for each OP
    for iOP,cd in enumerate(CD):
        filename='{:s}_Point{:02d}.csv'.format(baseName,iOP+1)
        filenames.append(filename)
        with open(filename, 'w') as f:
            Modes = cd['Modes']
            f.write(','.join(['Mode Number:,{:d},,,'.format(im+1) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Natural (undamped) frequency (Hz):, {:f},,,'.format(m['NaturalFreq_Hz']) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Damped frequency (Hz):, {:f},,,'.format(m['DampedFreq_Hz']) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Damping Ratio (-):, {:f},,,'.format(m['DampingRatio']) for im,m in enumerate(Modes)]) +'\n')
            f.write(','.join(['Mode {:d} state description,State has max at mode {:d},Mode {:d} signed magnitude,Mode {:d} phase (deg),'.format(im+1,im+1,im+1,im+1) for im,m in enumerate(Modes)]) +'\n')
            nComp = len(Modes[0]['DescStates'])
            for iC in np.arange(nComp):
                f.write(','.join(['{:s},{:d},{:f},{:f},'.format(m['DescStates'][iC].replace(',',' '),m['StateHasMaxAtThisMode'][iC],m['MagnitudePhase'][iC],m['PhaseDiff'][iC]) for im,m in enumerate(Modes)]) +'\n')
                   
    return filenames[0]


def postproMBC(xlsFile=None, csvModesIDFile=None, xlssheet=None, verbose=True, WS_legacy=None, suffix=''):
    """ 
    Generate Cambell diagram data from an xls file, or a set of csv files
    INPUTS:
      - xlsFile: path to an excel file, or, basename for a set of csv files generated by campbellFolderPostPro
      - csvModesIDFile: filename of csv file containing mode identification table.
           Its parent directory is noted csvBase.
           Other csv files will be asssumed to located in the same folder:
             - csvBase + 'Campbell_OP.csv'
             - csvBase + 'Campbell_PointI.csv' I=1..n_Op
             - (csvBase + 'Campbell_ModesID.csv')
      - vebose: more outputs to screen
      - WS_legacy: for OpenFAST 2.3, wind speed is unknown (NaN), a vector of operating point WS is needed
    OUTPUTS:
        - Freq: dataframe with columns [WS, RPM, Freq_Mode1,.., Freq_ModeN] for the N identified modes
        - Damp: dataframe with columns [WS, RPM, Damp_Mode1,.., Damp_ModeN] (damping ratios), for the N identified modes
        - UnMapped: dataframe with columns [WS, RPM, Freq, Damp] for all unidenfied modes
        - ModesData: low-level data, dictionaries for each OP with frequencies and damping
    """
    rpmSweep=False
    if xlsFile is not None:
        from pandas import ExcelFile
        # --- Excel file reading
        OPFileName=xlsFile;
        IDFileName=xlsFile;
        sheets=dict()
        # Reading all sheets
        if verbose:
            print('Reading Excel file: ',xlsFile)
        xls = pd.ExcelFile(xlsFile)
        dfs = {}
        for sheet_name in xls.sheet_names:
            # Reading sheet
            df = xls.parse(sheet_name, header=None)
            if df.shape[0]>0:
                sheets[sheet_name]=df
        OP   = sheets['OP']
        WS   = OP.iloc[1,1:].values
        RPM  = OP.iloc[2,1:].values
        ID   = None
        if any([s.find('mps')>0 for s in sheets.keys()]):
            rpmSweep  = False
            sweepVar  = WS
            sweepUnit = 'mps'
            xlssheet='WS_ModesID' if xlssheet is None else xlssheet
        else:
            rpmSweep  = True
            sweepVar  = RPM
            sweepUnit = 'rpm'
            xlssheet ='ModesID' if xlssheet is None else xlssheet
        if xlssheet not in sheets.keys():
            raise Exception('Mode identification sheet {} not found in excel file'.format(xlssheet))
        ID = sheets[xlssheet]
        # Storing data for each points, we try a bunch of keys since matlab script uses a loose num2str for now
        Points=dict()
        for i,v in enumerate(sweepVar):
            keys=[('{:.'+str(ires)+'f} {:s}').format(v,sweepUnit) for ires in [0,1,2,3,4]]
            for k in keys:
                try:
                    Points[i] = sheets[k]
                except:
                    pass
            if i not in Points.keys():
                raise Exception('Couldnf find sheet for operating point {:d}'.format(i))
    elif csvModesIDFile is not None:
        # --- csv file reading
        IDFileName=csvModesIDFile
        csvBase=os.path.join(os.path.dirname(csvModesIDFile),'')
        OPFileName=csvBase+'Campbell_OP{:}.csv'.format(suffix)
        if verbose:
            print('Reading csv file: ',OPFileName)
        OP      = pd.read_csv(OPFileName, sep = ',')
        if verbose:
            print('Reading csv file: ',IDFileName)
        ID      = pd.read_csv(IDFileName, sep = ',',header=None)
        nCol    = OP.shape[1]-1
        naCount = OP.isna().sum(axis=1)
        WS     = OP.iloc[0,1:].values
        RPM    = OP.iloc[1,1:].values
        del OP
        # Storing data for each points into a dict
        Points=dict()
        for i,v in enumerate(WS):
            OPFile = csvBase+'Campbell_Point{:02d}{:}.csv'.format(i+1,suffix)
            #print(OPFile, WS[i], RPM[i])
            Points[i] = pd.read_csv(OPFile, sep = ',', header=None)
    else:
        raise Exception('Provide either an Excel file or a csv (ModesID) file')
    # --- Mode Identification
    ID.iloc[:,0].fillna('Unknown', inplace=True) # replace nan
    ModeNames = ID.iloc[2: ,0].values
    ModeIDs   = ID.iloc[2: ,1:].values
    nModesIDd = len(ModeNames)  # Number of modes identified in the ID file

    if ModeIDs.shape[1]!=len(WS):
        print('OP Windspeed:',WS)
        raise Exception('Inconsistent number of operating points between OP ({} points) and ID ({} points) data.\nOP filename: {}\nID filename: {}\n'.format(ModeIDs.shape[1], len(WS), OPFileName, IDFileName))

    # --- Extract Frequencies and Damping from Point table
    ModeData=[]
    ioff=0
    coff=0
    for i,ws in enumerate(WS):
        P = Points[i]
        opData = dict()
        opData['Fnat']  = P.iloc[1+ioff, 1::5+coff].values[:].astype(float) # natural frequencies
        opData['Fdmp']  = P.iloc[2+ioff, 1::5+coff].values[:].astype(float) # damped frequencies
        opData['Damps'] = P.iloc[3+ioff, 1::5+coff].values[:].astype(float) # Damping values
        ModeData.append(opData)

    # ---Dealing with missing WS
    if np.any(np.isnan(np.array(WS).astype(float))) or np.all(np.array(WS).astype(float)==0):
        print('[WARN] WS were not provided in linearization (all NaN), likely old OpenFAST version ')
        if WS_legacy is not None:
            print('    > Using WS_legacy provided.')
            if len(WS_legacy)!=len(WS):
                raise Exception('WS_legacy should have length {} instead of {}'.format(len(WS),len(WS_legacy)))
            WS = WS_legacy
        else:
            print('[WARN] WS was replaced with index!')
            WS = np.arange(0,len(WS))

    # --- Creating a cleaner table of operating points
    OP = pd.DataFrame(np.nan, index=np.arange(len(WS)), columns=['WS_[m/s]', 'RotSpeed_[rpm]'])
    OP['WS_[m/s]']        = WS
    OP['RotSpeed_[rpm]'] = RPM

    UnMapped_WS   = []
    UnMapped_RPM  = []
    UnMapped_Freq = []
    UnMapped_Damp = []
    # --- Unidentified modes, before "nModes"
    for iOP,(ws,rpm) in enumerate(zip(WS,RPM)):
        m = ModeData[iOP]
        nModesMax = len(m['Fnat'])
        #nModes = min(nModesMax, nModesIDd) # somehow sometimes we have 15 modes IDd but only 14 in the ModeData..
        Indices = (np.asarray(ModeIDs[:,iOP])-1).astype(int)
        IndicesMissing = [i for i in np.arange(nModesMax) if i not in Indices]
        f   = np.asarray([m['Fnat'] [iiMode] for iiMode in IndicesMissing])
        d   = np.asarray([m['Damps'][iiMode] for iiMode in IndicesMissing])
        ws  = np.asarray([ws]*len(f))
        rpm = np.asarray([rpm]*len(f))
        UnMapped_Freq = np.concatenate((UnMapped_Freq, f))
        UnMapped_Damp = np.concatenate((UnMapped_Damp, d))
        UnMapped_WS   = np.concatenate((UnMapped_WS, ws))
        UnMapped_RPM  = np.concatenate((UnMapped_RPM, rpm))
    ## --- Unidentified modes, beyond "nModes"
    #for m,ws,rpm in zip(ModeData,WS,RPM):
    #   f   = m['Fnat'][nModesIDd:]
    #   d   = m['Damps'][nModesIDd:]
    #   ws  = np.asarray([ws]*len(f))
    #   rpm = np.asarray([rpm]*len(f))
    #   UnMapped_Freq = np.concatenate((UnMapped_Freq, f))
    #   UnMapped_Damp = np.concatenate((UnMapped_Damp, d))
    #   UnMapped_WS   = np.concatenate((UnMapped_WS, ws))
    #   UnMapped_RPM  = np.concatenate((UnMapped_RPM, rpm))

    # --- Put identified modes into a more convenient form
    cols=[ m.split('-')[0].strip().replace(' ','_') for m in ModeNames]
    cols=[v + str(cols[:i].count(v) + 1) if cols.count(v) > 1 else v for i, v in enumerate(cols)]
    Freq = pd.DataFrame(np.nan, index=np.arange(len(WS)), columns=cols)
    Damp = pd.DataFrame(np.nan, index=np.arange(len(WS)), columns=cols)
    for iMode in np.arange(nModesIDd):
        ModeIndices = (np.asarray(ModeIDs[iMode,:])-1).astype(int)
        ModeName= ModeNames[iMode].replace('_',' ')
        if ModeName.find('(not shown)')>0:
            f   = np.asarray([m['Fnat'][iiMode]  for m,iiMode in zip(ModeData,ModeIndices) if iiMode>=0 ])
            d   = np.asarray([m['Damps'][iiMode] for m,iiMode in zip(ModeData,ModeIndices) if iiMode>=0 ])
            ws  = np.asarray([ws  for ws,iiMode in zip(WS,ModeIndices)   if iiMode>=0 ])
            rpm = np.asarray([rpm for rpm,iiMode in zip(RPM,ModeIndices) if iiMode>=0 ])
            UnMapped_Freq = np.concatenate((UnMapped_Freq, f))
            UnMapped_Damp = np.concatenate((UnMapped_Damp, d))
            UnMapped_WS   = np.concatenate((UnMapped_WS, ws))
            UnMapped_RPM  = np.concatenate((UnMapped_RPM, rpm))
        else:
            if all(ModeIndices==-1):
                print('Skipping mode number ',iMode)
            else:
                f=[]
                d=[]
                #f2=np.asarray([m['Fnat'] [iiMode] if iiMode>=0 else np.nan for m,iiMode in zip(ModeData,ModeIndices)])
                #d2=np.asarray([m['Damps'][iiMode] if iiMode>=0 else np.nan for m,iiMode in zip(ModeData,ModeIndices)])
                for iOP, (m,iiMode) in enumerate(zip(ModeData, ModeIndices)):
                    if iiMode<0:
                        f.append(np.nan)
                        d.append(np.nan)
                    elif iiMode>=len(m['Fnat']):
                        print('[WARN] ID {} for mode `{}` at OP {} is beyond maximum number allowed'.format(iiMode+1, ModeName, iOP+1))
                        f.append(np.nan)
                        d.append(np.nan)
                    else:
                        f.append(m['Fnat'] [iiMode])
                        d.append(m['Damps'][iiMode])
                Freq.iloc[:, iMode]=f
                Damp.iloc[:, iMode]=d
    #  Removing modes that are full nan (not_shown ones)
    # NOTE: damgerous since OP is not part of it anymore
    # Freq.dropna(how='all',axis=0,inplace=True)
    # Freq.dropna(how='all',axis=1,inplace=True)
    # Damp.dropna(how='all',axis=0,inplace=True)
    # Damp.dropna(how='all',axis=1,inplace=True)

    # --- UnMapped modes into a dataframe
    M = np.column_stack((UnMapped_WS, UnMapped_RPM, UnMapped_Freq, UnMapped_Damp))
    UnMapped = pd.DataFrame(data=M, columns=['WS_[m/s]','RotSpeed_[rpm]','Freq_[Hz]','Damping_[-]'])

    return OP, Freq, Damp, UnMapped, ModeData


def campbellModeStyles(i, lbl):
    """ """
    import matplotlib.pyplot as plt
    FullLineStyles = [':', '-', '-+', '-o', '-^', '-s', '--x', '--d', '-.', '-v', '-+', ':o', ':^', ':s', ':x', ':d', ':.', '--','--+','--o','--^','--s','--x','--d','--.'];
    Markers    = ['', '+', 'o', '^', 's', 'd', 'x', '.']
    LineStyles = ['-', ':', '-.', '--'];
    Colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    MW_Light_Blue    = np.array([114,147,203])/255.
    MW_Light_Orange  = np.array([225,151,76])/255.
    MW_Light_Green   = np.array([132,186,91])/255.
    MW_LightLight_Green   = np.array([163,230,112])/255.
    MW_Light_Red     = np.array([226,115,115])/255.
    MW_Light_Gray    = np.array([128,133,133])/255.
    MW_Light_Purple  = np.array([144,103,167])/255.
    MW_Light_DarkRed = np.array([171,104,87])/255.
    MW_Light_Kaki    = np.array([204,194,16])/255.
    MW_Blue     =     np.array([57,106,177])/255.
    MW_Orange   =     np.array([218,124,48])/255.
    MW_Green    =     np.array([62,150,81])/255.
    MW_Red      =     np.array([204,37,41])/255.
    MW_Gray     =     np.array([83,81,84])/255.
    MW_Purple   =     np.array([107,76,154])/255.
    MW_DarkRed  =     np.array([146,36,40])/255.
    MW_Kaki     =     np.array([148,139,61])/255.

    lbl=lbl.lower().replace('_',' ')
    ms = 4
    c  = Colors[np.mod(i,len(Colors))]
    ls = LineStyles[np.mod(int(i/len(Markers)),len(LineStyles))]
    mk = Markers[np.mod(i,len(Markers))]
    # Color
    if any([s in lbl for s in ['1st tower']]):
        c=MW_Blue
    elif any([s in lbl for s in ['2nd tower']]):
        c=MW_Light_Blue
    elif any([s in lbl for s in ['1st blade edge','drivetrain']]):
        c=MW_Red
    elif any([s in lbl for s in ['1st blade flap']]):
        c=MW_Green
    elif any([s in lbl for s in ['2nd blade flap']]):
        c=MW_Light_Green
    elif any([s in lbl for s in ['3rd blade flap']]):
        c=MW_LightLight_Green
    elif any([s in lbl for s in ['2nd blade edge']]):
        c=MW_Light_Red
    elif any([s in lbl for s in ['1st blade torsion']]):
        c=MW_Purple
    # Line style
    if any([s in lbl for s in ['tower fa','collective','drivetrain','coll']]):
        ls='-'
    elif any([s in lbl for s in ['tower ss','regressive','bw']]):
        ls='--'
    elif any([s in lbl for s in ['tower ss','progressive','fw']]):
        ls='-.'
    # Marker
    if any([s in lbl for s in ['collective','coll']]):
        mk='2'; ms=8
    elif any([s in lbl for s in ['blade','tower','drivetrain']]):
        mk=''; 
    return c, ls, ms, mk

def plotCampbell(OP, Freq, Damp, sx='WS_[m/s]', UnMapped=None, fig=None, axes=None, ylim=None, legend=True, plotUnMapped=True, ps=[1,3,6,9]):
    """ Plot Campbell data as returned by postproMBC 

    INPUTS:
      - OP : dataframe of operating points (as returned by postMBC)
      - Freq : dataframe of Frequencies at each OP (as returned by postMBC)
      - OP : dataframe of damping ratios at each OP points (as returned by postMBC)
      - sx: label of the dataframes used for the "x" axis of the Campbell plot
    OPTIONAL INPUTS:
      - UnMapped: dataframe of UnMapped modes
      - fig, axes: optional fig and axes used for plotting (freq and damp)
      - ylim: limits for the frequency axis
      - ps: multiple of "p" (rotational speed) to plot in the background
    """
    import matplotlib.pyplot as plt


    # Init figure
    if fig is None:
        fig,axes_ = plt.subplots(1,2)
        # fig.set_size_inches(7,7.0,forward=True) # default is (6.4,4.8)
        fig.set_size_inches(13,7.0,forward=True) # default is (6.4,4.8)
        fig.subplots_adjust(top=0.78,bottom=0.11,left=0.04,right=0.98,hspace=0.06,wspace=0.16)
    if axes is None:
        axes=axes_

    # Estimating figure range
    FreqRange = [0                         , np.nanmax(Freq.values)*1.01]
    DampRange = [np.nanmin(Damp.iloc[:,2:]), np.nanmax(Damp.values)*1.01]

    if ylim is not None:
        FreqRange=ylim
    if DampRange[0]>0:
        DampRange[0]=0

    # Plot "background" 1p 3p 6p 9p. Potentiallymake this an option
    RPM     = OP['RotSpeed_[rpm]'].values
    omega   = RPM/60*2*np.pi
    freq_1p = omega/(2*np.pi)
    for p in ps:
        axes[0].plot(OP[sx].values,  p*freq_1p, ':',color=(0.7,0.7,0.7), lw=1.0)

    # Plot mapped modes
    Markers = ['', '+', 'o', '^', 's', 'd', 'x', '.']
    iModeValid=0
    xPlot=[]; yPlot=[]
    for iMode,lbl in enumerate(Freq.columns.values):
        if lbl.find('not_shown')>0:
            # TODO ADD TO UNMAPPED
            continue
        iModeValid+=1
        c, ls, ms, mk = campbellModeStyles(iModeValid, lbl)
        if len(RPM)==1 and len(mk)==0:
            mk = Markers[np.mod(iMode,len(Markers))]
            ms=5
        axes[0].plot(OP[sx].values, Freq[lbl].values, ls, marker=mk, label=lbl.replace('_',' '), markersize=ms, color=c)
        axes[1].plot(OP[sx].values, Damp[lbl].values, ls, marker=mk                            , markersize=ms, color=c)
        xPlot=np.concatenate((xPlot, OP[sx].values))
        yPlot=np.concatenate((yPlot, Freq[lbl].values))

    # Unmapped modes (NOTE: plotted after to over-plot)
    if plotUnMapped and UnMapped is not None:
        axes[0].plot(UnMapped[sx].values, UnMapped['Freq_[Hz]'  ].values, '.', markersize=2, color=[0.5,0.5,0.5])
        axes[1].plot(UnMapped[sx].values, UnMapped['Damping_[-]'].values, '.', markersize=1, color=[0.5,0.5,0.5])
    # Highligh duplicates (also after)
    Points=[(x,y) for x,y in zip(xPlot,yPlot)]
    Plot = pd.Series(Points)
    for xDupl,yDupl in Plot[Plot.duplicated()]:
        axes[0].plot(xDupl,yDupl, 'o',color='r')

    axes[0].set_xlabel(sx.replace('_',' '))
    axes[1].set_xlabel(sx.replace('_',' '))
    axes[0].set_ylabel('Frequencies [Hz]')
    axes[1].set_ylabel('Damping ratios [-]')
    if legend:
        axes[0].legend(bbox_to_anchor=(0., 1.02, 2.16, .802), loc='lower left', ncol=4, mode="expand", borderaxespad=0.)
    if not np.any(np.isnan(FreqRange)):
        axes[0].set_ylim(FreqRange)
    
    XLIM=axes[1].get_xlim()
    axes[1].plot(XLIM, [0,0],'-', color='k', lw=0.5)
    axes[1].set_xlim(XLIM)
    if not np.any(np.isnan(DampRange)):
        axes[1].set_ylim(DampRange)
    return fig, axes


def plotCampbellDataFile(xls_or_csv, ws_or_rpm='rpm', sheetname=None, ylim=None, WS_legacy=None, to_csv=False, suffix='', returnData=False, 
        fig=None, axes=None, legend=True, plotUnMapped=True):
    """ 
    Wrapper for plotCampbell, takes an Excel or csv file as argument. Returns a figure.

    INPUTS:


      - WS_legacy: for OpenFAST 2.3, wind speed is unknown (NaN), a vector of operating point WS is needed
    """
    if ws_or_rpm.lower()=='ws':
        sx='WS_[m/s]'
    else:
        sx='RotSpeed_[rpm]'

    ext     = os.path.splitext(xls_or_csv)[1].lower()
    baseDir = os.path.dirname(xls_or_csv)
    basename = os.path.splitext(os.path.basename(xls_or_csv))[0]

    # --- Read xlsx or csv filse
    if ext=='.xlsx':
        OP, Freq, Damp, UnMapped, ModeData =  postproMBC(xlsFile=xls_or_csv,xlssheet=sheetname, WS_legacy=WS_legacy, suffix=suffix)

    elif ext=='.csv':
        OP, Freq, Damp, UnMapped, ModeData =  postproMBC(csvModesIDFile=xls_or_csv, WS_legacy=WS_legacy, suffix=suffix)

        pass

    else:
        raise Exception('Extension should be csv or xlsx, got {} instead.'.format(ext),)

    # --- Plot
    fig, axes = plotCampbell(OP, Freq, Damp, sx=sx, UnMapped=UnMapped, ylim=ylim, fig=fig, axes=axes, legend=legend, plotUnMapped=plotUnMapped)
    figName = os.path.join(baseDir,basename+'_'+ws_or_rpm)

    if to_csv:
        Freq.to_csv(os.path.join(baseDir, 'freq.csv'), index=None, sep=' ', na_rep='NaN')
        Damp.to_csv(os.path.join(baseDir, 'damp.csv'), index=None, sep=' ', na_rep='NaN')
        OP.to_csv(os.path.join(baseDir, 'op.csv'), index=None, sep=' ', na_rep='NaN')

    if returnData:
        return fig, axes, figName, OP, Freq, Damp
    else:
        return fig, axes, figName

