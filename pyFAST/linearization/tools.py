""" 
Simple tools to assist in doing linearization analyses with OpenFAST
"""
import numpy as np
import re
import os
import glob

from pyFAST.linearization.mbc import fx_mbc3
from pyFAST.linearization.campbell_data import campbell_diagram_data_oneOP # 


def getFST_and_LinFiles(fstFile_or_linFiles, verbose=False):
    """ 
    Given a .fst or a list of .lin files, return both:
    - if a fst file is provided, .lin file next to it are sought for
    - if a list of line files are provided, return the .fst file from which they originated
    """
    if isinstance(fstFile_or_linFiles,str):
        # The user provided a string, we expect it's a .fst file
        fullpathbase, ext = os.path.splitext(fstFile_or_linFiles)
        if ext.lower()!='.fst':
            raise Exception('Provide either one fst file, or a list of .lin files')
        fstFile = fstFile_or_linFiles
        # --- Find available lin files
        linFiles = findLinFiles(fstFile, verbose=verbose)
    else:
        # the user provided a list (hopefully)
        fullpathbase, ext = os.path.splitext(fstFile_or_linFiles[0])
        if ext.lower()!='.lin':
            raise Exception('Provide either one fst file, or a list of .lin files')
        linFiles = fstFile_or_linFiles
        fullpathbase, ext = os.path.splitext(fullpathbase)
        fstFile = fullpathbase+'.fst'

    return fstFile, linFiles


def getCDDOP(fstFile_or_linFiles, verbose=False, BladeLen=None, TowerLen=None, **kwargs):
    """ 
    Return Campbell Data at one operating point from a .fst file or a list of lin files
    """
    # --- Figure out if the user provided a .fst file or a list of .lin files
    fstFile,linFiles =  getFST_and_LinFiles(fstFile_or_linFiles, verbose=verbose)

    # --- Open lin files for given OP/fst file, perform MBC
    MBCOP, matData = getMBCOP(fstFile=fstFile, linFiles=linFiles, verbose=verbose, **kwargs)

    # --- Estimate blade and tower length for scaling
    if BladeLen is None and TowerLen is None:
        BladeLen, TowerLen = estimateLengths(fstFile)

    # --- put data into "CampbellData" format
    CDDOP = campbell_diagram_data_oneOP(MBCOP, BladeLen, TowerLen)

    return CDDOP, MBCOP

def getMBCOP(fstFile, linFiles=None, verbose=False, **kwargs):
    """ 
    Run necessary MBC for an OpenFAST file (one operating point)
    """

    # --- Find available lin files
    if linFiles is None:
        linFiles = findLinFiles(fstFile, verbose=verbose)

    # --- Check if checkpoint file exists
    fullpathbase, ext = os.path.splitext(fstFile)
    fullpath_modes = None
    fullpath_chkp  = fullpathbase+ '.ModeShapeVTK.chkp'
    if os.path.exists(fullpath_chkp):
        fullpath_modes = fullpathbase+ '.ModeShapeVTK.pyPostMBC'

    # --- run MBC3 and campbell post_pro on lin files, generate postMBC file if needed 
    if len(linFiles)>0:
        MBC, matData = fx_mbc3(linFiles, modesFilename=fullpath_modes, verbose=False)
    else:
        MBC, matData = None, None

    # --- Write Viz file if checkpoint file present
    if fullpath_chkp is not None:
        vizfile = writeVizFile(fstFile, **kwargs)
        MBC['vizfile'] = vizfile
    else:
        MBC['vizfile'] = None

    return MBC, matData


def writeVizFile(fstFile, VTKLinModes=15, VTKLinScale=10, VTKLinTim=1, VTKLinTimes1=True, VTKLinPhase=0, VTKModes=None):
    fullpathbase, ext = os.path.splitext(fstFile)
    filebase          = os.path.basename(fullpathbase)
    fullpath_viz      = fullpathbase + '.ModeShapeVTK.viz'
    base_chkp         = filebase     + '.ModeShapeVTK'
    base_modes        = filebase     + '.ModeShapeVTK.pyPostMBC'

    if VTKModes is None:
        VTKModes=','.join((np.arange(VTKLinModes)+1).astype(str))

    with open(fullpath_viz, 'w') as f:
        f.write('------- OpenFAST MODE-SHAPE INPUT FILE -------------------------------------------\n');
        f.write('# Options for visualizing mode shapes\n');
        f.write('---------------------- FILE NAMES ----------------------------------------------\n');
        f.write('"{:s}"   CheckpointRoot - Rootname of the checkpoint file written when OpenFAST generated the linearization files (without the ".chkp" extension)\n'.format(base_chkp))
        f.write('"{:s}"   ModesFileName - Name of the mode-shape file (with eigenvectors)\n'.format(base_modes))
        f.write('---------------------- VISUALIZATION OPTIONS -----------------------------------\n')
        f.write('{:d}        VTKLinModes   - Number of modes to visualize (0 <= VTKLinModes <= NumModes)\n'.format(VTKLinModes))
        f.write('{:s}        VTKModes      - List of which VTKLinModes modes will be visualized (modes will be added sequentially from the last value entered)\n'.format(VTKModes))
        f.write('{:f}        VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n'.format(VTKLinScale)) 
        f.write('{:d}        VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n'.format(VTKLinTim))
        f.write('{}      VTKLinTimes1  - If VTKLinTim=2, visualize modes at LinTimes(1) only? (if false, files will be generated at all LinTimes)\n'.format(VTKLinTimes1))
        f.write('{:f}       VTKLinPhase   - Phase used when making one animation for all LinTimes together (used only when VTKLinTim=1)\n'.format(VTKLinPhase))
#             if isnan(opts.VTKLinScale)
#                 % Then user didn't specify it, we use some logic
#                 if CompElast==1 % ElastoDyn - VTKLinScale=10
#                     fprintf(fid,'10    VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n');
#                 elseif CompElast==2 % BeamDyn - VTKLinScale=0.1
#                     fprintf(fid,'0.1   VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n'); 
#                     fprintf(fid,'%f    VTKLinScale   - Mode shape visualization scaling factor (exaggerates mode shapes: try 10 for ElastoDyn; 0.1 for BeamDyn)\n',opts.VTKLinScale); 
#             end
#             if isnan(opts.VTKLinTim)
#                 % The user didn't specify this, we use some logic
#                 fprintf(fid,'2         VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n');
#             else
#                 if (RPM<1e-3) % When RPM =0, VTKLinTim=1 would only produce one VTK
#                     fprintf(fid,'2         VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n');
#                 else
#                     fprintf(fid,'%d        VTKLinTim     - Switch to make one animation for all LinTimes together (VTKLinTim=1) or separate animations for each LinTimes (VTKLinTim=2)\n',opts.VTKLinTim);
#                 end
#             end
#             if length(opts.VTKLinTimes1)==0
#                 % The user didn't specify this
#                 fprintf(fid,'true      VTKLinTimes1  - If VTKLinTim=2, visualize modes at LinTimes(1) only? (if false, files will be generated at all LinTimes)\n');
#             else
#                 fprintf(fid,'%s        VTKLinTimes1  - If VTKLinTim=2, visualize modes at LinTimes(1) only? (if false, files will be generated at all LinTimes)\n',opts.VTKLinTimes1);
#             end
#             fclose(fid);
        print('Written: ',fullpath_viz)
    return fullpath_viz


# ------------------------------------------------------------------------
# def writeModes(f0, fd, zeta, Q, modesFilename, format='OFBinary'):
def writeModes(VTK, modesFilename):
    """ 
    write binary file that will be read by OpenFAST to export modes to VTK
    """
    import struct
    def fwrite(fid, data, type):
        """ Mimic the matlab function fwrite"""
        # @ is used for packing in native byte order
        #  B - unsigned integer 8 bits
        #  h - integer 16 bits
        #  i - integer 32 bits
        #  f - float 32 bits
        #  d - float 64 bits
        fmt, _ = {'uint8': ('B', 1), 'int16':('h', 2), 'int32':('i', 4), 'float32':('f', 4), 'float64':('d', 8)}[type]
        if hasattr(data, '__len__'):
            data = data.flatten(order='F')
            n=len(data)
            fid.write(struct.pack('@'+str(n)+fmt, *data))
        else:
            fid.write(struct.pack('@'+fmt, data))

    fileFmt = 'float64' #8-byte real numbers
    nStates, nModes, nLinTimes = VTK['x_eig_magnitude'].shape

    #------- HACK
    #VTK['NaturalFreq_Hz'] =VTK['NaturalFreq_Hz'][:19]*0 + 1
    #VTK['DampingRatio']   =VTK['DampingRatio']  [:19]*0 + 2
    #VTK['DampedFreq_Hz']  =VTK['DampedFreq_Hz'] [:19]*0 + 3
    #     for iMode in range(nModes):
    #         VTK['x_eig_magnitude'][:,iMode,:] = np.zeros((nStates,nLinTimes)) + iMode+1
    #         VTK['x_eig_phase']    [:,iMode,:] = np.zeros((nStates,nLinTimes)) + iMode+1
    #         VTK['x_eig_magnitude'][2,iMode,:] = 12
    #         VTK['x_eig_phase']    [4,iMode,:] = 11
    #nModes=1
    # ------END HACK
    # --- Reduce differences python/Matlab by rounding
    #VTK['NaturalFreq_Hz']  = np.round(VTK['NaturalFreq_Hz']  *1000)/1000
    #VTK['DampingRatio']    = np.round(VTK['DampingRatio']    *1000)/1000
    #VTK['DampedFreq_Hz']   = np.round(VTK['DampedFreq_Hz']   *1000)/1000
    #VTK['x_eig_magnitude'] = np.round(VTK['x_eig_magnitude'] *1000)/1000
    #VTK['x_eig_phase'    ] = np.round(VTK['x_eig_phase']     *1000)/1000

    # --- Write to disk
    with open(modesFilename, 'wb') as fid:
        fwrite(fid, 1,        'int32' )# write a file identifier in case we ever change this format
        fwrite(fid, nModes,   'int32' )# number of modes (for easier file reading)
        fwrite(fid, nStates,  'int32' )# number of states (for easier file reading)
        fwrite(fid, nLinTimes,'int32' )# number of azimuths (i.e., LinTimes) (for easier file reading)
        # Freq and damping (not used in the FAST visualization algorithm)
        fwrite(fid, VTK['NaturalFreq_Hz'], fileFmt)
        fwrite(fid, VTK['DampingRatio'],   fileFmt)
        fwrite(fid, VTK['DampedFreq_Hz'],  fileFmt)
        # Writing data mode by mode
        for iMode in range(nModes):
            fwrite(fid, VTK['x_eig_magnitude'][:,iMode,:], fileFmt)
            fwrite(fid, VTK['x_eig_phase']    [:,iMode,:], fileFmt)
    print('Written: ', modesFilename)



def findLinFiles(fstFile, verbose=False):
    """ 
    Find .lin files given a .fst file
    """
    fullpathbase, ext = os.path.splitext(fstFile)
    # NOTE: the code below is problematic for module lin files ED.1.lin 
    # So we do a re search to filter these out
    # First use glob
    lin_file_fmt    = '{}.*.lin'.format(fullpathbase)
    lin_files       = glob.glob(lin_file_fmt)
    # Then use re for stricter search
    lin_file_fmt_re    = r'.*\.[0-9]+\.lin'
    lin_files = glob_re(lin_file_fmt_re, lin_files)
    if verbose:
        print('       Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
    if verbose:
        print('[WARN] Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
    return lin_files


def estimateLengths(fstFile):
    if os.path.exists(fstFile):
        # try to read BladeLen and TowerLen from fst file
        # TODO: can be done with pyFAST.io or AeroElasticSE
        # The interface is very similar
        from pyFAST.input_output.fast_input_deck import FASTInputDeck
        fst = FASTInputDeck(fstFile, 'ED')
        print(fst)
        ED = fst.fst_vt['ElastoDyn']
        if ED is None:
            raise Exception('Unable to infer BladeLen and TowerLen, ElastoDyn file not found.')
        BladeLen = ED['TipRad'] - ED['HubRad']
        TowerLen = ED['TowerHt']
    else:
        raise Exception('Provide `BladeLen` and `TowerLen`, or, an existing fst and ED file')
    return BladeLen, TowerLen

def glob_re(pattern, strings):
    """ Apply a pattern to a list of strings 
    Typically used as "glob" and "re" to select a list of files matting a given mattern"""
    return list(filter(re.compile(pattern).match, strings))
