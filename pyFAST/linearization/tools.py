""" 
Simple tools for linearization analyses with OpenFAST
"""
import numpy as np
import re
import os
import glob

try:
    import pyFAST.linearization.mbc.mbc3 as mbc
except ImportError:
    import weis.control.mbc.mbc3 as mbc


def MBC_OF(fstfile, verbose=False, **kwargs):
    """ 
    Run necessary MBC for an OpenFAST file (one operating point)
    """
    fullpathbase, ext = os.path.splitext(fstfile)
    filebase = os.path.basename(fullpathbase)



    # --- Find available lin files
    # NOTE: the code below is problematic for module lin files ED.1.lin 
    # So we do a re search to filter these out
    # First use glob
    lin_file_fmt    = '{}.*.lin'.format(fullpathbase)
    lin_files       = glob.glob(lin_file_fmt)
    # Then use re for stricter search
    lin_file_fmt_re    = r'.*\.[0-9]+\.lin'
    lin_files = glob_re(lin_file_fmt_re, lin_files)

    # --- Check if checkpoint file exists
    fullpath_modes = None
    fullpath_chkp  = fullpathbase+ '.ModeShapeVTK.chkp'
    if os.path.exists(fullpath_chkp):
        fullpath_modes = fullpathbase+ '.ModeShapeVTK.pyPostMBC'

    # --- run MBC3 and campbell post_pro on lin files, generate postMBC file if needed 
    if len(lin_files)>0:
        if verbose:
            print('       Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
        MBC, matData, FAST_linData = mbc.fx_mbc3(lin_files, modesFilename=fullpath_modes, verbose=False)
    else:
        if verbose:
            print('[WARN] Lin. files: {} ({})'.format(lin_file_fmt, len(lin_files)))
        MBC, matData, FAST_linData = None, None, None



    # --- Write Viz file if checkpoint file present
    if fullpath_chkp is not None:
        vizfile = writeVizFile(fstfile, **kwargs)
        MBC['vizfile'] = vizfile
    else:
        MBC['vizfile'] = None


    return MBC, matData, FAST_linData


def writeVizFile(fstfile, VTKLinModes=15, VTKLinScale=10, VTKLinTim=1, VTKLinTimes1=True, VTKLinPhase=0, VTKModes=None):
    fullpathbase, ext = os.path.splitext(fstfile)
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






def glob_re(pattern, strings):
    return list(filter(re.compile(pattern).match, strings))
