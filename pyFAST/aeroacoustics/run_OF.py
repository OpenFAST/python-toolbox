import os
import platform
import subprocess

FAST_exe        = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'openfast' + os.sep + 'build' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'openfast'
FAST_directory  = os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath(__file__) ) ) ) ) + os.sep + 'openfast' + os.sep + 'reg_tests' + os.sep + 'r-test' + os.sep + 'glue-codes' + os.sep + 'openfast' + os.sep + 'IEA_LB_RWT-AeroAcoustics'
FAST_InputFile  = FAST_directory + os.sep + 'IEA_LB_RWT-AeroAcoustics.fst'

if platform.system() == 'Windows':
    FAST_exe = FAST_exe + '.exe'

exec_str = []
exec_str.append(FAST_exe)
exec_str.append(FAST_InputFile)

subprocess.call(exec_str)