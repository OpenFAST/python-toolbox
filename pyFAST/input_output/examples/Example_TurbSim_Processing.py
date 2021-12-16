from pyFAST.input_output import FASTInputFile
from pyFAST.case_generation import runner
from pyFAST.input_output import TurbSimFile
import pandas as pd

def main():
    """Modify TurbSim parameters and write"""
    filename = '../tests/example_files/FASTIn_TurbSim.inp' # Name of TurbSim's input file
    f = FASTInputFile(filename)
    f['WrBLFF'] = False
    f['WrADFF'] = True
    f['WrADTWR'] = True
    f['NumGrid_Z'] = 15
    f['NumGrid_Y'] = 15
    f['AnalysisTime'] = 600
    f.write('../tests/example_files/FASTIn_TurbSim_change.inp')  # Modified file name

    """Run TurbSim"""
    TurbSim_FILE = '../tests/example_files/FASTIn_TurbSim_change.inp'  # Input file
    Turbsim_EXE = 'TurbSim_x64.exe'  # Change to the path of the TurbSim executable
    runner.run_cmd(TurbSim_FILE, Turbsim_EXE, wait=True, showOutputs=False, showCommand=True)

    """Open the turbulence box, containing the wind speed in 3 directions"""
    ts = TurbSimFile('../tests/example_files/FASTIn_TurbSim_change.bts')  # Output file
    print(ts.keys())
    print(ts['info'])
    print(ts['u'].shape)

    """Save wind speed data in CSV format"""
    fgp_data = pd.DataFrame(ts['u'][:,:,0,0]).T   # 3D turbulence data for the first grid point
    fgp_data.columns = ['streamwise wind speed', 'transverse wind speed', 'vertical wind speed'] # add columns name
    fgp_data.to_csv('../tests/example_files/wind_data_of_first_grid_point.csv')  # save data in CSV format

if __name__ == "__main__":
    main()
