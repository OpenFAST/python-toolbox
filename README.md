# pyfast

[![Build status](https://github.com/openfast/python-toolbox/workflows/Development%20Pipeline/badge.svg)](https://github.com/OpenFAST/python-toolbox/actions?query=workflow%3A%22Development+Pipeline%22)
[![Python: 3.6+](https://img.shields.io/badge/python-3.6%2B-informational)](https://www.python.org/)

Python package to work with NREL-supported [OpenFAST](https://github.com/OpenFAST/openfast) tool.
This repository intends to provide simple scripts to help OpenFAST users setup models, run simulations and postprocess the results. 


## Installation and testing

```bash
git clone http://github.com/OpenFAST/python-toolbox
cd python-toolbox
python -m pip install -e .
pytest
```



## Subpackages

The repository contains a set of small packages:

- input\_output: read/write OpenFAST/FAST.Farm/OLAF input and output files (see [README](pyFAST/input_output)) and postprocess OpenFAST outputs (see [examples](pyFAST/input_output/examples))
- linearization: tools to deal with OpenFAST linearization, e.g. generate a Campbell diagram (see [examples](pyFAST/linearization/examples/))
- aeroacoustics: tools for aeroacoustics (generate BL files and plot outputs)
- case\_generation: tools to generate and run a set of input of OpenFAST input files (see [examples](pyFAST/case_generation/examples))


## QuickStart and main usage

### Read and write files
Find examples scripts in this [folder](pyFAST/input_output/examples) and the different fileformats [here](pyFAST/input_output). 

Read an AeroDyn file (or any OpenFAST input file), modifies some values and write the modified file:
```python
from pyFAST.input_output import FASTInputFile
filename = 'AeroDyn.dat'
f = FASTInputFile(filename)
f['TwrAero'] = True
f['AirDens'] = 1.225
f.write('AeroDyn_Changed.dat')
```

Read an OpenFAST binary output file and convert it to a pandas DataFrame
```python
from pyFAST.input_output import FASTOutputFile
df = FASTOutputFile('5MW.outb').toDataFrame()
time  = df['Time_[s]']
Omega = df['RotSpeed_[rpm]']
```

Read a TurbSim binary file, modify it and write it back
```python 
from pyFAST.input_output import TurbSimFile
ts = TurbSimFile('Turb.bts')
print(ts.keys())
print(ts['u'].shape)  
ts['u'][0,:,:,:] += 1 # Adding 1 m/s in the streamwise
tw.write('NewTurbulenceBox.bts')
```

### Polar/airfoil manipulation
Find examples scripts in this [folder](pyFAST/polar/examples).


Read a CSV file with `alpha, Cl, Cd, Cm`, and write it to AeroDyn format (also computes unsteady coefficients)
```python 
from pyFAST.airfoils.Polar import Polar
polar = Polar('pyFAST/airfoils/data/DU21_A17.csv', fformat='delimited')
ADpol = polar.toAeroDyn('AeroDyn_Polar_DU21_A17.dat')
```

### Write a set of OpenFAST input file for multiple simulations
Find examples scripts in this [folder](pyFAST/case_generation/examples).



## Future work and friend projects

This repository intends to provide simple scripts to help users of OpenFAST. 
The repo is still in its early phase, so you may find that functionalities are missing or not bullet proof. 
Your contributions would be much appreciated, feel free to post issues and pull-requests. We will thrive to provide tests and examples.

In the meantime, you can also find relevant python scripts in the following repositories:

- [WISDEM](https://github.com/WISDEM/WISDEM): models for assessing overall wind plant cost of energy (COE), also contains file IO, (DLC) case generation, polar manipulations, visualization, and much more! 
- [weio](https://github.com/ebranlard/weio) : reader and writer for typical files used by the wind energy community
- [welib](https://github.com/ebranlard/welib): misc tools for wind energy applications (BEM, FEM, polars, ...)
- [ROSCO_toolbox](https://github.com/NREL/ROSCO_toolbox): tools to work with the [ROSCO](https://github.com/NREL/ROSCO) controller that is supported by OpenFAST
- [windtools](https://github.com/NREL/windtools): toolsfor wind simulation setup, data processing and analysis
- [pyDatView](https://github.com/ebranlard/pyDatView): cross-platform visualization and processing program for input and output files

Matlab script are found in the [matlab-toolbox]((https://github.com/OpenFAST/matlab-toolbox)).

Open-source OpenFAST wind turbine models can be found here:
- [openfast-turbine-models](https://github.com/NREL/openfast-turbine-models): open source wind turbine models (in development)
- [r-test](https://github.com/OpenFAST/r-test): regression tests for OpenFAST, contains models for OpenFAST and its drivers (AeroDyn, SubDyn, HydroDyn, etc.). This repository is not intended to be used as a "database" of models, but it has the advantage that the input files are alwasy up to date with the latest [format specifications](https://openfast.readthedocs.io/en/master/source/user/api_change.html)


General documentation for OpenFAST is found on its [readthedocs](https://openfast.readthedocs.io/) page.



## Contributing

This repository is still work in progress, thank you for your understanding.
Any contribution is much appreciated, feel free to post issues or pull-requests.



