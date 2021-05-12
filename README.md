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



