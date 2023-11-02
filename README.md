
# organic_molecules_simulations

## Installation
Needed python packages can be found in environment.yaml. On top of these, when creating a project, 
the code will look for gromacs and ambertools in $PATH.

A typical installation procedure would be:

` conda create -n myenv python=3.8 numpy pandas scipy jupyter conda-forge::parmed conda-forge::mdanalysis conda-forge::matplotlib conda-forge::rdkit conda-forge::nglview 
 conda activate myenv`

Then package can be installed locally in a conda environment running the following command from the directory where the `setup.py` file is  

`pip -e install .  `

## Code and Documentation

Source code is in `sim_launch_py/`  
You can access up-to-date html documentation opening `docs/build/hml/index.html` from an internet browser  

## Examples

Jupyter notebooks can be found in `examples/notebooks`

## Known Issues

Sometimes it happen that Python throws the following error when importing numpy:

`ImportError: Error importing numpy: you should not try to import numpy from
        its source directory; please exit the numpy source tree, and relaunch
        your python interpreter from there.`  

This is due to a modification of PYTHONPATH when sourcing AmberTools. Running

`export PYTHONPATH="" `

solves it.

