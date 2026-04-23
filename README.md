
# organic_molecules_simulations

#### Last update 23/04/2026

##### Bug fixes:

`System.delete_molecule()` had an issue in the numbering of the indexes
`System._update_composition()` now removes from the `System.species` `dict` species that don't have any molecule with that resname (e.g. following deletion or renaming)

##### Added:

`System.rename_molecules(self, idx_list:list, new_name_list: list)` to easily rename selected `Molecule` objects.
`System.reorder_molecules()` to reorder the molecules in the `System`.

## Installation
Needed python packages can be found in environment.yaml. 
On top of these, when creating a project, 
the code will look for Gromacs and Ambertools in $PATH, for their installation refer to the respective online documentation.

A typical installation procedure would be:

```
conda create -n myenv python=3.12 mamba 

conda activate myenv 

mamba install -c conda-forge numpy pandas scipy jupyter parmed mdanalysis matplotlib rdkit nglview ipykernel  
```

Then package can be installed locally in a conda environment running the following command from the directory where the `setup.py` file is  

```
pip install -e .  
```

## Code and Documentation

Source code is in `sim_launch_py/`  
You can access up-to-date html documentation opening `docs/build/hml/index.html` from an internet browser  

## Examples

Jupyter notebooks can be found in `examples/notebooks`

## Known Issues

Sometimes it happens that Python throws the following error when importing numpy:

```
ImportError: Error importing numpy: you should not try to import numpy from
        its source directory; please exit the numpy source tree, and relaunch
        your python interpreter from there.  
```

This is due to a modification of PYTHONPATH when sourcing AmberTools. Running

```
export PYTHONPATH="" 
```


solves it.

