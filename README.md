
# organic_molecules_simulations

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

## Updates

##### Last update 20/07/2026

##### Added

- Added the possibility to solvate the box with 3-point water using the gmx solvate tool. This is available through System.solvate(initial_conf: str, final_conf:str). To this end, the directories structures/solvents has been added to the library with inside a structure of a single 3-point water molecule (water.pdb) that is used by the code. This can also be used as an example/guideline to expand the code to use other solvents (e.g. 4-point water). 

##### 05/05/2026

##### Fixed

- Bug in cell rotation that was introducing a misalignment between the coordinates of the atoms and the box, affecting PBCs in structures that span the whole box (e.g. bulk crystals or membranes)

#### 30/04/2026

##### Added

- `System.rotate_cell(angles: list, recenter=True, rotation_center=None, degrees=True)` to arbitrarily rotate the cell. Final orientation will have the `a` component aligned along the `x` axis.
- `System.write_ndx(filename:str = 'index.ndx', overwrite=False)` to create GROMACS compatible index (.ndx) files from the groups defined in the system.

##### 23/04/2026

##### Bug fixes:

- `System.delete_molecule()` had an issue in the numbering of the indexes
- `System._update_composition()` now removes from the `System.species` `dict` species that don't have any molecule with that resname (e.g. following deletion or renaming)
- Fixed bug that was not allowing to update the absolute indexes of the atoms, with repercussions of groups creation.

##### Added:

- `System.rename_molecules(idx_list:list, new_name_list: list)` to easily rename selected `Molecule` objects.
- `System.reorder_molecules()` to reorder the molecules in the `System`.


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
