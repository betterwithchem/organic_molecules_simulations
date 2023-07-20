import os
import sys
from scripts import ffparms
from scripts import toolbox

import numpy as np
import pandas as pd

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-system',dest='systems_file',type=str, 
                    help='file with systems definitions')
parser.add_argument('-molecules',dest='molecules_file', type=str,
                    help='file with molecules definitions')
args = parser.parse_args()

# check for availability of needed commands

molecules_file=os.path.abspath(args.molecules_file) #'/home/matteo/Work/solvent_concentration/somemolecules.list'
systems_file=os.path.abspath(args.systems_file) #'/home/matteo/Work/solvent_concentration/systems.list'

# create new project and change directory to the project
project=toolbox.Project(title='test',path='test_project',overwrite=False, continuation=True)

# copy input files
toolbox.copyFiles(molecules_file,systems_file,path=project.path,overwrite=False)

# Read species 
species=pd.read_csv(molecules_file, header=0, sep='\s+')


# Add species to project and create local directories for structure and parameter files
project.addSpecies(species)

# compute gaff parameters if needed
for i,s in enumerate(project.species):
    print("\nNow creating parameters for molecule {0} ({1})".format(s,getattr(project,s).resname))
    ffparms.gaff(s, getattr(project,s).structure, getattr(project,s).path, res_name=getattr(project,s).resname, generate_charges='bcc', atomtype='gaff2', overwrite=False)
    getattr(project,s).top=getattr(project,s).path+'/'+s+'.top'


atomtypes={}
for imol,mol in enumerate(project.species):
    print("#### {} ####".format(imol))
    atomtypes=toolbox.extract_atomtypes_from_gmx_top(getattr(project,mol).top,atomtypes)

project.atomtypes=atomtypes

# create include (itp) files for each species in the respective directories
project.createIncludeFiles()

"""
## TO DO
for each species:
- create itp file with molecule definition
"""

# add system to project
systems=pd.read_csv(systems_file, header=1, sep='\s+')

project.addSystems(systems)

for s in project.systems:

    nmol_solute=50
    volume=toolbox.vol_from_molconc(nmol_solute,getattr(project,s).solute_conc,getattr(project,getattr(project,s).solute).mw)
    box_vector=float(volume**(1/3))
    
    getattr(project,s).addBox(box_vector, shape='cubic')

    nsolv= int(getattr(project,s).solvent_dens * 6.022/10 / getattr(project,getattr(project,s).solvent).mw * volume)

    getattr(project,getattr(project,s).solute).nmol=nmol_solute
    getattr(project,getattr(project,s).solvent).nmol=nsolv

    print(box_vector)
        
exit()

"""
## TO DO ##
for each system:

- create simulation box
- create topology (header + include files + footer)
- create mdp input files
"""



exit()

# create run input files


# run simulations



                                
                
