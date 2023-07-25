from sim_launch_py.classes import Project
from sim_launch_py.ffparms import gaff
import numpy as np
import pandas as pd
import os
import sim_launch_py.utilities as util

##### create a new project
project=Project.new_project(name='test',path='./new_test_project',overwrite=True)

##### save it
#project.save()

##### load existing project
#project=Project.load_project("./new_test_project")

##### add species that will be used
 
# how to import the list of molecules to be used in the project is up to the user
# it could be a pandas dataframe from a text file or a dict or lists...
# this is done in order to keep a degree of flexibility in the user experience

molecules=pd.read_csv('somemolecules.list',sep='\s+',header=0)
for i,name in enumerate(molecules.molname):
    project.add_molecule(name=name,
                         resname=molecules.loc[i].at['resname'],
                         structure=molecules.loc[i].at['path']+'/'+name+'.pdb')

##### if necessary, compute force field parameters
##### this will create also include topology (itp) files and a file with all atom types
for i,mol in enumerate(project.molecules):
    #gaff(mol.name, mol.structure_path, project.topology_path,
    #     res_name=mol.resname, generate_charges='bcc', atomtype='gaff2',
    #     overwrite=True)
    gaff(mol, project.topology_path,
         res_name=mol.resname, generate_charges='bcc', atomtype='gaff2',
         overwrite=True)    
    mol.mw=util.molecularWeightFromTop(mol.topology_path)


#project.save()
#exit()
#project=Project.load_project("./new_test_project")


##### add systems to be simulated
systems=pd.read_csv('systems.list',sep='\s+',header=1)

for i,name in enumerate(systems.mol_1):
    project.add_system(name='s{}'.format(i))

for i,sys in enumerate(project.systems):
    sys.temperature=systems.loc[i].at['temperature']
    sys.add_molecule(systems.loc[i].at['mol_1'],moltype='solvent',knownmolecules=project.molecules)
    sys.add_molecule(systems.loc[i].at['mol_2'],moltype='solute',knownmolecules=project.molecules)
    sys.box=5

project.save()
#project=Project.load_project("./new_test_project")

exit()
n=0
for i,sys in enumerate(project.systems):

    print(sys.path)

    for mol in sys.molecules:
        if 'solvent' in mol.mol_attributes:
            sys.createSolventBox(mol,output_structure='solvent_box.pdb',density=1000)

    for mol in sys.molecules:
        if 'solute' in mol.mol_attributes:
            sys.insertSolute(mol,solvent_box='solvent_box.pdb',concentration=systems.loc[i].at['conc_2'],output_structure='start.pdb')
            

    



    # add solute removing solvent

    # write topology

project.save()
#project=Project.load_project("./new_test_project")



#for sys in project.systems:
#    print(sys.name,':',sys.temperature,[mol.name for mol in sys.molecules],sys.box)
    

