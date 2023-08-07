from sim_launch_py.classes import Project
from sim_launch_py.ffparms import gaff, getTop
import numpy as np
import pandas as pd
import os
import sim_launch_py.utilities as util


molecules=pd.read_csv('molecules.list',sep='\s+',header=0)
systems=pd.read_csv('systems.list',sep='\s+',header=1)

ppath='./two_conc'

##### create a new project
project=Project.new_project(name='two_conc',path=ppath,overwrite=True)


##### save it
#project.save()

##### load existing project
#project=Project.load_project(ppath)

##### add species that will be used
 
# how to import the list of molecules to be used in the project is up to the user
# it could be a pandas dataframe from a text file or a dict or lists...
# this is done in order to keep a degree of flexibility in the user experience


for i,name in enumerate(molecules.molname):
    project.add_molecule(name=name,
                         resname=molecules.loc[i].at['resname'],
                         structure=molecules.loc[i].at['path']+'/'+name+'.pdb')


##### this will create also include topology (itp) files and a file with all atom types
for i,mol in enumerate(project.molecules):
    # if necessary, compute force field parameters
    #gaff(mol, os.path.abspath(project.topology_path),
    #     res_name=mol.resname, generate_charges='bcc', atomtype='gaff2',
    #     overwrite=False)

    # or copy them from a defined location
    getTop(mol,fromPath="/home/ucecmpa/Scratch/organic_molecules_simulations/Topologies" ,toPath=project.topology_path)

    # in any case compute molecular weight and number of atoms
    mol.mw=util.molecularWeightFromTop(mol.topology_path)
    mol.natoms=util.numberOfAtomsFromTop(mol.topology_path)
   



##### add systems to be simulated


for i,name in enumerate(systems.mol_1):
    project.add_system(name='s{}'.format(i))

for i,sys in enumerate(project.systems):
    sys.temperature=systems.loc[i].at['temperature']
    sys.add_molecule(systems.loc[i].at['mol_1'],moltype='solvent',knownmolecules=project.molecules)
    sys.add_molecule(systems.loc[i].at['mol_2'],moltype='solute',knownmolecules=project.molecules)
    sys.box=systems.loc[i].at['side']



n=0

for i,sys in enumerate(project.systems):

    #################################
    ## all the computation inside this for loop could be included in a method for Gromacs simulations
    ## something to think about...
    #################################
    
    # build initial configurations:

    # in this example we have a solute in a solvent. We first create a box and fill it with
    # solvent molecules at a given concentration
    # then we add solute molecules at a given concentration and remove overlapping solvent molecules
    # (both done by exploiting gmx tools)
    
    for mol in sys.molecules:
        if 'solvent' in mol.mol_attributes:
            solvent=mol
            print(sys.gromacs)
            sys.createSolventBox(solvent,output_structure='solvent_box.pdb',density=systems.loc[i].at['conc_1'])

    for mol in sys.molecules:
        if 'solute' in mol.mol_attributes:
            solute=mol
            sys.insertSolute(solute,solvent,solvent_box='solvent_box.pdb',concentration=systems.loc[i].at['conc_2'],output_structure='start.pdb')

    # get final number of molecules of each species in order to create the topology file:
    for mol in sys.molecules:
        mol.nmols=util.check_number_molecules(sys.path+'/start.pdb',mol)
    
    # write topology
    sys.writeTop(project.topology_path+'/atomtypes.itp',solvent,solute)

    # we can now create the files needed for the simulations:
    # for unbiased simulations we just need mdp files

project.save()

#project=Project.load_project(ppath)

project.job_script_path='sim_launch_py/job_scripts'
mdpdir='sim_launch_py/mdp/'

for sys in project.systems:

    sys.new_simulation('em',name='em',mdrun_options='-v -nsteps 500',start_coord=sys.path+'/start.pdb', mdp=mdpdir+'em.mdp',gmxbin=project.gromacs)
    sys.new_simulation('md',name='npt',mdrun_options='-v -nsteps 100000',mdp=mdpdir+'mdvvberendsen.mdp',maxwarn=1, gmxbin=project.gromacs)
    sys.new_simulation('md',name='md',mdrun_options='-v -nsteps 10000000',mdp=mdpdir+'mdvvparrinello.mdp', gmxbin=project.gromacs)

    #sys.print_command('run.sh')
        
project.write_sub_command('launch_jobs.sh',system='myriad')


    
    
    

