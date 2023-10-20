from sim_launch_py.classes import Project
from sim_launch_py.utilities import get_distance
import numpy as np
import os

# Create a new Project

pname='sample_project'
ppath='./{}'.format(pname)

p=Project.new_project(pname,ppath,overwrite=True)

# Create a system

p.add_system('UNOGIN_eq1_bulk')

s=p.systems[0]

s.add_molecule('olanzapine', structure_file='/home/matteo/Work/organic_molecules_simulations/UNOGIN_eq1.pdb')

s.replicate_cell(repl=[10,10,10])

s.save_pdb('replicated_cell.pdb')


# center of the box position
cob=np.array([0. for i in range(3)])

n=0
for m in s.molecules:
    for a in m.atoms:
            cob+=a.coordinates
            n+=1

cob/=n

# lets keep only the molecules that have COM within 2nm of the center of the box

delete_mols=[]
cutoff=2 #nm
for m in s.molecules:
    if get_distance(m.com,ref=cob,box=s.box)>cutoff:
        delete_mols.append(m.index)

s.delete_molecule(delete_mols)
s.save_pdb('remaining_mols.pdb')

s.add_box(6,shape='cubic')
s.center_box()
s.save_pdb('centered.pdb')

#get center of the seed:
seed_com=np.array([0., 0., 0.])

n=0
for m in s.molecules:
    for a in m.atoms:
        seed_com+=a.coordinates
        n+=1
    
seed_com/=n

s.insert_molecules('olanzapine','/home/matteo/Work/organic_molecules_simulations/Structures/olanzapine/olanzapine.pdb', initial_conf='centered.pdb', final_conf='added_olanzapine.pdb', nmol=6)

# estimate with some criterion the number of solvent molecules to add
nwat_to_add=int(1000/18*6.022/10*(s.box[0]**3))-300

# add solvent molecules
s.insert_molecules('water','/home/matteo/Work/organic_molecules_simulations/Structures/water/water.pdb',initial_conf='added_olanzapine.pdb',final_conf='with_water.pdb',nmol=nwat_to_add)

# once water molecules have been inserted, remove those that could be within the seed

water_mols=s.find_molecule_by_resname('WAT')

delete_mols=[]
for wm in s.molecules[water_mols[0]:water_mols[-1]]:
    if get_distance(wm.com,ref=seed_com, box=s.box)<=cutoff:
        delete_mols.append(wm.index)

s.delete_molecule(delete_mols)

# rinse and repeat if we need a precise number of solvent molecules. something like "while [nmols solvent < target nmols]; insert molecules; remove those that are in forbidden positions"

s.save_pdb('removed_water_in_seed.pdb')

s.species['OLA']['top']='/home/matteo/Work/organic_molecules_simulations/Topologies/olanzapine.top'
s.species['WAT']['top']='/home/matteo/Work/organic_molecules_simulations/Topologies/water.top'

# create a .top file, this will read the .top files of each species, extract the atomtypes and the molecule definitions and compose them in a single .top file.
s.create_topology('topol.top')

# for each simulation, define a dict with the options for the simulation (some parameters, if not defined, will use default values)

em_dict={'path_mdp':'/home/matteo/Work/organic_molecules_simulations/sim_launch_py/mdp/em.mdp',
         'maxwarn':2, 'nsteps':1000, 'coordinates':s.last_saved_structure}

s.add_simulation('em','em', simulation_dict=em_dict)

nvt_dict={'path_mdp': '/home/matteo/Work/organic_molecules_simulations/sim_launch_py/mdp/nvt.mdp',
          'maxwarn':2,
          'coordinates':'{}/{}.gro'.format(s.simulations[-1].path,s.simulations[-1].name),
          'posre':'{}/{}.gro'.format(s.simulations[-1].path,s.simulations[-1].name),
          'nsteps':50000}

s.add_simulation('nvt','posre',simulation_dict=nvt_dict)

npt_dict={'path_mdp': '/home/matteo/Work/organic_molecules_simulations/sim_launch_py/mdp/npt.mdp',
          'maxwarn':2,
          'coordinates':'{}/{}.gro'.format(s.simulations[-1].path,s.simulations[-1].name),
          'posre':'{}/{}.gro'.format(s.simulations[-1].path,s.simulations[-1].name),
          'nsteps':50000}

s.add_simulation('npt','posre',simulation_dict=npt_dict)

md_dict={'path_mdp': '/home/matteo/Work/organic_molecules_simulations/sim_launch_py/mdp/mdparrinello.mdp',
         'maxwarn':2, 'coordinates':'{}/{}.gro'.format(s.simulations[-1].path,s.simulations[-1].name),
         'nsteps':5000000}

s.add_simulation('md','md', simulation_dict=md_dict)

# options for myriad
platform_dict={'wallclock':'06:00:00',
               'job_name':s.name,
               'mpi':16,
               'omp':4}

s.create_run_script('run.job',platform='myriad',platform_dict=platform_dict)


                 

