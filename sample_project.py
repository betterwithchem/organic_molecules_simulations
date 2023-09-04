from sim_launch_py.classes import Project
from sim_launch_py.ffparms import gaff, getTop
import numpy as np
import pandas as pd
import os
import sim_launch_py.utilities as util
import sim_launch_py.plumed as plumed

from sim_launch_py.external.molecules import findTorsionalAngles

molecules=pd.read_csv('molecules.list',sep='\s+',header=0)
systems=pd.read_csv('systems.list',sep='\s+',header=1)

pname='single_mol'
ppath='./{}'.format(pname)


##### create a new project
project=Project.new_project(pname,ppath,overwrite=True)


##### save it
#project.save()

##### load existing project
#project=Project.load_project(ppath)

##### add species that will be used
 
# how to import the list of molecules to be used in the project is up to the user
# it could be a pandas dataframe from a text file or a dict or lists...
# this is done in order to keep a degree of flexibility in the user experience


for i,name in enumerate(molecules.molname):
    project.add_molecule(name,
                         resname=molecules.loc[i].at['resname'],
                         structure=molecules.loc[i].at['path']+'/'+name+'.pdb')


##### this will create also include topology (itp) files and a file with all atom types
for i,mol in enumerate(project.molecules):

    # if necessary, compute force field parameters
    #gaff(mol, os.path.abspath(project.topology_path),
    #     res_name=mol.resname, generate_charges='bcc', atomtype='gaff2',
    #     overwrite=False)

    # or copy them from a defined location
    getTop(mol,fromPath=os.path.abspath("Topologies") ,toPath=project.topology_path)


##### add systems to be simulated

for i,name in enumerate(systems.mol_1):
    project.add_system('s{}'.format(i))

for i,sys in enumerate(project.systems):
    sys.temperature=systems.loc[i].at['temperature']
    sys.add_molecule(systems.loc[i].at['mol_1'],moltype='solvent',knownmolecules=project.molecules)
    sys.add_molecule(systems.loc[i].at['mol_2'],moltype='solute',knownmolecules=project.molecules)
    


for i,sys in enumerate(project.systems):
    
    # build initial configurations:

    # in this example we have a solute in a solvent. We first create a box, then we fill it with
    # solvent molecules at a given concentration. Finally, we add solute molecules at a given
    # concentration and remove overlapping solvent molecules

    sys.addBox(systems.loc[i].at['side'],shape='dodecahedron')
    
    for mol in sys.molecules:
        if 'solvent' in mol.mol_attributes:
            solvent=mol
            print(sys.gromacs)
            sys.createSolventBox(solvent,output_structure='{}/solvent_box.pdb'.format(sys.path),density=systems.loc[i].at['conc_1'])

    for mol in sys.molecules:
        if 'solute' in mol.mol_attributes:
            solute=mol
            sys.insertSolute(solute,solvent,solvent_box='{}/solvent_box.pdb'.format(sys.path),concentration=systems.loc[i].at['conc_2'],output_structure='{}/start.pdb'.format(sys.path))
    
    # write topology
    sys.writeTop(project.topology_path+'/atomtypes.itp') #,solvent,solute)

    # we can now create the files needed for the simulations:
    # for unbiased simulations we just need mdp files

#project.save()


#project=Project.load_project(ppath)


project.job_script_path=os.path.abspath('sim_launch_py/job_scripts')
mdpdir=os.path.abspath('sim_launch_py/mdp/')

for sys in project.systems:

    sys.add_simulation('em','em',mdrun_options='-v -nsteps 500',start_coord=sys.path+'/start.pdb', mdp="{}/em.mdp".format(mdpdir))
    sys.add_simulation('npt','md',mdrun_options='-v -nsteps 100000',mdp="{}/mdvvberendsen.mdp".format(mdpdir),maxwarn=1)

    md=sys.add_simulation('md','md',mdrun_options='-v -nsteps 10000000',mdp="{}/mdparrinello.mdp".format(mdpdir), plumed="plumed.dat" )

    mol=sys.molecules[-1]
    dihangles=findTorsionalAngles(mol.structure_path)

    md.add_cv('ene', 'energy')
    
    for iangle,angle in enumerate(dihangles):
        md.add_cv('dih_{}'.format(iangle), 'torsion', atoms=[mol.atoms[i].atomID for i in angle])   

        md.add_bias('metad_dih_{}'.format(iangle),'metad','dih_{}'.format(iangle),sigma=5*np.pi/180,
                    height=1.5, temp=300, pace=500, hills_file='HILLS_dih{}'.format(iangle),
                    biasfactor=5, grid_min='-pi', grid_max='pi', grid_spacing=2.5*np.pi/180)

    plumed.writePlumedFile("{}/plumed.dat".format(sys.path),md,colvar="COLVAR",printstride=50)

    sys.setSimsToRun(sys.simulations)
    

project.save()
    
project.write_sub_command('launch_jobs.sh',system='myriad')


    
    
    

