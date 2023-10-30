from sim_launch_py.classes import Project
from sim_launch_py.external import molecules as mols
from sim_launch_py.plumed import write_plumed_file as wrtplumed
import numpy as np

pname='sample_metad'
p=Project.new_project(name=pname,path='./{}'.format(pname),overwrite=True)


solutes={"sulfamerazine":{"pdb":"/home/matteo/Work/organic_molecules_simulations/Structures/sulfamerazine/sulfamerazine.pdb", "top":"/home/matteo/Work/organic_molecules_simulations/Topologies/sulfamerazine.top","resname":"SUR"},
         "sulfadiazine":{"pdb":"/home/matteo/Work/organic_molecules_simulations/Structures/sulfadiazine/sulfadiazine.pdb", "top":"/home/matteo/Work/organic_molecules_simulations/Topologies/sulfadiazine.top","resname":"SUD"},
         "ccdc-xxxii":{"pdb":"/home/matteo/Work/organic_molecules_simulations/Structures/ccdc-xxxii/ccdc-xxxii.pdb", "top":"/home/matteo/Work/organic_molecules_simulations/Topologies/ccdc-xxxii.top","resname":"CCD"}}

solvents={"water":{"pdb":"/home/matteo/Work/organic_molecules_simulations/Structures/water/water.pdb", "top":"/home/matteo/Work/organic_molecules_simulations/Topologies/water.top", "density":1000., "MW":18., "resname":"WAT"},
          "1-butanol":{"pdb":"/home/matteo/Work/organic_molecules_simulations/Structures/1-butanol/1-butanol.pdb", "top":"/home/matteo/Work/organic_molecules_simulations/Topologies/1-butanol.top", "density":810.,"MW":74.123, "resname":"BNL"},
          "toluene":{"pdb":"/home/matteo/Work/organic_molecules_simulations/Structures/toluene/toluene.pdb", "top":"/home/matteo/Work/organic_molecules_simulations/Topologies/toluene.top", "density":870.,"MW":92.141,"resname":"TOL"}}

side=3.5 # nm
V=side*side*side*0.97

for solute in solutes:
    solutes[solute]["rotatables"]=mols.findTorsionalAngles(solutes[solute]['pdb'])

for solute in solutes:
    for solvent in solvents:

        name='{}_{}'.format(solute,solvent)
        p.add_system(name)

        s=p.find_system_by_name(name)

        s.add_box(side,shape='cubic')
        s.insert_molecules(solute,solutes[solute]['pdb'], nmol=1, final_conf='{}.pdb'.format(solute))

        s.center_box()
        
        nmol_solv= int(solvents[solvent]["density"]/solvents[solvent]["MW"]*6.022*1e23/1e24*V)

        s.insert_molecules(solvent,solvents[solvent]['pdb'],nmol=nmol_solv, initial_conf=s.last_saved_structure, final_conf='{}.pdb'.format(name))

        s.species[solutes[solute]["resname"]]['top']=solutes[solute]['top']
        s.species[solvents[solvent]["resname"]]['top']=solvents[solvent]['top']

        s.create_topology('topol.top')

        em_dict={'path_mdp':'/home/matteo/Work/organic_molecules_simulations/sim_launch_py/mdp/em.mdp',
                 'maxwarn':2, 'nsteps':1000, 'coordinates':s.last_saved_structure}

        s.add_simulation('em','em',simulation_dict=em_dict)

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
                 'nsteps':5000000, 'plumed':'plumed.dat'}
                 
        s.add_simulation('md','md',simulation_dict=md_dict)

        sim=s.simulations[s.find_simulation_by_name('md')]

        for irot,rot in enumerate(solutes[solute]['rotatables']):
            sim.add_cv('dih_{}'.format(irot),'torsion',cv_dict={'atoms':np.array(rot)+1})

        wrtplumed('{}/plumed.dat'.format(sim.path),sim,colvar='COLVAR')
        
        myriad_dict={'wallclock':'12:00:00',
                     'job_name':s.name,
                     'mpi':12,
                     'omp':6,
                     'gmx_bin':'gmx_mpi'}

        s.create_run_script('run.myriad',platform='myriad',platform_dict=myriad_dict)

        
        
        
p.save()        
        
