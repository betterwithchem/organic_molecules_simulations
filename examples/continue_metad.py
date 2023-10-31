from sim_launch_py.classes import Project
from sim_launch_py.external.gaussianFit import gaussianFit as gfit
import numpy as np
from copy import deepcopy
from sim_launch_py.plumed import write_plumed_file as wrtplumed

p=Project.load_project('sample_metad')

s=p.systems[-1]
sim=s.simulations[-1]

ndih=len(sim.cvs)

centers,sigmas=[],[]

for col in range(ndih):

    cen,sig=gfit(sim.path+'/COLVAR',datacol=col+1,plotFit=False)

    centers.append(cen)
    sigmas.append(sig)
    
# only sigmas are important in this example:
# sigmas are used to define the SIGMA for metad simulations
    
metad_dict={'path_mdp': '/home/matteo/Work/organic_molecules_simulations/sim_launch_py/mdp/mdparrinello.mdp',
            'maxwarn':2, 'coordinates':'{}/{}.gro'.format(s.simulations[-1].path,s.simulations[-1].name),
            'nsteps':5000000, 'plumed':'plumed.dat'}

s.add_simulation('metad','md',simulation_dict=metad_dict)
sim=s.find_simulation_by_name('metad')

prev_sim=s.find_simulation_by_name('md')


for idih,dih_sigma in enumerate(sigmas):
    
    sig=np.round(min(dih_sigma)/2,decimals=3)

    print(idih,sig, sig*180/np.pi)

    # duplicate dihedral cvs used in MD simulation
    newcv=deepcopy(prev_sim.cvs[idih])
    sim.cvs.append(newcv)

    metadbias_dict={'sigma':sig,
                    'height':2.5,
                    'biasfactor':5,
                    'grid_min':'-pi',
                    'grid_max':'pi',
                    'grid_spacing':np.round(np.pi/180*2,decimals=4),
                    'pace':500,
                    'hills_file':'HILLS_dih_{}'.format(idih)}

    sim.add_bias('metad_dih_{}'.format(idih),'METAD',sim.cvs[-1].name,bias_dict=metadbias_dict)
    
wrtplumed('{}/plumed.dat'.format(sim.path),sim,colvar='COLVAR')

#p.save()
