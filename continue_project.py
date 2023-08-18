from sim_launch_py.classes import Project
import sim_launch_py.plumed as plumed
from sim_launch_py.external.gaussianFit import gaussianFit as gaussfit
import os
import numpy as np


"""In this example we load an existing project after running an unbiased MD simulation
   and estimate parameters to run a Metadynamics simulation of Torsional angles. Then
   we prepare the inputs for it.
"""

ppath='./single_mol'
project=Project.load_project(ppath)

project.job_script_path=os.path.abspath('sim_launch_py/job_scripts')
mdpdir=os.path.abspath('sim_launch_py/mdp/')

for isys,sys in enumerate(project.systems):

    # reset bash command for the system
    sys.bash_command=""
    
    # select the last MD simulation that has been performed
    md=sys.simulations[-1]

    # create the MetaD simulation
    metad=sys.add_simulation('metad','md',mdrun_options='-v -nsteps 10000000',mdp="{}/mdparrinello.mdp".format(mdpdir), plumed="metad.dat")
    
    # data file to be used to estimate the parameters
    colvar="{}/COLVAR".format(md.path)

    icol=0
    iangle=0
    
    for cv in md.cvs:
        if cv.cvtype=='TORSION':

            centers,sigma=gaussfit(colvar,datacol=icol+1)

            sigma_cv=min(sigma)/2

            metad.add_cv('dih_{}'.format(iangle),'torsion',atoms=cv.atoms)
            metad.add_bias('metad_dih_{}'.format(iangle), 'metad', 'dih_{}'.format(iangle),
                           sigma=np.min(sigma_cv), height=1.5, temp=300, pace=500, hills_file='HILLS_dih{}'.format(iangle),
                           biasfactor=10, grid_min='-pi', grid_max='pi', grid_spacing=min(2.5*np.pi/180,np.min(sigma_cv)/2))

            iangle+=1
            icol+=1

    plumed.writePlumedFile("{}/metad.dat".format(metad.path),metad,colvar="METAD_COLVAR",printstride=50)
    sys.setSimsToRun([sys.simulations[-1]])

project.save()

project.write_sub_command('launch_metad.sh',system='myriad')

    
