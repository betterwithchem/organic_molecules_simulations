import shutil
import os
import sim_launch_py.plumed as plumed

class _Simulation():

        def __init__(self, name: str, mdrun_options: str, topology: str,
                 path: str):
            """Simulation Class Constructor

            Args:
                name (str): Name of the simulation
                mdrun_options (str): options for running MD simulation
                topology (str): name of the topology file
                path (str): path definition
            """
            self._name=name
            self._mdrun_options=mdrun_options
            self._topology=topology
            #self._run_options=run_options
            self._path=path

        @property
        def name(self):
            return self._name

        @property
        def mdrun_options(self):
            return self._mdrun_options

        @property
        def topology(self):
            return self._topology

        #@property
        #def run_options(self):
        #    return self._run_options

        @property
        def path(self):
            return self._path

        @name.setter
        def name(self,n):
            self._name=n


class MD(_Simulation):
    
    def __init__(self,  name: str, mdrun_options='', coord='', topology='',
                 path_mdp='', path_input='', path_output='', temperature=300, thermostat='vv',
                 pressure=1, barostat='Parrinello-Rahman', nsteps=100,
                 print_bash=False, maxwarn=0, gmxbin='gmx_mpi',
                 plumed=None):
        """MD simulation class constructor

        Args:
            name (str): Name of the system / project
            mdrun_options (str, optional): MD simulation options. Defaults to ''.
            coord (str, optional): coordinate file name. Defaults to ''.
            topology (str, optional): topology file name. Defaults to ''.
            path_mdp (str, optional): path to mdp. Defaults to ''.
            path_input (str, optional): path to inout file. Defaults to ''.
            path_output (str, optional): path to output file. Defaults to ''.
            temperature (int, optional): temperature of the MD simulation. Defaults to 300.
            thermostat (str, optional): thermostat. Defaults to 'vv'.
            pressure (int, optional): pressure. Defaults to 1.
            barostat (str, optional): barostat. Defaults to 'Parrinello-Rahman'.
            nsteps (int, optional): number of simulation steps. Defaults to 100.
            plumed (str, optional): plumed file. Defaults to None.
            print_bash (bool, optional): print a bash string. Defaults to False.
            maxwarn (int, optional): number of warnings tollerated. Defaults to 0.
            gmxbin (str, optional): gromacs executable name. Defaults to 'gmx_mpi'.
        """

        # Inherit properties of the _Simulation Class
        super().__init__(name, mdrun_options, topology, path_output)

        self._name=name
        self._cvs=[]
        self._biases=[]

        # copy mdp files
        shutil.copy(path_mdp,path_input)

        # Print bash string
        if print_bash:
            mdp=os.path.basename(path_mdp)

            grompp = "{0} grompp -f {1} -p {2} -c {3} -o {4}.tpr -maxwarn {5}".format(gmxbin,
                                                                                          mdp,
                                                                                          topology,
                                                                                          coord,
                                                                                          name,
                                                                                          maxwarn)
            mdrun = "{0} mdrun -deffnm {1} {2}".format(gmxbin,name,mdrun_options)
                       
            if plumed is not None:
                    mdrun +=" -plumed {}".format(plumed)

            self.bash_command="{}\n{}\n".format(grompp,mdrun)


    @property
    def cvs(self):
        return self._cvs

    @property
    def biases(self):
        return self._biases

    def add_cv(self,name,cvtype,**kwargs):
        
        if cvtype.upper()=='TORSION':
            new_cv=plumed.Torsion(name,kwargs['atoms'])
            self._cvs.append(new_cv)
        else:
            print("Error: for the moment only TORSION is supported as collective variable... sorry")
            exit()
    
    def add_bias(self,name,biastype,**kwargs):

        if biastype.upper()=='METAD':
            new_bias=plumed.Metad(name,**kwargs)
            self._biases.append(new_bias)
        else:
            print("Error: for the moment only METAD is supported as bias... sorry")
            exit()
                        
            
            

                


class EnergyMinimization(_Simulation):

    def __init__(self,name: str,
                 mdrun_options='', coord='start.pdb', topology='topol.top',
                 path_mdp='', path_input='',path_output='',
                 print_bash=False,maxwarn=0,gmxbin='gmx_mpi'):
        """Energy minimization class

        Args:
            name (str): Name of the system / project
            mdrun_options (str, optional): MD simulation options. Defaults to ''.
            coord (str, optional): coordinate file name. Defaults to ''.
            topology (str, optional): topology file name. Defaults to ''.
            path_mdp (str, optional): path to mdp. Defaults to ''.
            path_input (str, optional): path to inout file. Defaults to ''.
            path_output (str, optional): path to output file. Defaults to ''.
            print_bash (bool, optional): print a bash string. Defaults to False.
            maxwarn (int, optional): number of warnings tollerated. Defaults to 0.
            gmxbin (str, optional): gromacs executable name. Defaults to 'gmx_mpi'.
        """        
        # Inherit properties of the _Simulation Class
        super().__init__(name, mdrun_options, topology, path_output)

        self._name=name

        # copy mdp files
        shutil.copy(path_mdp,path_input)

        # Print bash string
        if print_bash:
            mdp=os.path.basename(path_mdp)
            self.bash_command="{5} grompp -f {0} -p {1} -c {2} -o {3}.tpr -maxwarn {4}\n\n".format(mdp,topology,coord,name,maxwarn,gmxbin)+\
                "{2} mdrun -deffnm {0} {1}\n".format(name,mdrun_options,gmxbin)
                
            

        



"""
class Metad(MD, _Simulation):

    def __init__(self,plumedfile='plumed.dat'):

        super().__init__()

        self._plumedfile=plumedfile

    @property
    def plumedfile(self):
        return self._plumedfile

    @plumedfile.setter
    def plumedfile(self,p):
        self._plumedfile=p
        
    def add_cv(self):

        newcv=CV()
        self.append(newcv)


    def add_metad(self):

        newmetad=Metad()
        self.append(newmetad)

    def writePlumedFile(self):

        with open(self.plumedfile,'w') as f:
            f.write("# do nothing")

    def writeRunFile(self,runfile):

        with open(runfile,'w') as f:
            f.write("gmx mdrun -deffnm md -v -plumed {0} ".format())
"""


            
        
        


