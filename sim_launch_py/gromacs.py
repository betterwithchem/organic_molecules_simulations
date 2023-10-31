import shutil
import os
import sim_launch_py.plumed as plumed

class _Simulation():

        def __init__(self, name: str, index: int):
                """Simulation Class Constructor

                :param name: Name of the simulation.
                :type name: str
                :param index: Index of the simulation.
                :type index: int 

                """
                self._name=name
                self._index=index
                self._path=None
                self._simtype=None
                self._state=None
                self._gromacs=None                
                self._mdrun_options=None
                # if running also grompp
                self._topology=None
                self._coordinates=None
                self._mdp=None
                self._maxwarn=None
                # if running directly mdrun
                self._tpr=None
                self._nsteps=None
                self._cvs=list()
                                
        @property
        def name(self):
                return self._name

        @property
        def index(self):
                return self._index

        @property
        def path(self):
                return self._path

        @property
        def state(self):
                return self._state

        @property
        def gromacs(self):
                return self._gromacs

        @property
        def mdrun_options(self):
                return self._mdrun_options

        @property
        def topology(self):
                return self._topology

        @property
        def coordinates(self):
                return self._coordinates

        @property
        def mdp(self):
                return self._mdp

        @property
        def maxwarn(self):
                return self._maxwarn

        @property
        def tpr(self):
                return self._tpr

        @property
        def simtype(self):
                return self._simtype

        @property
        def nsteps(self):
                return self._nsteps

        @property
        def cvs(self):
                return self._cvs

        @name.setter
        def name(self,n):
                self._name=n

        @index.setter
        def index(self,ndx):
                self._index=ndx

        @path.setter
        def path(self,p):
                self._path=p

        @state.setter
        def state(self,s):
                self._state=s

        @gromacs.setter
        def gromacs(self,gmx):
                self._gromacs=gmx

        @mdrun_options.setter
        def mdrun_options(self,run_opt):
                self._mdrun_options=run_opt

        @topology.setter
        def topology(self,top):
                self._topology=top

        @coordinates.setter
        def coordinates(self,x):
                self._coordinates=x

        @mdp.setter
        def mdp(self,m):
                self._mdp=m

        @maxwarn.setter
        def maxwarn(self,mw):
                self._maxwarn=mw

        @tpr.setter
        def tpr(self,t):
                self._tpr=t

        @simtype.setter
        def simtype(self,t):
                self._simtype=t

        @nsteps.setter
        def nsteps(self,n):
                self._nsteps=n

        def _print_gmx_command(self):
                """Create the list of commands used to preprocess (grompp) and run (mdrun) the simulation 

                :returns gmx_commands: Gromacs commands.
                :rtype gmx_commands: list

                """

                gmx_commands=[]
                if self.mdp:

                        gmx_commands.append('{0} grompp -f {1} -o {2}.tpr -maxwarn {3} -p {4} -c {5}'.format(self.gromacs,self.mdp,self.name,self.maxwarn,self.topology,self.coordinates))

                gmx_commands.append('{0} mdrun -deffnm {1} {2}'.format(self.gromacs,self.name,self.mdrun_options))
                
                return gmx_commands


        def add_cv(self, name: str, cv_type: str, cv_dict: dict=None):
                """Add a Collective Variable to the Simulation.

                :param name: Name of the new collective variable.
                :type name: str
                :param cv_type: Type of collective variable. See documentation of sim_launch_py.plumed for details.
                :type cv_type: str
                :param cv_dict: dict with the parameters for the collective variable. Defaults to None.
                :type cv_dict: dict, optional
                :returns: 

                """

                import sim_launch_py.plumed as plumed

                available_cv_types=['TORSION','ENERGY']

                if cv_type.upper() not in available_cv_types:
                    print("!"*20)
                    print("! ERROR: given cv type ({}) is not a valid type. Acceptable types are {}.".format(cv_type.upper(),available_cv_types))
                    print("!"*20)
                    return

                if cv_type.upper()=='TORSION':

                    new_cv=plumed.Torsion(name,cv_dict['atoms'])

                elif cv_type.upper()=='ENERGY':

                    new_cv=plumed.PotentialEnergy(name)

                self.cvs.append(new_cv)
        

        def add_bias(self, name: str, biastype: str, cv: object , bias_dict: dict=None):

                """ Add a bias to the simulation

                :param name: name of the bias
                :type name: str
                :param biastype: type of bias to apply, available types are:  'METAD'
                :type biastype: str
                :param cv: name of the collective variable to bias
                :type cv: object
                :param bias_dict : dict of parameters for the bias. Defaults to None.
                :type bias_dict: dict, optional.
                
                """
                
                supported=["METAD"] #,"UPPER_WALLS","LOWER_WALLS"]
        
                if biastype.upper()=='METAD':
                    new_bias=plumed.Metad(name,cv,bias_dict)            
                #elif biastype.upper()=="UPPER_WALLS":
                #    new_bias=plumed.UpperWalls(name,cv,bias_dict)
                #elif biastype.upper()=="LOWER_WALLS":
                #    new_bias=plumed.LowerWalls(name,cv,bias_dict)                            
                else:
                    print("Error: for the moment only {} are supported as bias... sorry".format(supported))
                    exit()

                self._biases.append(new_bias)

                #return new_bias
                

class MD(_Simulation):
        
        def __init__(self,  name: str, index: int):
                """MD simulation class constructor

                :param name: Name of the simulation.
                :type name: str
                :param index: Index of the simulation.
                :type index: int

                """
                # Inherit properties of the _Simulation Class
                super().__init__(name, index)

                self._cvs=[]
                self._biases=[]

                self._simtype='md'
        
        @property
        def cvs(self):
            return self._cvs

        @property
        def biases(self):
            return self._biases

            
class EnergyMinimization(_Simulation):

        def __init__(self, name: str, index: int):
                """Energy minimization class

                :param name: Name of the simulation.
                :type name: str
                :param index: Index of the simulation
                :type index: int

                """
                # Inherit properties of the _Simulation Class
                super().__init__(name, index)

                self._simtype='em'


class Posre_MD(_Simulation):

        def __init__(self, name: str, index: int):
                """Position Restrained Molecular Dynamics Simulation

                :param name: Name of the simulation.
                :type name: str
                :param index: Index of the simulation
                :type index: int

                """
                # Inherit properties of the _Simulation Class
                super().__init__(name, index)

                self._cvs=[]
                self._biases=[]
                
                self._simtype='posre'

        @property
        def cvs(self):
            return self._cvs

        @property
        def biases(self):
            return self._biases                

                


            
        
        


