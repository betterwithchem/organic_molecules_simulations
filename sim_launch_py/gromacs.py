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

        def add_cv(self, name: str ,cvtype: str , **kwargs):

            """ Add a collective variable to the simulation

            :param name: name of the collective variable
            :type name: str
            :param cvtype: type of collective variable, available types are: 'TORSION', 'ENERGY'.
            :type cvtype: str
            :param **kwargs : keyword arguments for the chosen type of variable
            :returns new_cv: the new collective variable
            :rtype new_cv: CV object 

            """
        
            supported=['TORSION','ENERGY']

            if cvtype.upper()=='TORSION':
                    new_cv=plumed.Torsion(name,kwargs['atoms'])
                    self._cvs.append(new_cv)
            elif cvtype.upper()=='ENERGY':
                    new_cv=plumed.PotentialEnergy(name)
                    self._cvs.append(new_cv)
            else:
                    print("Error: for the moment only {} are supported as collective variables... sorry".format(supported))
                    exit()

            return new_cv
    
        def add_bias(self, name: str, biastype: str, cv: object ,**kwargs):

                """ Add a bias to the simulation

                :param name: name of the bias
                :type name: str
                :param biastype: type of bias to apply, available types are:  'METAD', 'UPPER_WALLS', 'LOWER_WALLS'
                :type biastype: str
                :param cv: name of the collective variable to bias
                :type cv: object
                :param **kwargs : keyword arguments that depend on the specific bias
                :returns new_bias: the new bias object
                :rtype object: 

                """
                
                supported=["METAD","UPPER_WALLS","LOWER_WALLS"]
        
                if biastype.upper()=='METAD':
                    new_bias=plumed.Metad(name,cv,**kwargs)            
                elif biastype.upper()=="UPPER_WALLS":
                    new_bias=plumed.UpperWalls(name,cv,**kwargs)
                elif biastype.upper()=="LOWER_WALLS":
                    new_bias=plumed.LowerWalls(name,cv,**kwargs)                            
                else:
                    print("Error: for the moment only {} are supported as bias... sorry".format(supported))
                    exit()

                self._biases.append(new_bias)

                return new_bias
                

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

                


            
        
        


