


class _CV():

    def __init__(self, name: str, print_cv=True):
        """CV Class Constructur
   
        Args:
           name (str): Name of the collective variable
           print_cv (bool, optional): Print the value of the variable? Defaults to True.
        """
        
        self._name=name
        self._print_cv=print_cv

    @property
    def name(self):
        return self._name

    @property
    def print_cv(self):
        return self._print_cv

    @name.setter
    def name(self,n):
        self._name=n

    @print_cv.setter
    def print_cv(self,p):
        self._print_cv=p

class Torsion(_CV):

    def __init__(self,name: str, atoms: list()):
        """Torsion collective variable class constructor
 
        Args:
            name(str) : Name of the collective variable. 
            atoms(list) : List of atom ids that compose the dihedral angle. 
        """

        # Inherit properties of the _CV Class
        super().__init__(name)

        self._name=name
        
        if len(atoms)!=4:
            print("ERROR: Exactly 4 atoms are needed to define a Torsion angle. {} given {}.".format(len(atoms),atoms))
            exit()
        
        self._atoms=atoms

        self._directive="{}: TORSION ATOMS={},{},{},{}".format(name,
                                                              atoms[0],
                                                              atoms[1],
                                                              atoms[2],
                                                              atoms[3])

        

    @property
    def atoms(self):
        return self._atoms

    @property
    def directive(self):
        return self._directive


def writePlumedFile(plumed_file: str, simulation: object, colvar=None,printstride=500):

    """Write Plumed File
    
    Args:       
       simulation (simulation object): simulation with groups, cvs, and biases to use.
       plumed_file (str): path of the output plumed file.
colvar (str, optional): name of the colvar file where values of CVs will be saved. Defaults to None.
       printstride (int, optional): stride for output of colvar file. Defaults to 500 steps.
    """

    import os
    import shutil
    
    if os.path.isfile(plumed_file):
        location=os.path
        shutil.copy(plumed_file,"{}.bkp".format(plumed_file))
        
    with open(plumed_file,'w') as f:

        if hasattr(simulation,'groups'):
            f.write("# Groups section\n\n")
            for igroup,group in enumerate(simulation.groups):
                f.write("{}\n".format(group.directive))

        if hasattr(simulation,'cvs'):
            f.write("\n# CV section\n\n")
            for icv,cv in enumerate(simulation.cvs):
                f.write("{}\n".format(cv.directive))

        if hasattr(simulation,'biases'):
            f.write("\n# Bias section\n\n")
            for ibias,bias in enumerate(simulation.biases):
                f.write("{}\n".format(cv.biases))

        if  hasattr(simulation,'cvs') and colvar is not None:
            f.write("\n# Print section\n\n")
            f.write("PRINT ...\n"
                    "\tFILE={0}\n"
                    "\tSTRIDE={1}\n"
                    "\tARG=".format(colvar,printstride))
            for icv,cv in enumerate(simulation.cvs):
                f.write("{},".format(cv.name))
            f.write("\n")
            f.write("... PRINT\n")

            
                

    

"""      
class Group():

    def __init__(self, name: str, atoms=[],
                 universe=None, selection=''):
        
        ""Group Class Constructor

        Args:
            name (str): Name of the simulation
            atoms (list, optional): list of atoms in the group. Defaults to [].
            universe (MDAnalysis.Universe, optional): Universe of the system. Defaults to None.
            selection (str, optional): selection for the Universe. Defaults to ''.
        ""
        
        self._name=name
        self._atoms=atoms
        self._universe=universe
        self._selection=selection

        self._natoms=len(atoms)

    @property
    def name(self):
        return self._name

    @property
    def atoms(self):
        return self._atoms

    @property
    def universe(self):
        return self._universe

    @property
    def selection(self):
        return self._selection

    @property
    def natoms(self):
        return self._natoms

    @atoms.setter
    def atoms(self,a):
        self._atoms=a
        self._natoms=len(a)

    @universe.setter
    def universe(self,u):
        self._universe=u
        
    @selection.setter
    def selection(self,s):
        self._selection=s


class _Bias():

    def __init__(self,name: str):
#        ""Bias Class Constructor
#
#        Args:
#           name (str): Name of the bias attribute
#        ""

        self._name=name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,n):
        self._name=n



class UpperWalls(_Bias):

    def __init__(self, name: str, cv: str, kharm: float, at: float, exp=2, eps=1, offset=0.):
        ""UpperWalls bias class constructor

        Args:
           name (str): Label for the bias.
           cv (str): Collective variables to which the bias is applied.
           kharm (float): Force constant for the walls.
           at (float): Position of the wall.
           exp (float, optional): Powers for the walls. Defaults to 2 for each cv.
           eps (float, optional): Normalization for the walls. Defaults to 1 for each cv.
           offset (float, optional): Offset for the start of the wall. Defaults to 0 for each cv.
        "" 

        # Inherit properties of the _Bias Class
        super().__init__(name)

        self._name=name

        if isinstance(cv,list):
            n_walls = len(cv)
        elif isinstance(cv,str):
            n_walls = 1
            cv = [ cv ]

        
"""
        

        
        
