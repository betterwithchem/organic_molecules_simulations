
class _CV():

    def __init__(self, name: str, print_cv=True):
        """CV Class Constructur
   
        Args:
           name (str): Name of the collective variable
           print_cv (bool, optional): Print the value of the variable? Defaults to True.
        """
        
        self._name=name
        self._print_cv=print_cv
        self._cvtype=None

    @property
    def name(self):
        return self._name

    @property
    def print_cv(self):
        return self._print_cv

    @property
    def cvtype(self):
        return self._cvtype

    @name.setter
    def name(self,n):
        self._name=n

    @print_cv.setter
    def print_cv(self,p):
        self._print_cv=p

    @cvtype.setter
    def cvtype(self,cv):
        self._cvtype=cv

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
        self._cvtype='TORSION'
        
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


class PotentialEnergy(_CV):

    def __init__(self, name: str):
        """Potential Energy collective variable class constructor

        Args:
           name (str) : Name of the collective variable.
        """

        # Inherit properties of the _CV Class
        super().__init__(name)

        self._name=name
        self._cvtype='ENERGY'

        self._directive="{}: ENERGY".format(name)

    @property
    def directive(self):
        return self._directive
    

class _Bias():

    def __init__(self,name: str):
        """Bias Class Constructor

        Args:
           name (str): Name of the bias attribute
        """

        self._name=name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,n):
        self._name=n

        
class Metad(_Bias):

    def __init__(self, name: str, cv: str, sigma=None, height=None, temp=300, pace=500,
                 hills_file='HILLS', biasfactor=None,
                 grid_min=None, grid_max=None, grid_spacing=None):
        """Metadynamics bias class constructor
     
        Args:
           name (str): label for the metadynamics directive
           cv (str or list of str): collective variables to bias
           hills_file (str, optional): 
        """

        # Inherit properties of the _Bias Class
        super().__init__(name)

        #self._name=name

        if isinstance(cv,str):
            ncvs=1
            cv=[cv]
        elif isinstance(cv,list):
            ncvs=len(cv)

        if ncvs==0:
            print("Error: you need at least one CV to apply bias!")
            exit()
            
        if isinstance(sigma,float) or isinstance(sigma,int):
            if ncvs==1:
                sigma=[sigma]

        if len(cv)!=len(sigma):
            print("Error: The number of sigma values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} sigma values ({})".format(len(cv),cv,len(sigma),sigma))
            exit()

        do_grid=False
        if any([grid_min, grid_max, grid_spacing]):
            if all([grid_min, grid_max, grid_spacing]):
                do_grid=True
            else:
                print("Error: To use grids during a Metadynamics simulation all grid-related variables need to be provided (grid_min, grid_max, grid_spacing).")
                exit()
        
        
        self.directive="METAD ...\n"\
            "    ARG="
        for var in cv:
            self.directive+="{},".format(var)
        self.directive+="\n"\
            "    FILE={}\n".format(hills_file)
        self.directive+="    SIGMA="
        for s in sigma:
            self.directive+="{},".format(s)
        self.directive+="\n"\
            "    HEIGHT={}\n"\
            "    PACE={}\n".format(height,pace)

        if biasfactor is not None:
            self.directive+="    BIASFACTOR={}\n"\
                "    TEMP={}\n".format(biasfactor,temp)
        if do_grid:
            self.directive+="    GRID_MIN={}\n"\
            "    GRID_MAX={}\n"\
            "    GRID_SPACING={}\n".format(grid_min,grid_max,grid_spacing)
    
        self.directive+="... METAD\n"


class UpperWalls(_Bias):

    def __init__(self, name: str, cv: str, kappa=None, at=None, exp=None, eps=None, offset=None):
        """Upper Walls bias class constructor

        Args:
           name (str) : label of the Upper Walls directive.
           cv (str or list of str): collective variables to bias.          
           at (float or list of float) : position of the wall.
           kappa (float or list of float) : energy constant of the wall. 
           exp (float or list of float, optional) : exponent determining the power law. Defaults to 2 for all walls.
           eps (float or list of float ,optional) : rescaling factor. Defaults to 1 for all walls.
           offset (float or list of float, optional) : offset for the definition of the wall. Defaults to 0 for all walls.
        """

        # Inherit properties of the _Bias Class
        super().__init__(name)

        #self._name=name

        if isinstance(cv,str):
            ncvs=1
            cv=[cv]
        elif isinstance(cv,list):
            ncvs=len(cv)

        if ncvs==0:
            print("Error: you need at least one CV to apply bias!")
            exit()

        if isinstance(at,float) or isinstance(at,int):
            if ncvs==1:
                at=[at]

        if len(cv)!=len(at):
             print("Error: The number of wall position values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} at values ({})".format(len(cv),cv,len(at),at))
             exit()

        if isinstance(kappa,float) or isinstance(kappa,int):
            if ncvs==1:
                kappa=[kappa]

        if len(cv)!=len(kappa):
             print("Error: The number of kappa values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} kappa values ({})".format(len(cv),cv,len(kappa),kappa))
             exit()

        if exp==None:
            exp=[2 for c in cv]
        else:
            if isinstance(exp,float) or isinstance(exp,int):
                if ncvs==1:
                    exp=[exp]

            if len(cv)!=len(exp):
                print("Error: The number of exp values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} exp values ({})".format(len(cv),cv,len(exp),exp))
                exit()

        if eps==None:
            eps=[1 for c in cv]
        else:
            if isinstance(eps,float) or isinstance(eps,int):
                if ncvs==1:
                    eps=[eps]

            if len(cv)!=len(eps):
                print("Error: The number of epsilon values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} epsilon values ({})".format(len(cv),cv,len(eps),eps))
                exit()

        if offset==None:
            offset=[0 for c in cv]
        else:
            if isinstance(offset,float) or isinstance(offset,int):
                if ncvs==1:
                    offset=[offset]

            if len(cv)!=len(offset):
                print("Error: The number of offset values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} offset values ({})".format(len(cv),cv,len(offset),offset))
                exit()

        self.directive="UPPER_WALLS ...\n"
        self.directive+="    LABEL={}\n"\
            "    ARG=".format(name)
        for c in cv:
            self.directive+="{},".format(c)
        self.directive+="\n"
        self.directive+="    AT="
        for a in at:
            self.directive+="{},".format(a)
        self.directive+="\n"
        self.directive+="    KAPPA="
        for k in kappa:
            self.directive+="{},".format(k)
        self.directive+="\n"
        self.directive+="    EXP="
        for e in exp:
            self.directive+="{},".format(e)
        self.directive+="\n"
        self.directive+="    EPS="
        for e in eps:
            self.directive+="{},".format(e)
        self.directive+="\n"
        self.directive+="    OFFSET="
        for o in offset:
            self.directive+="{},".format(o)
        self.directive+="\n... UPPER_WALLS\n\n"


class LowerWalls(_Bias):

    def __init__(self, name: str, cv: str, at=None, kappa=None,
                 exp=None, eps=None, offset=None):
        """Lower Walls bias class constructor

        Args:
           name (str) : label of the Lower Walls directive.
           cv (str or list of str): collective variables to bias.          
           at (float or list of float) : position of the wall.
           kappa (float or list of float) : energy constant of the wall. 
           exp (float or list of float, optional) : exponent determining the power law. Defaults to 2 for all walls.
           eps (float or list of float ,optional) : rescaling factor. Defaults to 1 for all walls.
           offset (float or list of float, optional) : offset for the definition of the wall. Defaults to 0 for all walls.
        """

        # Inherit properties of the _Bias Class
        super().__init__(name)

        #self._name=name

        if isinstance(cv,str):
            ncvs=1
            cv=[cv]
        elif isinstance(cv,list):
            ncvs=len(cv)

        if ncvs==0:
            print("Error: you need at least one CV to apply bias!")
            exit()

        if isinstance(at,float) or isinstance(at,int):
            if ncvs==1:
                at=[at]

        if len(cv)!=len(at):
             print("Error: The number of wall position values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} at values ({})".format(len(cv),cv,len(at),at))
             exit()

        if isinstance(kappa,float) or isinstance(kappa,int):
            if ncvs==1:
                kappa=[kappa]

        if len(cv)!=len(kappa):
             print("Error: The number of kappa values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} kappa values ({})".format(len(cv),cv,len(kappa),kappa))
             exit()

        if exp==None:
            exp=[2 for c in cv]
        else:
            if isinstance(exp,float) or isinstance(exp,int):
                if ncvs==1:
                    exp=[exp]

            if len(cv)!=len(exp):
                print("Error: The number of exp values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} exp values ({})".format(len(cv),cv,len(exp),exp))
                exit()

        if eps==None:
            eps=[1 for c in cv]
        else:
            if isinstance(eps,float) or isinstance(eps,int):
                if ncvs==1:
                    eps=[eps]

            if len(cv)!=len(eps):
                print("Error: The number of epsilon values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} epsilon values ({})".format(len(cv),cv,len(eps),eps))
                exit()

        if offset==None:
            offset=[0 for c in cv]
        else:
            if isinstance(offset,float) or isinstance(offset,int):
                if ncvs==1:
                    offset=[offset]

            if len(cv)!=len(offset):
                print("Error: The number of offset values need to be equal to the number of collective variables."\
                  "Given {} cv values ({}) and {} offset values ({})".format(len(cv),cv,len(offset),offset))
                exit()

        self.directive="LOWER_WALLS ...\n"
        self.directive+="    LABEL={}\n"\
            "    ARG=".format(name)
        for c in cv:
            self.directive+="{},".format(c)
        self.directive+="\n"
        self.directive+="    AT="
        for a in at:
            self.directive+="{},".format(a)
        self.directive+="\n"
        self.directive+="    KAPPA="
        for k in kappa:
            self.directive+="{},".format(k)
        self.directive+="\n"
        self.directive+="    EXP="
        for e in exp:
            self.directive+="{},".format(e)
        self.directive+="\n"
        self.directive+="    EPS="
        for e in eps:
            self.directive+="{},".format(e)
        self.directive+="\n"
        self.directive+="    OFFSET="
        for o in offset:
            self.directive+="{},".format(o)
        self.directive+="\n... LOWER_WALLS\n\n"

        

            

    
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
                f.write("{}\n".format(bias.directive))

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

            
