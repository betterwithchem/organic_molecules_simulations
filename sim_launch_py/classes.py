import os
import sim_launch_py.utilities as util
import numpy as np

class Project():
    """
    The project class that stores and manage all the information and methods.

    Attributes:
    - name : name of the project. This defines the directory where the project will be stored and most of path variables. 
    - project_path : absolute path of the project. The name of the directory is derived from the name of the project.
    - systems_path : directory where the data of the simulated systems are stored. This is the parent directory of the systems. 
    - topology_path : directory where the topology files (.top and .itp) of the molecules used in the project are stored and retrieved to be used in each system.
    - init_struct_path : path where initial structures of the molecules used in the project are stored and retrieved to be used in each system.
    - mdp_path : path where default .mdp input files for gromacs are stored.
    - pickle_path : path where the project is saved in .pkl format.
    - job_script_path : path where the template files for submission scripts are stored.
    - systems : list of System() objects in the project.
    - molecules : list of Molecule() objects in the project.
    - gromacs : path to gromacs binary
    - ambertools : path to ambertools binaries directory
 
    Methods:
    - help() : print the help for this class.
    - add_molecule(self, name: str, resname='UNK', structure=None) : add and initialize a Molecule() object to the Project() object.
    - add_system(self,name: str) : add and initialize System() to the Project()
    - new_project(name=None, path=None, overwrite=False) : create new Project()
    - save(self) : save Project() in Project().pickle_path
    - load_project(project_folder: str) : load Project() stored in project_folde/.multisim.pkl. 
    - write_sub_command(self,scriptname: str,hpc_system: str, template: str): write job submission/run scripts for the given HPC system for each System() and a global bash script to call all the job scripts.
    

     #### TODO METHODS ####
    - change_project_path(self,newpath: str) : when implemented, it will allow to copy the project to a new path and
change all the pertinent path attributes.

    """

    def __init__(self,name=None, path=None):
        """Project Class Constructor

        Args:
            name (str, optional): Project name. Defaults to None.
            path (str, optional): Project path. Defaults to None.
        """
        path=os.path.abspath(path)

        if path.rstrip().endswith("/"):
            path = path.rstrip()[:-1]
        
        self._name=name
        self._project_path=os.path.abspath(path)
        self._systems_path="{}/Systems".format(path)
        self._topology_path="{}/Topologies".format(path)
        self._init_struct_path="{}/Initial_structures".format(path)
        self._mdp_path="{}/mdp".format(path)        
        self._pickle_path="{}/.multisim.pkl".format(path)
        self._job_script_path=None
        #self._logfile="{}/project.log".format(path)
        self._systems=list()
        self._molecules=list()        
        self._gromacs=None
        self._ambertools=None

        # look for gromacs and ambertools in the $PATH
        self._checkGromacs()
        self._checkAmberTools()
       
    @property
    def project_path(self):
        return self._project_path

    @property
    def mdp_path(self):
        return self._mdp_path

    @property
    def init_struct_path(self):
        return self._init_struct_path

    @property
    def topology_path(self):
        return self._topology_path

    @property
    def job_script_path(self):
        return self._job_script_path

    @property
    def name(self):
        return self._name

    @property
    def systems(self):
        return self._systems
    
    @property
    def molecules(self):
        return self._molecules

    @property
    def logfile(self):
        return self._logfile

    @property
    def systems_path(self):
        return self._systems_path

    @property
    def gromacs(self):
        return self._gromacs

    @property
    def ambertools(self):
        return self._ambertools
   
    @job_script_path.setter
    def job_script_path(self,sp):
        self._job_script_path=sp

    @logfile.setter
    def logfile(self,log):
        self._logfile=log

    @systems_path.setter
    def systems_path(self,p):
        self._systems_path=p

    @gromacs.setter
    def gromacs(self,gmx):
        self._gromacs=gmx

    @ambertools.setter
    def ambertools(self,amber):
        self._ambertools=amber
        
    @staticmethod
    def help():
        print("""Attributes:
    - name : name of the project. This defines the directory where the project will be stored and most of path variables. 
    - project_path : absolute path of the project. The name of the directory is derived from the name of the project.
    - systems_path : directory where the data of the simulated systems are stored. This is the parent directory of the systems. 
    - topology_path : directory where the topology files (.top and .itp) of the molecules used in the project are stored and retrieved to be used in each system.
    - init_struct_path : path where initial structures of the molecules used in the project are stored and retrieved to be used in each system.
    - mdp_path : path where default .mdp input files for gromacs are stored.
    - pickle_path : path where the project is saved in .pkl format.
    - job_script_path : path where the template files for submission scripts are stored.
    - systems : list of System() in the project.
    - molecules : list of Molecule() in the project.
    - gromacs : path to gromacs binary
    - ambertools : path to ambertools binaries directory
     
    Methods:
    - help() : print the help for this class.
    - add_molecule(self,name=None,resname='UNK', structure=None) : add and initialize Molecule() to the Project()
    - add_system(self,name=None) : add and initialize System() to the Project()
    - new_project(name=None, path=None, overwrite=False) : create new Project()
    - save(self) : save Project() in Project().pickle_path
    - load_project(project_folder: str) : load Project() stored in project_folde/.multisim.pkl. 
    - write_sub_command(self,scriptname: str,hpc_system: str, template: str): write job submission/run scripts for the given HPC system for each System() and a global bash script to call all the job scripts.
    

     #### TODO METHODS ####
    - change_project_path(self, newpath: str) : when implemented, it will allow to copy the project to a new path and change all the pertinent path attributes.
    """)

    def add_molecule(self, name: str, resname='UNK', structure=None):
        """Add molecule to project

        Args:
            name (str): molecule name. Defaults to None.
            resname (str, optional): residue name. Defaults to 'UNK'.
            structure (str, optional): molecular structure file. Defaults to None.
        Returns: 
            newmolecule (object): new molecule.
        """

        # restrict resname to a string of capital letters of length 3
        resname=resname.upper()[0:3]
        
        import shutil
        for mol in self._molecules:
            if name == mol._name:                
                print("Error: a molecule with name {} ({}) already exists in the project!".format(name,resname))
                exit()

        newmolecule=Molecule(name=name,resname=resname,structure=os.path.abspath(structure))

        self._molecules.append(newmolecule)        
        shutil.copy(structure, self._init_struct_path)

        return newmolecule

    def add_system(self, name: str):
        """Add system to project

        Args:
            name (str): system name.
        Returns:
            newsystem (object): new system.
        """

        from sim_launch_py.utilities import create
        
        for sys in self._systems:
            if name==sys.name:
                print("Error: System {} already exists!".format(name))
                exit()

        syspath=self._systems_path+'/'+name
        if os.path.isdir(syspath) is False:
            create(syspath, arg_type='dir')
        
        newsystem=System(name=name,path=syspath,gromacs=self._gromacs)
        self._systems.append(newsystem)

        return newsystem
       

    def new_project(name: str, path: str, overwrite=False):
        """Add new project

        Args:
            name (str): Project Name.
            path (str): Project path.
            overwrite (bool, optional): Overwrite previous project? Defaults to False.

        Returns:
            nproject (object): new project.
        """
        
        from sim_launch_py.utilities import create
        
        path=os.path.abspath(path)
        nproject = Project(path=path, name=name)

        if not os.path.exists(nproject._project_path):
            print("New Project: {}".format(nproject._name))
            print("=" * 100)
            create(nproject._project_path, arg_type='dir', backup=False)
        else:
            if overwrite:
                print("New Project: {}".format(nproject._name))
                print("=" * 100)
                create(nproject._project_path, arg_type='dir', backup=False)
            else:
                print("Error: Folder already exists.\n "
                      "You can change directory name or overwrite with 'overwrite=True' "
                      "(everything inside will be deleted)")
                exit()
        create(nproject._systems_path, arg_type='dir')
        create(nproject._init_struct_path, arg_type='dir')
        create(nproject._topology_path, arg_type='dir')
        create(nproject._mdp_path, arg_type='dir')
        
        return nproject


    def save(self):
        """
        Save project to project folder.
        """
        print("Saving Project...", end="")
        import pickle
        import os
        if os.path.exists(self._pickle_path):
            os.rename(self._pickle_path, self._project_path+ "/.multisim.bck.pkl")
        with open(self._pickle_path, "wb") as file_pickle:
            pickle.dump(self, file_pickle)
        print("done")



    def load_project(project_folder: str): 
        """Load an existing project. The Project object is saved in the project directory every time the command Project.save()
        is used.
        
        Args:
           project_folder (str): location of the project to be loaded.

        Returns:
           project (object): loaded project.
        """

        import pickle
        project_folder = os.path.realpath(project_folder)
        file_pickle = project_folder + "/.multisim.pkl"
        if os.path.exists(file_pickle):
            project = pickle.load(open(file_pickle, "rb"))
            print("Loading Project Name: {}\n".format(project._name))
            if os.path.realpath(project._project_path) != project_folder:
                project.projecy_path = project_folder
            return project
        else:
            print("No project found in '{}'. Use the 'Project.new_project' module to create a new project."
                  "".format(project_folder))

    def write_sub_command(self,scriptname: str,system='myriad', template=None):
        """For each system Write a script to run the simulations and write a global bash script to initiate all systems.

        Args:
           scriptname (str): name of the global bash script.
           system (str, optional): type of script for the systems. Supported are 'bash' and 'myriad' (GridEngine format on Myriad cluster at UCL). Defaults to 'myriad'.
           template (str, optional): custom template to be used for the local script. Defaults to None. 
        """
        import shutil
        recognised_systems=['bash','myriad']

        scriptname=self.project_path+'/'+scriptname
        
        if system=='bash':
            for sys in self.systems:
                print("Sorry... for the moment this doesn't do anything...")                

        elif system=='myriad':

            if template==None:
                template=self.job_script_path+'/myriad.job'
            else:
                template=os.path.abspath(template)

            filename=os.path.basename(template)

            for sys in self.systems:
                shutil.copy(template,sys.path)
                with open(sys.path+'/'+filename,'a') as f:
                    f.write("\n\n")
                    for s in sys.simulations:
                        if s.run:
                            f.write("# {}\n".format(s.name))
                            f.write(s.bash_command+'\n\n')
                sys.run_command='qsub -N {0} {1}\n'.format(self.name+'_'+sys.name, filename)

            with open(scriptname,'w') as f:
                f.write("#!/bin/bash\n\n")
                for sys in self.systems:
                    f.write("cd {}\n".format(sys.path))
                    f.write(sys.run_command)
                    f.write("cd {}\n\n".format(self.project_path))

        else:
            print("ERROR: Type of submission script '{}' not recognized. acceptable values are: ".format(system))
            for rs in recognised_systems:
                print("- {}".format(rs))
            exit()

    def _checkGromacs(self):
        """Look for Gromacs binary and add it to the project.
           This method looks for standard names of gromacs binaries (i.e. 'gmx' and 'gmx_mpi'). If either of these is not found,
           it will prompt the user to input manually the absolute path of the gromacs executable.
           The value of the attribute is the absolute path of the binary (/path/to/gmx or /path/to/gmx_mpi or /path/to/gmx_bin_custom_name).
        """
        import shutil

        if self.gromacs is None:
            if shutil.which('gmx') is not None:
                self.gromacs=shutil.which('gmx')
            elif shutil.which('gmx_mpi') is not None:
                self.gromacs=shutil.which('gmx_mpi')
            else:
                self.gromacs=input("I couldn't find gromacs command, you can add it manually: ")
        elif self.gromacs=='':
            exit("Error: gromacs command is not in the PATH and is not added to the project")

    def _checkAmberTools(self):
        """Look for AmberTools binaries directory and add it to the project.
           This method looks for the antechamber binary. If it is not found, it will prompt the user to input manually the absolute path of the ambertools binaries directory.
           The value of the attribute is the absolute path of the directory (/path/to/bin).
        """
        
        import shutil
        if self.ambertools is None:
            if shutil.which('antechamber') is not None:
                self.ambertools=shutil.which('antechamber')
            else:
                self.ambertools=input("I couldn't find ambertools binaries directory, you can add it manually: ")
        elif self.ambertools=='':
            exit("Error: ambertools binaries is not in the PATH and is not added to the project")

        
class System():
    """
    The system class that stores and manage all the information and methods.

    Attributes:
    - name : name of the system. This is used to create the directory of the system.
    - path : absolute path of the system.
    - molecules : list of Molecule() objects belonging to the system.
    - box : size of the box (of the initial configuration).
    - boxshape : shape of the simulation box.
    - simulations : list of Simulation() objects of the system.
    - run_command : command lines to be used to run the simulations.
    - gromacs : gromacs binary path, inherited from the Project(). 
    - composition : composition of the system. Order of values is the same of the order in self.molecules.
    - atoms : list of Atom() objects in the system.
    - natoms : total numer of atoms in the systems.

    Methods:
    - help() : print the help for this class
    - add_molecule(self, name: str, moltype=None, knownmolecules=None) : add species to the system.
    - add_box(self, box_side: float, shape='cubic') : create simulation box.
    - createSolventBox(self, solvent: object, output_structure="solvent_box.pdb", density=None, nmols=None): add solvent molecules to the simulation box.  
    - insertSolute(self, solute: object, solvent: object, solvent_box="solvent_box.pdb", concentration=0, output_structure="start.pdb"): add solute molecules to the box and remove excess solvent molecules.
    - writeTop(self, atomtypes_path: str, *molecules: objects): write the topology file for the system in gromacs format.
    - add_simulation(self, simtype: str, mdrun_options='', mdp='', print_bash=True, name='',maxwarn=0, start_coord='',gmxbin=''): add simulation to the system.
    - print_command(self, bash_file): print the bash script file to run the simulations.
   
    """

    def __init__(self,name,path=None,gromacs=None):
        """System Class Constructor

        Args:
            name (str): Systen name
            path (str, optional): System path. Defaults to None.
            gromacs (str, optional): Gromacs binary name. Defaults to None.
        """

        self._name=name
        self._path=path
        self._molecules=list()
        #self._temperature=0
        self._box=list()
        self._boxshape=None
        self._simulations=list()
        self._run_command=None
        self._gromacs=gromacs
        self.composition=list()
        self._atoms=list()
        self._natoms=0        

    @property
    def name(self):
        return self._name

    @property
    def path(self):
        return self._path

    @property
    def molecules(self):
        return self._molecules

    #@property
    #def temperature(self):
    #    return self._temperature

    @property
    def box(self):
        return self._box

    @property
    def boxshape(self):
        return self._boxshape

    @property
    def simulations(self):
        return self._simulations

    @property
    def run_command(self):
        return self._run_command

    @property
    def gromacs(self):
        return self._gromacs

    @property
    def atoms(self):
        return self._atoms

    @property
    def natoms(self):
        return self._natoms

    @gromacs.setter
    def gromacs(self,gmx):
        self._gromacs=gmx
    
    #@temperature.setter
    #def temperature(self,T):
    #    self._temperature=T

    @natoms.setter
    def natoms(self,n):
        self._natoms=n

    @box.setter
    def box(self,box):
        self._box=box

    @boxshape.setter
    def boxshape(self,shape):
        self._boxshape=shape

    @run_command.setter
    def run_command(self,command):
        self._run_command=command
    
    def add_molecule(self, name: str, moltype=None, knownmolecules=None):
        """Add new molecule type to the system.

        Args:
           name (str): name of the molecule.
           moltype (str, optional): generic attribute of the molecule. Defaults to None.
           knownmolecules (list of Molecule(), optional): molecules already present in the project.
        Return: 
           newmolecule (object): new molecule
        """

        import copy
        
        for mol in knownmolecules:
            if name==mol.name:
                newmolecule=copy.deepcopy(mol)
                newmolecule._mol_attributes=moltype
                self._molecules.append(newmolecule)
                return newmolecule

        print("Error: Couldn't add molecule {} to system {}. Molecule unknown.".format(name,self.name))
        exit()
        

    def addBox(self,box_side: float, shape='cubic'):
        """ Create simulation box.
 
        Args: 
           box_side (float): side of the box.
           shape (str, optional): shape of the box. Defaults to 'cubic'.
        """

        self.boxshape=shape
        
        if shape=='cubic':
            self.box=[box_side]*3+[90.0]*3
        elif shape=='dodecahedron':
            self.box=[box_side, box_side, 0.5*np.sqrt(2)*box_side] + [60.0]*3 # see Gromacs manual
        elif shape=='octahedron':
            self.box=[box_side, 2/3*np.sqrt(2)*box_side, 1/3*np.sqrt(6)*box_side]+[71.53, 109.47, 71.53] # see Gromacs manual
           

    def createSolventBox(self, solvent: object, output_structure="solvent_box.pdb", density=None, nmols=None):
        """ Add solvent molecules to the simulation box.
 
        Args: 
           solvent (Molecule() object) : molecule species to be used as solvent.
           output_structure (str, optional) : name of the structure file with the solvent box. Defaults to 'solvent_box.pdb'.
           density (float, optional) : desired density of solvent molecules in [g/L]. Defaults to None.
           nmols (int, optional) : desired number of solvent molecules. Defaults to None. If both density and number of molecules are provided, the size of the box is changed.
        """
        if (density is not None) and (nmols is not None):
            
            volume_box=solvent.mw*10/(density*6.022)*nmol
            self.box[0]=volume_box**(1/3)
            print("Both density and nmol have been defined: changing side of the box to {}.".format(self.box[0]))

        else:
            if self.boxshape=='cubic':
                volume_box=self.box[0]**3
            elif self.boxshape=='dodecahedron':
                volume_box=0.5*np.sqrt(2)*self.box[0]**3
            elif self.boxshape=='octahedron':
                volume_box=4/9*np.sqrt(3)*self.box[0]**3

        nmols=util.estimate_n_molecules(volume_box,solvent.mw,density)-1
        
        os.system("{0} -nobackup editconf -f {1} -o {2} -box {3} -bt {4}  -c".format(self.gromacs,solvent.structure_path,
                                                                                     output_structure,
                                                                                     self.box[0],
                                                                                     self.boxshape))
               
        os.system("{0} -nobackup insert-molecules -f {1} -o {1} -ci {1} -nmol {2} -try 20000".format(self.gromacs,output_structure,
                                                                        nmols))

        inserted_mols=util.check_number_molecules(output_structure,solvent)

        if inserted_mols!=(nmols+1):
            print("Error! Inserted number of mol{} molecules required, but after {} trials only {} where inserted!".format(nmols, 20000, inserted_mols))
            exit()

        solvent.nmols=nmols

        self._updateComposition()


        
    def insertSolute(self, solute: object ,solvent: object, solvent_box="solvent_box.pdb", concentration=0.,
                     output_structure="start.pdb"):
        """ Insert solute molecules in the simulation box.

        Args: 
           solute (Molecule() object) : molecule type to be used as solute.
           solvent (Molecule() object) : molecule type to be used as solvent. 
           solvent_box (str, optional) : structure file where to add solute molecules. Defaults to 'solvent_box.pdb'.
           concentration (float, optional) : desired concentration of solute molecules in [g/L]. A density of 0 is a flag for simulations of dilute systems with only 1 solute molecule.
           output_structure (str, optional) : structure file where to save the final configuration. Defaults to 'start.pdb'.
        """

        from sim_launch_py.utilities import check_number_molecules as cnm
        
        if concentration == 0:
            # then single molecule
            nmols=1
        elif concentration > 0:
            from math import ceil
            nmols=util.estimate_n_molecules( self.box[0]**3,solute.mw,concentration )

        os.system("{0} insert-molecules -f {1} -o {2} -ci {3} -nmol {4} -try 10000 -replace {5} ".format( self.gromacs,
                                                                                                          solvent_box,
                                                                                                          output_structure,
                                                                                                          solute.structure_path,
                                                                                                          nmols,
                                                                                                          solvent.resname))

        inserted_mols=util.check_number_molecules(output_structure,solute)

        if inserted_mols!=nmols:
            print("Error! Inserted number of mol{} molecules required, but after {} trials only {} where inserted!".format(nmols, 10000, inserted_mols))
            exit()

        solvent.nmols=cnm(output_structure,solvent)
        solute.nmols=nmols
        self._updateComposition()

    def writeTop(self,atomtypes_path: str): #, *molecules: object):
        """ Write the topology file in Gromacs format)

        Args: 
            atomtypes_path (str) : path of the file with the atomtypes definition in Gromacs format.
        """

        with open(self.path+'/topol.top','w') as top:

            top.write(";\n\
; Topology for system {}\n\
;\n\n".format(self.name))

            top.write("[ defaults ]\n\
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n\
1               2               yes             0.5          0.83333333  \n\
\n\
[ atomtypes ]\n\
; name    at.num    mass    charge ptype  sigma      epsilon\n")

            with open(atomtypes_path,'r') as af:
                for line in af:
                    top.write(line)

            top.write("\n\n")

            for mol in molecules:
                with open(mol.include_path,'r') as itp:
                    top.write(";\n; {}\n;\n".format(mol.name))
                    for line in itp:
                        top.write(line)
                top.write("\n\n")

            top.write("[ system ]\n; Name\nA bunch of molecules floating around without rhyme or reason\n\n")

            top.write("[ molecules ]\n")
            
            for mol in self.molecules:
                top.write("{0}\t{1}\n".format(mol.resname,mol.nmols))

            

    def add_simulation(self, name: str, simtype: str, mdrun_options='', mdp='', print_bash=True,
                       maxwarn=0, start_coord='', plumed=None):
        """ Add a simulation to the system.
         
        Args: 
           name (str) : name of the simulation.
           simtype (str) : type of simulation to be run.
           mdrun_options (str, optional) : optional flags to be added to the mdrun command. Defaults to ''.
           mdp (str, optional) : name of the custom mdp file to use. Defaults to ''.
           print_bash (bool, optional) : generate a bash command to launch the simulation. Defaults to True.
           maxwarn (int, optional) : maximum number of warnings for grompp. Defaults to 0.
           start_coord (str, optional) : starting configuration of the simulation. Defaults to ''.
           plumed (str, optional) : plumed file to be used in the simulation. Defaults to None.

        Returns:
          sim : Simulation() object.
        """
       
        import sim_launch_py.gromacs as gmx

        accepted_types=['em','md']
        
        if start_coord=='':
            if len(self.simulations)>0:
                prev_sim=self.simulations[-1].name
            else:
                print("ERROR: 	For the first simulation of the system you need to specify the name of the starting configuration!")
                exit()
            start_coord=self.path+'/'+prev_sim+'.gro'
        else:
            start_coord=os.path.abspath(start_coord)
                
        if simtype=='md':
            if len(self.simulations)>0:
                prev_sim=self.simulations[-1].name
            elif len(self.simulations)==0:
                prev_sim='start'
            sim=gmx.MD(name,
                       mdrun_options=mdrun_options, coord=start_coord, topology=self.path+'/topol.top',
                       path_mdp=mdp, maxwarn=maxwarn,
                       path_input=self.path,path_output=self.path,print_bash=True,gmxbin=self.gromacs,plumed=plumed)
            
        elif simtype=='em':
            sim=gmx.EnergyMinimization(name,
                                       mdrun_options=mdrun_options, coord=start_coord, topology=self.path+'/topol.top',
                                       path_mdp=mdp, maxwarn=maxwarn,
                                       path_input=self.path,path_output=self.path,print_bash=True, gmxbin=self.gromacs)
        #elif simtype=='wtmetad':
        #    sim=gmx.WTMD(name='wtmd')
        else:
            print("Error: simulation type not recognized, accepted types are {}, your input was {}.".format(accepted_types,simtype))
            exit()

        self.simulations.append(sim)

        return sim

    def setSimsToRun(self, sims_to_run: list):
        """ Define which simulations in a system are to be run.

        Args:
           sims_to_run (Simulation objects): list of simulations to run.
        """
        for sim in self.simulations:
            sim.run=False

        for sim in sims_to_run:
            sim.run=True
            
    
    def print_command(self, bash_file: str):
        """ Print bash script to run the simulations of the system.

        Args:
           bash_file (str): name of the bash file.
        """

        bash_file=self.path+'/'+bash_file

        with open(bash_file,'w') as f:
            f.write("#!/bin/bash\n\n")
            for s in self.simulations:
                f.write("# {}\n".format(s.name))
                f.write(s.bash_command+'\n\n')


    def _updateComposition(self):
        """ Update the composition attribute, the atoms and the number of atoms in the system.
        """

        import copy
        
        self.composition=[0 for mol in self.molecules]
        self._atoms=[]

        for imol,mol in enumerate(self.molecules):
            self.composition[imol]=mol.nmols
            for i in range(mol.nmols):
                for iatom,atom in enumerate(mol.atoms):
                    new_atom=copy.deepcopy(mol.atoms[iatom])
                    self._atoms.append(new_atom)

        self._natoms=len(self._atoms)
        self._updateAtomID()
        

    def _updateAtomID(self):
        """ Update the atom ID of the atoms composing the system.
        """

        iatom=0
        molatom=0
        
        for imol,mol in enumerate(self.molecules):
            for n in range(self.composition[imol]):
                for i in range(mol._natoms):
                    self._atoms[iatom]._atomID=iatom+1
                    mol._atoms[molatom]._atomID=iatom+1
                    iatom+=1
                    molatom+=1
                molatom=0
                
            

class Molecule():
    """
    The molecule class that stores and manage all the information and methods.

    Attributes:


    Methods:


    """


    def __init__(self,name=None, resname='UNK', structure=None):

        self._name=name
        self._resname=resname
        self._structure_path=structure
        self._topology_path=None
        self._include_path=None
        self._mw=0
        self._mol_attributes=[]
        self._nmols=0
        self._natoms=0
        self._atoms=list()

        if structure is not None:

            import MDAnalysis as mda
            import warnings
            warnings.filterwarnings('ignore')

            u=mda.Universe(structure)
            self._natoms=u.atoms.n_atoms
            for iatom in range(self._natoms):

                new_atom=Atom(u.atoms.names[iatom],atomID=u.atoms.ids[iatom],resname=resname,
                         mass=u.atoms.masses[iatom],element=u.atoms.types[iatom])
                                
                self._mw+=u.atoms.masses[iatom]
                self._atoms.append(new_atom)


    @property
    def name(self):
        return self._name

    @property
    def resname(self):
        return self._resname

    @property
    def structure_path(self):
        return self._structure_path

    @property
    def topology_path(self):
        return self._topology_path

    @property
    def include_path(self):
        return self._include_path

    @property
    def mw(self):
        return self._mw

    @property
    def mol_attributes(self):
        return self._mol_attributes

    @property
    def natoms(self):
        return self._natoms

    @property
    def nmols(self):
        return self._nmols

    @property
    def atoms(self):
        return self._atoms
    
    @resname.setter
    def resname(self,resname):
        self._resname=resname

    @structure_path.setter
    def structure_path(self,path):
        self._structure_path=path

    @topology_path.setter
    def topology_path(self,path):
        self._topology_path=path

    @include_path.setter
    def include_path(self,path):
        self._include_path=path

    @mw.setter
    def mw(self,mw):
        self._mw=mw

    @mol_attributes.setter
    def mol_attributes(self,attr):
        self._mol_attributes.append(attr)

    @natoms.setter
    def natoms(self,n):
        self._natoms=n

    @nmols.setter
    def nmols(self,n):
        self._nmols=n

    @staticmethod
    def help():
        print("""Class for Molecule objects
Each Molecule() has the following attributes:
-name
-resname 
-structure_path : name and position of the structure file 
-topology_path : name and position of the topology (.top) file with all force field parameters, including atom types definition
-include_path : name and position of the include topology (.itp) file with the definition of the molecule
-mw : molecular weight in g/mol""")

 
class Atom():

    def __init__(self,name: str, mass=0, atomtype=None, atomID=0, resnum=0, resname=None, element=None):
        """Atom Class Constructor

        Args:
            name (str): Atom name. Defaults to None.
            mass (float, optional): Mass of the atom. Defaults to 0.
            atomtype (str, optional): AtomType of the atom. Defaults to None.
            atomID (int, optional): AtomID of the atom. Defaults to 0.
            resnum (int, optional): Residue Number of the atom. Defaults to 0.
            resname (str, optional): Residue Name of the atom. Defaults to None.
            element (str, optional): Element of the atom. Defaults to None.
        """

        self._name=name
        self._mass=mass
        self._atomtype=atomtype
        self._atomID=atomID
        self._resnum=resnum
        self._resname=resname
        self._element=element

    @property
    def name(self):
        return self._name

    @property
    def mass(self):
        return self._mass

    @property
    def atomtype(self):
        return self._atomtype

    @property
    def atomID(self):
        return self._atomID

    @property
    def resnum(self):
        return self._resnum

    @property
    def resname(self):
        return self._resname

    @property
    def element(self):
        return self._element

    @name.setter
    def name(self,n):
        self._name=n

    @mass.setter
    def mass(self,mw):
        self._mass=mw

    @atomtype.setter
    def atomtype(self,atype):
        self._atomtype=atype

    @atomID.setter
    def atomID(self,aid):
        self._atomID=aid

    @resnum.setter
    def resnum(self,rn):
        self._resnum=rn

    @resname.setter
    def resname(self,rn):
        self._resname=rn

    @element.setter
    def element(self,el):
        self._element=el

