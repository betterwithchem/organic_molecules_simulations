import os
import sim_launch_py.utilities as util

class Project():
    """
    The project class that stores and manage all the information and methods.

    Attributes:
    - name : name of the project. This defines the directory where the project will be stored and most of path variab
les. 
    - project_path : absolute path of the project. The name of the directory is derived from the name of the project.
    - systems_path : directory where the data of the simulated systems are stored. This is the parent directory of th
e systems. 
    - topology_path : directory where the topology files (.top and .itp) of the molecules used in the project are sto
red 
                      and retrieved to be used in each system.
    - init_struct_path : path where initial structures of the molecules used in the project are stored and retrieved
                         to be used in each system.
    - mdp_path : path where default .mdp input files for gromacs are stored.
    - pickle_path : path where the project is saved in .pkl format.
    - job_script_path : path where the template files for submission scripts are stored.
    - systems : list of System() in the project.
    - molecules : list of Molecule() in the project.
    - version : version of the code (given in form of a YYYYMMDD date)
    - gromacs : path to gromacs binary
    - ambertools : path to ambertools binaries directory

    #### TODO ATTRIBUTES ####
    ** these path could be added in an automatic fashion by using shutil.which() on "gmx","gmx_mpi","tleap",... and f
ail in case they are not found (or give the possibility to add them manually at prompt). I would personally go for th
e automated fashion (with clear indication in a log file), because, in any case, they need to be in the environment o
f the system with all the loaded libraries (so, sourcing need to be done BEFORE than running the project). Also, in t
he cases where these programs are NOT needed (e.g. when adding molecules to the project/systems and systems to the pr
oject, it is not important that they are defined at runtime).
 
    Methods:
    - help() : print the help for this class.
    - add_molecule(self,name=None,resname='UNK', structure=None) : add and initialize Molecule() to the Project()
    - add_system(self,name=None) : add and initialize System() to the Project()
    - new_project(name=None, path=None, overwrite=False) : create new Project()
    - save(self) : save Project() in Project().pickle_path
    - load_project(project_folder: str) : load Project() stored in project_folde/.multisim.pkl. 
    - write_sub_command(self,scriptname: str,hpc_system: str, template: str): write job submission/run scripts for th
e given HPC system for each System() and a global bash script to call all the job scripts.
    - addGromacs

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
        self._mdp_path="{}/mdp".format(path)
        self._init_struct_path="{}/Initial_structures".format(path)
        self._topology_path="{}/Topologies".format(path)
        self._systems_path="{}/Systems".format(path)
        self._pickle_path="{}/.multisim.pkl".format(path)
        self._logfile="{}/project.log".format(path)
        self._name=name
        self._systems=list()
        self._molecules=list()
        self._job_script_path=None
        self._gromacs=None
        self._ambertools=None

        self.version = '20230804'

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
    - topology_path : directory where the topology files (.top and .itp) of the molecules used in the project are stored 
                      and retrieved to be used in each system.
    - init_struct_path : path where initial structures of the molecules used in the project are stored and retrieved
                         to be used in each system.
    - mdp_path : path where default .mdp input files for gromacs are stored.
    - pickle_path : path where the project is saved in .pkl format.
    - job_script_path : path where the template files for submission scripts are stored.
    - systems : list of System() in the project.
    - molecules : list of Molecule() in the project.
    - version : version of the code (given in form of a YYYYMMDD date)

    #### TODO ATTRIBUTES ####
 
    Methods:
    - help() : print the help for this class.
    - add_molecule(self,name=None,resname='UNK', structure=None) : add and initialize Molecule() to the Project()
    - add_system(self,name=None) : add and initialize System() to the Project()
    - new_project(name=None, path=None, overwrite=False) : create new Project()
    - save(self) : save Project() in Project().pickle_path
    - load_project(project_folder: str) : load Project() stored in project_folde/.multisim.pkl. 
    - write_sub_command(self,scriptname: str,hpc_system: str, template: str): write job submission/run scripts for the given HPC system for each System() and a global bash script to call all the job scripts.

     #### TODO METHODS ####
    - _check_needed_binaries(self) : this method, when implemented, looks for gromacs, ambertools (and plumed? and lammps?) binaries and add them to the project. If they are not found, a prompt could ask the user to give them. 
    - change_project_path(self,newpath: str) : when implemented, it will allow to copy the project to a new path and change all the pertinent path attributes.
    

    """)

    def add_molecule(self,name,resname='UNK', structure=None):
        """Add molecule to project

        Args:
            name (str): molecule name. Defaults to None.
            resname (str, optional): residue name. Defaults to 'UNK'.
            structure (str, optional): molecular structure file. Defaults to None.
        """

        # restrict resname to a string of capital letters of length 3
        resname=resname.upper()[0:3]
        
        import shutil
        for mol in self._molecules:
            if name == mol._name:                
                print("Ignoring molecule {} ({}), it already exists a molecule with the same name in the project.".format(name,resname))
                return

        newmolecule=Molecule(name=name,resname=resname,structure=os.path.abspath(structure))

        self._molecules.append(newmolecule)        
        shutil.copy(structure, self._init_struct_path)


    def add_system(self,name=None):
        """Add system to project

        Args:
            name (str, optional): system name. Defaults to None.
        Returns:
            newsystem (class): new system
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
        


    def new_project(name=None, path=None, overwrite=False):
        """Add new project

        Args:
            name (str, optional): Project Name. Defaults to None.
            path (str, optional): Project path. Defaults to None.
            overwrite (bool, optional): Overwrite previous project? Defaults to False.

        Returns:
            nproject (class): new project
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
        :return:
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
        """
        Load an existing project. The Project object is saved in the project directory every time the command Project.save()
        is used.
        :param project_folder: project folder specified in the new_project function
        :param use_backup: Use the Project object for the previous save
        :return: Project object
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

    def write_sub_command(self,scriptname,system='', template=None):
        """
        For each system Write a script to run the simulations and write a global bash script to initiate all systems.
        :param scriptname: name of the global bash script.
        :param system: type of script for the systems. Supported are 'bash' and 'myriad' (GridEngine format on Myriad cluster at UCL).
        :param template: custom template to be used for the local script. 
        """
        import shutil
        recognised_systems=['bash','myriad']

        scriptname=self.project_path+'/'+scriptname

        
        if system=='bash':
            for sys in self.systems:
                print("hey!")

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
        """                                                                                                          
        Look for Gromacs binary and add it to the project                                                            
        The value of the attribute is the actual path of the binary (/path/to/gmx or /path/to/gmx_mpi or /path/to/gm\
x_bin_custom_name)                                                                                                   
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
        """                                                                                                          
        Look for AmberTools binaries path and add it to the project                                                  
        The value of the attribute is the path to the directory where the binaries are (/path/to/amber/bin)          
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
    - molecules : list of Molecule() belonging to the system.
    - box : size of the box (of the initial configuration).
    - simulations : list of Simulation() performed for the system.
    - run_command : command lines used to run the simulations.

    Methods:
    - help() : print the help for this class
    - add_molecule(self,name,knownmolecules=None,mol_attr: list())
    - createSolventBox(self,solvent,output_structure="solvent_box.pdb",density=None,nmols=None):  
    - insertSolute(self,solute,solvent,solvent_box="solvent_box.pdb",concentration=0, output_structure="start.pdb"):
    - writeTop(self,atomtypes_path,*molecules):
    - add_simulation(self, simtype: str, mdrun_options='', mdp='', print_bash=True, name='',maxwarn=0,start_coord='',
gmxbin=''):
    - print_command(self, bash_file):
   
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
        self._temperature=0
        self._box=list()
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

    @property
    def temperature(self):
        return self._temperature

    @property
    def box(self):
        return self._box

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
    
    @temperature.setter
    def temperature(self,T):
        self._temperature=T

    @natoms.setter
    def natoms(self,n):
        self._natoms=n

    @box.setter
    def box(self,box):
        if isinstance(box,float) or isinstance(box,int) or len(box)==1:
            self._box=[box, box, box, 90, 90, 90]
        elif len(box)==6:
            self._box=box

    @run_command.setter
    def run_command(self,command):
        self._run_command=command

    
    def add_molecule(self,name,moltype=None,knownmolecules=None):

        """Add new molecule type to the system.

        Args:
           name (str): name of the molecule.
           moltype (str, optional): generic attribute of the molecule. Defaults to None.
           knownmolecules (list of Molecule(), optional): molecules already present in the project.
        Return: 
           newmolecule (class): new molecule
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
        


    def createSolventBox(self,solvent,output_structure="solvent_box.pdb",density=None,nmols=None):
        
        if (density is not None) and (nmols is not None):
            
            volume_box=solvent.mw*10/(density*6.022)*nmol
            self.box[0]=volume_box**(1/3)
            print("Both density and nmol have been defined: changing side of the box to {}.".format(self.box[0]))

        else:
            volume_box=self.box[0]**3

        nmols=util.estimate_n_molecules(volume_box,solvent.mw,density)-1
        
        os.system("{0} -nobackup editconf -f {1} -o {2} -box {3} {3} {3} -angles 90 90 90 -c".format(self.gromacs,solvent.structure_path,
                                                                                           output_structure,
                                                                                                     self.box[0]*(1.05)))
               
        os.system("{0} -nobackup insert-molecules -f {1} -o {1} -ci {1} -nmol {2} -try 20000".format(self.gromacs,output_structure,
                                                                        nmols))

        print("Checking number of molecules..")
        inserted_mols=util.check_number_molecules(output_structure,solvent)

        if inserted_mols!=(nmols+1):
            print("Error! Inserted number of mol{} molecules required, but after {} trials only {} where inserted!".format(nmols, 20000, inserted_mols))
            exit()

        solvent.nmols=nmols

        self._updateComposition()


        
    def insertSolute(self,solute,solvent,solvent_box="solvent_box.pdb",concentration=0, output_structure="start.pdb"):        

        from sim_launch_py.utilities import check_number_molecules as cnm
        
        if concentration == 0:
            # then single molecule
            nmols=1
        elif concentration > 0:
            from math import ceil
            nmols=util.estimate_n_molecules( self.box[0]**3,solute.mw,concentration )

        #print(self.name,solute.name,solvent.name,concentration,nmols)

        os.system("{0} insert-molecules -f {1} -o {2} -ci {3} -nmol {4} -try 10000 -replace {5} ".format( self.gromacs,
                                                                                                          solvent_box,
                                                                                                          output_structure,
                                                                                                          solute.structure_path,
                                                                                                          nmols,
                                                                                                          solvent.resname))

        print("Checking number of molecules..")
        inserted_mols=util.check_number_molecules(output_structure,solute)

        if inserted_mols!=nmols:
            print("Error! Inserted number of mol{} molecules required, but after {} trials only {} where inserted!".format(nmols, 10000, inserted_mols))
            exit()

        solvent.nmols=cnm(output_structure,solvent)
        solute.nmols=nmols
        self._updateComposition()

    def writeTop(self,atomtypes_path,*molecules):

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

            top.write("[ system ]\n; Name\nA bunch of molecules floating around without rhyme or reason\n")

            top.write("[ molecules ]\n")
            
            for mol in molecules:
                top.write("{0}\t{1}\n".format(mol.resname,mol.nmols))

            

    def add_simulation(self, simtype: str, mdrun_options='', mdp='', print_bash=True, name='',maxwarn=0,start_coord='',gmxbin='',plumed=None):
       
        import sim_launch_py.gromacs as gmx

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
                       path_input=self.path,path_output=self.path,print_bash=True,gmxbin=gmxbin)
            
        elif simtype=='em':
            sim=gmx.EnergyMinimization(name,
                                       mdrun_options=mdrun_options, coord=start_coord, topology=self.path+'/topol.top',
                                       path_mdp=mdp, maxwarn=maxwarn,
                                       path_input=self.path,path_output=self.path,print_bash=True, gmxbin=gmxbin)
        #elif simtype=='wtmetad':
        #    sim=gmx.WTMD(name='wtmd')
        else:
            print("What are you trying to do?")
            exit()

        self.simulations.append(sim)

        return sim

    
    def print_command(self, bash_file):

        bash_file=self.path+'/'+bash_file

        with open(bash_file,'w') as f:
            f.write("#!/bin/bash\n\n")
            for s in self.simulations:
                f.write("# {}\n".format(s.name))
                f.write(s.bash_command+'\n\n')


    def _updateComposition(self):

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

    


    
    
    
