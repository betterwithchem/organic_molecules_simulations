import os
import sim_launch_py.utilities as util
import numpy as np


class Project():
    """The project class that stores and manage all the information and methods.\n

    Attributes:\n
    - name : name of the project. This defines the directory where the project will be stored and most of path variables. \n
    - project_path : absolute path of the project. The name of the directory is derived from the name of the project. \n
    - systems_path : directory where the data of the simulated systems are stored. This is the parent directory of the systems. \n
    - topology_path : directory where the topology files (.top and .itp) of the molecules used in the project are stored and retrieved to be used in each system.\n
    - init_struct_path : path where initial structures of the molecules used in the project are stored and retrieved to be used in each system.\n
    - mdp_path : path where default .mdp input files for gromacs are stored.\n
    - pickle_path : path where the project is saved in .pkl format.\n
    - job_script_path : path where the template files for submission scripts are stored.\n
    - systems : list of System() objects in the project.\n
    - molecules : list of Molecule() objects in the project.\n
    - gromacs : path to gromacs binary. \n
    - ambertools : path to ambertools binaries directory. \n
 
    Methods:\n
    - help() : print the help for this class.\n
    - add_system(self,name: str) : add and initialize System() to the Project(). \n
    - new_project(name=None, path=None, overwrite=False) : create new Project(). \n
    - save(self) : save Project() in Project().pickle_path. \n
    - load_project(project_folder: str) : load Project() stored in project_folde/.multisim.pkl. \n
    - write_sub_command(self,scriptname: str,hpc_system: str, template: str): write job submission/run scripts for the given HPC system for each System() and a global bash script to call all the job scripts. \n
    
     #### TODO METHODS #### \n
    - change_project_path(self,newpath: str) : when implemented, it will allow to copy the project to a new path and change all the pertinent path attributes.\n

    """
    
    def __init__(self, name: str, path: str):
        
        """Project Class Constructor

        :param name: name of the project.
        :type name: str
        :param path: path of the project
        :type path: str
        
        """
        
        if path.rstrip().endswith("/"):
            path=path.rstrip()[:-1]
        
        self._name=name
        self._systems=list()
        self._path=os.path.abspath(path)
        self._systems_path='./Systems'
        self._defaults_path='./Defaults'
        self._pickle_path="{}/.multisim.pkl".format(path)
        self._gromacs=None
        self._ambertools=None

        # look for gromacs and ambertools in the $PATH
        self._checkGromacs()
        self._checkAmberTools()

        

    @property
    def name(self):
        return self._name

    @property
    def systems(self):
        return self._systems

    @property
    def path(self):
        return self._path

    @property
    def systems_path(self):
        return self._systems_path

    @property
    def defaults_path(self):
        return self._defaults_path

    @property
    def pickle_path(self):
        return self._pickle_path

    @property
    def gromacs(self):
        return self._gromacs

    @property
    def ambertools(self):
        return self._ambertools

    @path.setter
    def path(self,p):
        self._path=p
        return self._path
    
    @systems_path.setter
    def systems_path(self,sp):
        self._systems_path=sp
        return self._systems_path

    @defaults_path.setter
    def defaults_path(self,dp):
        self._defaults_path=dp
        return self._defaults_path

    @pickle_path.setter
    def pickle_path(self,pp):
        self._pickle_path=pp
        return self._pickle_path

    @gromacs.setter
    def gromacs(self,gp):
        self._gromacs=gp
        return self._gromacs

    @ambertools.setter
    def ambertools(self,ap):
        self._ambertools=ap
        return self._ambertools

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
    - add_system(self,name: str) : add and initialize System() to the Project()
    - new_project(name=None, path=None, overwrite=False) : create new Project()
    
    

     #### TODO METHODS ####
    - change_project_path(self, newpath: str) : when implemented, it will allow to copy the project to a new path and change all the pertinent path attributes.
    - save(self) : save Project() in Project().pickle_path
    - load_project(project_folder: str) : load Project() stored in project_folde/.multisim.pkl. 
    - write_sub_command(self,scriptname: str,hpc_system: str, template: str): write job submission/run scripts for the given HPC system for each System() and a global bash script to call all the job scripts.""")

    def new_project(name: str, path: str, overwrite=False):
        
        """Add new project

        :param name: Project Name.
        :type name: str
        :param path: Project path.
        :type path: str
        :param overwrite: Overwrite previous project? Defaults to False.
        :type overwrite: bool, optional
        :returns nproject: new project

        """
        
        path=os.path.abspath(path)
        nproject = Project(path=path, name=name)

        if not os.path.exists(nproject._path):
            print("New Project: {}".format(nproject._name))
            print("=" * 100)
            util.create(nproject._path, arg_type='dir', backup=False)
        else:
            if overwrite:
                print("New Project: {}".format(nproject._name))
                print("=" * 100)
                util.create(nproject._path, arg_type='dir', backup=False)
            else:
                print("ERROR: Folder already exists.\n "
                      "You can change directory name or overwrite with 'overwrite=True' "
                      "(everything inside will be deleted)")
                
        util.create(nproject._path+'/'+nproject._systems_path, arg_type='dir')
        #util.create(nproject._init_struct_path, arg_type='dir')
        #util.create(nproject._topology_path, arg_type='dir')
        #util.create(nproject._mdp_path, arg_type='dir')

        os.chdir(path)
        
        return nproject

    def save(self):
        """
	Save project to project folder.
	"""
        print("\nSaving Project...", end="")
        import pickle
        import os
        if os.path.exists(self._pickle_path):
            os.rename(self._pickle_path, self._path+ "/.multisim.bck.pkl")
        with open(self._pickle_path, "wb") as file_pickle:
            pickle.dump(self, file_pickle)
            print("done")

    def load_project(project_folder: str):
        """Load an existing project. The Project object is saved in the project directory every time the command Project.save()

        :param project_folder: Location of the project to be loaded.
        :type project_folder: str
        :returns project: Loaded project.
        :rtype project: object

        """
        import pickle
        project_folder = os.path.realpath(project_folder)
        file_pickle = project_folder + "/.multisim.pkl"
        if os.path.exists(file_pickle):
            project = pickle.load(open(file_pickle, "rb"))
            print("Loading Project Name: {}\n".format(project._name))
            if os.path.realpath(project._path) != project_folder:
                project.path = project_folder
            os.chdir(project.path)
            return project
        else:
            print("No project found in '{}'. Use the 'Project.new_project' module to create a new project."
                  "".format(project_folder))


    def change_path(self, newpath: str, reset_program_paths: bool=False):
        """Change path of the project 
        :param newpath: new position of the project. The new position of the project will be newpath/Project.name 
        :type newpath: str
        :param reset_program_paths: Update programs (Gromacs, Ambertools, ...) paths? Defaults to False
        :type reset_program_paths: bool, optional

        """
        import shutil

        if newpath.rstrip().endswith("/"):
            newpath=newpath.rstrip()[:-1]
            
        newpath='{}/{}'.format(os.path.abspath(newpath),self._name)

        print('-'*50)
        print('Changing the path of the project\n from {}\n to   {}.'.format(self._path,newpath))
        print('This changes ONLY the path attribute of the project. It does NOT copy files to the new path!')

        self._path=newpath
        self._pickle_path="{}/.multisim.pkl".format(newpath)

        if reset_program_paths:
            print("Updating gromacs and antechamber binary paths.")
            self._checkGromacs()
            self._checkAmberTools()



    def add_system(self, name:str):
        
        """Add system to project

        :param name: system name.
        :type name: str
     
        """

        for s in self.systems:
            if name==s.name:
                print("ERROR: a system with name {} already exists".format(name))
                return

        util.create(self._systems_path+'/{}'.format(name), arg_type='dir')
        
        newsystem=System(name, path=self._systems_path+'/{}'.format(name))

        if self.gromacs:
            newsystem.gromacs=self.gromacs
        if self.ambertools:
            newsystem.ambertools=self.ambertools

        newsystem.index=len(self.systems)
            
        self.systems.append(newsystem)

        print("*"*50)
        print("Created system {}".format(newsystem.name))
        print("*"*50)

    def delete_system(self, delete_list):
        
        """Delete a system from project and update the index of remaining systems.

        :param delete_list: list of system indexes
        :type param: list

        """
        
        kept_systems=[]

        for sys in self.systems:
            if sys.index not in delete_list:
                kept_systems.append(sys)

        self.systems=kept_systems

        self._renumber_systems()

    def find_system_by_name(self,name):

        for s in self.systems:
            if s.name==name:
                return s

        print("System with name {} not found in the project.".format(name))
        return None

    def _renumber_systems(self):
        
        """Renumber starting from 0 the index of the systems in a project.

        """
        
        for i,s in enumerate(self.systems):
            s.index=i

        
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
           The value of the attribute is the absolute path of antechamber (/path/to/bin/antechamber).

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
    The system class that stores and manage all the information and methods.\n

    Attributes:\n
    - name : name of the system. This is used to create the directory of the system.\n
    - path : absolute path of the system.\n
    - molecules : list of Molecule() objects belonging to the system.\n
    - box : size of the box (of the initial configuration).\n
    - boxshape : shape of the simulation box.\n
    - simulations : list of Simulation() objects of the system.\n
    - run_command : command lines to be used to run the simulations.\n
    - gromacs : gromacs binary path, inherited from the Project(). \n
    - composition : composition of the system. Order of values is the same of the order in self.molecules.\n
    - atoms : list of Atom() objects in the system.\n
    - natoms : total numer of atoms in the systems.\n
    
    Methods:\n
    - help() : print the help for this class\n
    - add_molecule(self, name: str, moltype=None, knownmolecules=None) : add species to the system.\n
    - add_box(self, box_side: float, shape='cubic') : create simulation box.\n
    - createSolventBox(self, solvent: object, output_structure="solvent_box.pdb", density=None, nmols=None): add solvent molecules to the simulation box.  \n
    - insertSolute(self, solute: object, solvent: object, solvent_box="solvent_box.pdb", concentration=0, output_structure="start.pdb"): add solute molecules to the box and remove excess solvent molecules.\n
    - writeTop(self, atomtypes_path: str, molecules: objects): write the topology file for the system in gromacs format.\n
    - add_simulation(self, simtype: str, mdrun_options='', mdp='', print_bash=True, name='',maxwarn=0, start_coord='',gmxbin=''): add simulation to the system.\n
    - print_command(self, bash_file): print the bash script file to run the simulations.\n
   
    """
    
    def __init__(self, name: str, path=None):
        """System Class Constructor

        :param name: System name
        :type name: str
        :param path: System path. Defaults to None.
        :type path: str, optional

        """
        
        self._name=name
        self._index=None
        self._molecules=list()
        self._simulations=list()
        self._groups=list()
        self._box=None
        self._path=path
        self._topology=None
        self._gromacs=None
        self._ambertools=None
        self._atomtypes={}
        self._species={}
        self._last_saved_structure=None

    @property
    def name(self):
        return self._name

    @property
    def index(self):
        return self._index
    
    @property
    def molecules(self):
        return self._molecules

    @property
    def simulations(self):
        return self._simulations

    @property
    def groups(self):
        return self._groups

    @property
    def box(self):
        return self._box

    @property
    def path(self):
        return self._path

    @property
    def topology(self):
        return self._topology

    @property
    def gromacs(self):
        return self._gromacs

    @property
    def ambertools(self):
        return self._ambertools

    @property
    def atomtypes(self):
        return self._atomtypes

    @property
    def species(self):
        return self._species

    @property
    def last_saved_structure(self):
        return self._last_saved_structure
    
    @name.setter
    def name(self,n):
        self._name=n
    
    @index.setter
    def index(self,ndx):
        self._index=ndx

    @molecules.setter
    def molecules(self,ms):
        self._molecules=ms

    @box.setter
    def box(self,box_vector):
        self._box=box_vector

    @path.setter
    def path(self,p):
        self._path=p

    @topology.setter
    def topology(self,tp):
        self._topology=tp

    @gromacs.setter
    def gromacs(self,gp):
        self._gromacs=gp

    @ambertools.setter
    def ambertools(self,ap):
        self._ambertools=ap

    @last_saved_structure.setter
    def last_saved_structure(self,pdb):
        self._last_saved_structure=pdb
        

    @staticmethod
    def help():
        print("""The system class that stores and manage all the information and methods.

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
    - print_command(self, bash_file): print the bash script file to run the simulations.""")
    
    def add_molecule(self, name: str, structure_file: str=None, keep_coordinates: bool=True, keep_box: bool=True):
        """Add new molecules to a system starting from a structure file (e.g. a PDB file).

        :param name: name of the molecule. 
        :type name: str
        :param structure_file: name of the structure file to be read. Defaults to None.
        :type structure_file: str, optional
        :param keep_coordinates: keep the coordinates read from the structure file. Defaults to True.
        :type keep_coordinates: bool, optional
        :param keep_box: read the box parameters from the structure file and use them to define/update the box attribute of the system. Defaults to True.
        :type keep_box: bool, optional

        """
        import subprocess
        
        if not structure_file:
            print('Error: a structure file is required')
            return

        curdir=os.path.abspath(os.path.curdir)
        os.chdir(self.path)
                
        # convert structure_file to mol2

        extension=structure_file.split('.')[-1]
        basename_structure_file=os.path.basename(os.path.splitext(structure_file)[0])
                
        mol2file=self.path+'/{}'.format(basename_structure_file+'.mol2')

        result=subprocess.run('{3} -i {0} -fi {1} -o {2} -fo mol2 -j 5 -dr no -at gaff2'.format(structure_file,
                                                                                                extension,
                                                                                                basename_structure_file+'.mol2',
                                                                                                self.ambertools),
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True)


        os.chdir(curdir)

        self._loadfrommol2(name,mol2file,keep_coordinates=keep_coordinates)

        found_box=False
        if keep_box:
            if extension=='pdb':
                with open(structure_file,'r') as f:
                    for line in f:
                        if line.startswith('CRYST1'):
                            found_box=True
                            self.box=[0. for i in range(6)]
                            self.box[0]=float(line[6:15])/10
                            self.box[1]=float(line[15:24])/10
                            self.box[2]=float(line[24:33])/10
                            self.box[3]=float(line[33:40])
                            self.box[4]=float(line[40:47])
                            self.box[5]=float(line[47:54])
                            break
            else:
                print('ERROR: only PDB files are supported to include box data from structure files.' \
                      'You can add manually box data with system.add_box().')
        

        for m in self.molecules:
            if m.resname not in self.species:
                self.species[m.resname]={}

        self._renumber_atoms()



    def insert_molecules(self, name: str, molstruct: str, initial_conf: str=None, final_conf: str='inserted.pdb', nmol: int=1):
        """Insert molecules in random positions.

        :param name: name of the molecule.
        :type name: str
        :param molstruct: structure file of the molecule to insert. PN: it supports only structures with a single molecule. 
        :type molstruct: str
        :param initial_conf: structure file into which the new molecules are inserted. Defaults to None. If no initial configuration is given, an empty box is created and populated. System.box need to have been already defined.
        :type initial_conf: str, optional
        :param final_conf: final structure file (in PDB format) with the inserted molecules. Defaults to inserted.pdb.
        :type final_conf: str, optional
        :param nmol: number of molecules to insert. Defaults to 1.
        :type nmol: int, optional
        
        """
        if len(self.molecules)>0 and initial_conf==None:
            print("ERROR: you are trying to insert new molecules without taking into account molecules" \
                  "that are already in the system.")

            return

        
        final_conf_ext=final_conf.split('.')[-1]
        if final_conf_ext!='pdb':
            print('ERROR while inserting new molecules: Only PDB format is accepted for final_conf')
            exit()

        directory,filename=os.path.split(final_conf)
        if os.path.isabs(final_conf):
            if directory==self.path:
                pass
            else:
                print('Warning: path given for the final configuration is not in {0}. Path has been changed to {0}.'.format(self.path))
                final_conf=self.path+'/{}'.format(filename)
        else:
            final_conf=self.path+'/{}'.format(filename)

        
        if not initial_conf:

            initial_conf=self.path+'/{}'.format('TEMP_EMPTY_BOX.pdb')
            if self.box==None:
                print('ERROR: when adding molecules without an initial configuration, a box must have been already defined')
                return

            with open('{}'.format(initial_conf),'w') as f:
                f.write('TITLE	MOLECULES IN A BOX\n')
                f.write('{:6s}{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}'.format("CRYST1",
                                                                                 self.box[0]*10,  # sides in A
                                                                                 self.box[1]*10,
                                                                                 self.box[2]*10,
                                                                                 self.box[3],	  # anglse in degrees
                                                                                 self.box[4],
                                                                                 self.box[5]))
        else:

            if os.path.isabs(initial_conf):
                if os.path.isfile(initial_conf):
                    pass
                else:
                    print('ERROR: file {} does not exist'.format(initial_conf))
                    return
            elif not os.path.isabs(initial_conf):
                if os.path.isfile(self.path+'/{}'.format(initial_conf)):
                    initial_conf=self.path+'/{}'.format(initial_conf)
                elif  os.path.isfile(initial_conf):
                    pass
                else:
                    print('ERROR: file {} does not exist'.format(initial_conf))
                    return
                    

            

        import subprocess        
        result=subprocess.run("{0} insert-molecules -f {1} -ci {2} -nmol {3} -o {4} -try 10000".format("gmx",
                                                                                            initial_conf,
                                                                                            molstruct,
                                                                                            nmol,
                                                                                            final_conf),
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True)

        # now check that all molecules have been inserted
        added_line=result.stdout.decode()[result.stdout.decode().find("Added "):result.stdout.decode().find("Added ")+result.stdout.decode()[result.stdout.decode().find("Added "):].find(')')+1]

        added=int(added_line.split()[1])
        if added != nmol:
            print("ERROR! {}".format(added_line))
            exit()

        self._last_saved_structure=final_conf


        # now add/create Molecule and Atom objects

        index_first_new_mol=len(self.molecules)
        index_last_new_mol=index_first_new_mol+nmol-1

        #print(len(self.molecules))
        self.add_molecule(name,structure_file=molstruct,keep_coordinates=False,keep_box=False)
        #print(len(self.molecules))

        if nmol>1:
            lastmol=len(self.molecules)-1
            import copy
        
            for i in range(1,nmol):
                new_molecule=copy.deepcopy(self.molecules[lastmol])
                self.molecules.append(new_molecule)

        else:
            new_molecule=self.molecules[-1]

        self._update_molecule_indexes()
        self._update_composition()

        # and update coordinates

        #print('Update coordinates for molecules from {} to {} reading from file {}.'.format(index_first_new_mol,index_last_new_mol, final_conf))        
        self._update_coordinates(final_conf, start=index_first_new_mol, end=index_last_new_mol)
        self._renumber_atoms()
        
        print("{} molecules of residue {} have been added. The file {} has been created.".format(nmol,new_molecule.resname, os.path.relpath(final_conf)))


        
    def create_group(self, name: str, atoms: list=None,molecules: list=None):
        """Create a new group of atoms

        :param name: Name of the new group 
        :type name: str
        :param atoms: list of atoms to add to the group. Defaults to None 
        :type atoms: list, optional
        :param molecules: list of molecules to add to the group (its atoms will be added to the group). Defaults to None. 
        :type molecules: list, optional
        :returns: 

        """

        for g in self.groups:
            if g.name==name:
                print("ERROR: group with name {} already exists in system {}".format(name,self.name))
                return
        
        newgroup=Group(name)
        if atoms is not None:
            for atom in atoms:
                newgroup.atoms.append(atom)

        if molecules is not None:
            for molecule in molecules:
                for atom in molecule.atoms:
                    newgroup.atoms.append(atom)

        if newgroup.natoms>0:
            newgroup._check_duplicates()
        
        self.groups.append(newgroup)
        

    def add_to_group(self, name: str, atoms: list=None, molecules: list=None):
        """Add atoms to an existing group.

        :param name: Name of the group to which atoms will be added. 
        :type name: str
        :param atoms: List of atoms to add. Defaults to None.
        :type atoms: list, optional
        :param molecules: List of molecules to add. Defaults to None.
        :type molecules: list, optional
        :returns: 

        """

        gindex=self.find_group_by_name(name)
        g=self.groups[gindex]

        if atoms is not None:
            for atom in atoms:
                g.atoms.append(atom)

        if molecules is not None:
            for molecule in molecules:
                for atom in molecule.atoms:
                    g.atoms.append(atom)

        g._check_duplicates()

    


    def find_group_by_name(self, name: str):
        """Get the index of a group in a system given its name.

        :param name: Name of the index
        :type name: str
        :returns group_index: index of the group
        :rtype group_index: int 

        """

        for group_index,g in enumerate(self.groups):
            if g.name==name:
                return group_index

        print("Warning: Group with name {} not found in system {}.".format(name,self.name))
        return None


    def delete_group(self, name: str):
        """Delete a group from the system

        :param name: Name of the group to eliminate.
        :type name: str
        :returns: 

        """

        kept_groups=[]
        for g in self.groups:
            if g.name!=name:
                kept_groups.append(g)

        self.groups=kept_groups


    def replicate_cell(self,repl: list=[1, 1, 1]):
        """Replicate the given box in the 3 directions.

        :param repl:  number of times the cell is replicated in each direction. Defaults to [1, 1, 1]
        :type repl: list, optional

        """

        import copy
        
        M=self._box_matrix(self.box)
        M_inv=np.linalg.inv(M)

        newmolecules=[]

        for xmult in range(0,repl[0]):
            for ymult in range(0,repl[1]):
                for zmult in range(0,repl[2]):
                    if any(np.array(repl)>0):
                        factor=[xmult, ymult, zmult]
                        if np.any([xmult,ymult,zmult]):
                            for im,m in enumerate(self.molecules):
                                new_mol=copy.deepcopy(m)
                                for ia,a in enumerate(new_mol.atoms):

                                    fract_coords=M_inv.dot(a.coordinates)
                                    new_fract_coords=[factor[i]*fract_coords[i] for i in range(3)]

                                    #new_cart_coords=M.dot(new_fract_coords)
                                    new_cart_coords=[0., 0., 0.]
                                    for idim,rf in enumerate(factor):
                                        #print(im,ia,factor)
                                        new_cart_coords[idim]=a.coordinates[idim]+\
                                            xmult*self.box[0]*M[0][idim]+\
                                            ymult*self.box[1]*M[1][idim]+\
                                            zmult*self.box[2]*M[2][idim]

                                    a.coordinates=new_cart_coords
                                    #print(im,ia,a.coordinates)
                                    #print(new_mol.atoms[ia].coordinates)

                                newmolecules.append(new_mol)

        for i in range(3):
            self.box[i]*=repl[i]
        
        
        self.molecules+=newmolecules
        self._update_molecule_indexes()
        self._update_composition()
        

    def create_topology(self, topology: str='topology.top'):
        """Create topology for the current system

        :param topology: name of the resulting topology file. Defaults to 'topology.top'.
        :type topology: str, optional

        """

        #1)

        from sim_launch_py import ffparms
        
        for species in self.species:
           
            atomtypes=ffparms.extract_atomtypes_from_gmx_top(self.species[species]['top'])
            for atype in atomtypes:
                if atype not in self.atomtypes:
                    self.atomtypes.update({atype:atomtypes[atype]})
                else:
                    print('atom type {} already existing, continuing...'.format(atype))

            self.species[species]['itp']=ffparms.extract_molecule_from_gmx_top(self.species[species]['top'])
            
        #2)
        # CAVEAT: the order of the molecules in the structure file is NOT checked.
        self._update_composition()

        #3)
        topology='{}/{}'.format(self.path,topology)
        with open(topology,'w') as top:
            top.write(";\n\
; Topology for system {}\n\
;\n\n".format(self.name))

            top.write("[ defaults ]\n\
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n\
1               2               yes             0.5          0.83333333  \n\
\n\
[ atomtypes ]\n\
; name    at.num    mass    charge ptype  sigma      epsilon\n")

            for atype in self.atomtypes:
                top.write('{}\t{}\n'.format(atype,self.atomtypes[atype]['parms']))

            top.write('\n')

            for species in self.species:
                top.write(self.species[species]['itp'])
                top.write('\n\n')

            top.write('[ system ]\n{}\n\n'.format(self.name))

            top.write('[ molecules ]\n')

            for species in self.species:
                top.write('{}\t{}\n'.format(species,self.species[species]['nmols']))

        self.topology='{}'.format(topology)

        
    def _renumber_atoms(self):

        """Assign the absolute index to each atom of the system

        """

        index=0
        for im in range(len(self.molecules)):
            for a in self.molecules[im].atoms:
                a.absindex=index
                index+=1

         

    def _update_composition(self):

        """Update the 'nmols' attribute of self.species dictionary by checking the residue names of the molecules in the system

        """
        found_species={}
        for m in self.molecules:
            if m.resname not in found_species:
                found_species[m.resname]=1
            else:
                found_species[m.resname]+=1

        for sp in found_species:
            self.species[sp].update({'nmols':found_species[sp]})

        
         
        
        
    def add_simulation(self, name: str, simtype: str, simulation_dict: dict=None):

        """Add a simulation to the current system object.

        :param name: Name of the simulation object
        :type name: str
        :param simtype: Type of simulation. Available types are 'em' and 'md'.
        :type simtype: str
        :param simulation_dict: dict with simulation specific parameters. Defaults to None. Values in the dict are: 'coordinates', 'mdrun_options', 'topology', 'path_mdp', 'maxwarn', 'tpr', 'nsteps', 'plumed'.
        :type simulation_dict: dict, optional
 
        """
        import sim_launch_py.gromacs as sim
        import shutil

        accepted_simtypes=['em','md','posre']

        if simtype not in accepted_simtypes:
            print("ERROR: provided simtype, {}, is not valid. Acceptable values are {}.".format(simtype,accepted_simtypes))
            return

        util.create('{}/{}'.format(self.path,name), arg_type='dir', backup=False)

        

        newsim_index=len(self.simulations)
        newsim_path='{}/{}'.format(self.path,name)
        
        
        if simtype=='em':
            newsim=sim.EnergyMinimization(name, newsim_index)

        elif simtype=='md':

            newsim=sim.MD(name, newsim_index)
            newsim.plumed=simulation_dict.pop('plumed',None)

        elif simtype=='posre':
            
            newsim=sim.Posre_MD(name, newsim_index)
            try:
                newsim.posre=simulation_dict.pop('posre')
            except KeyError:
                print('ERROR: you tried to add a simulation with position restraints but no'\
                      'reference structure for the restraint has been passed')
                return
            newsim.posre=os.path.relpath(newsim.posre,start=newsim_path)


        
        newsim.path=newsim_path
        newsim.state='Pending'
        
        
        try:
            newsim.coordinates=simulation_dict.pop('coordinates')
        except KeyError:
            print("ERROR: unknown starting configuration, please specify a structure file.")
            return
        else:
            print("Simulation {}(type {}) will start from configuration {}.".format(name,simtype,newsim.coordinates))      
        newsim.mdrun_options=simulation_dict.pop('mdrun_options','-v')
        newsim.topology=simulation_dict.pop('topology',self.topology)
        newsim.mdp=simulation_dict.pop('path_mdp',None)
        newsim.maxwarn=simulation_dict.pop('maxwarn',0)
        newsim.tpr=simulation_dict.pop('tpr',None)
        newsim.nsteps=simulation_dict.pop('nsteps',None)

        shutil.copy(newsim.mdp,newsim.path)
        newsim.mdp=os.path.basename(newsim.mdp)

        shutil.copy(newsim.topology,newsim.path)
        newsim.topology=os.path.basename(newsim.topology)

        #shutil.copy(newsim.coordinates,newsim.path)
        newsim.coordinates=os.path.relpath(newsim.coordinates,start=newsim.path)

        if newsim.simtype=='md' and newsim.plumed:
            newsim.mdrun_options+=' -plumed {}'.format(newsim.plumed)

        
        if len(simulation_dict)>0:
            print("Warning: While creating the simulation object, there are remaining keyword arguments that have not been used: {}".format(simulation_dict.keys()))

        self.simulations.append(newsim)

        print("Added simulation {} ({}) to system {} ({}).".format(name,newsim_index,self.name,self.index))


    def create_run_script(self, scriptname: str, simulations: list=None, platform: str='bash', platform_dict: dict=None):

        """Create the script to run the simulations.

        :param scriptname: name of the final script to launch the simulations.
        :type scriptname: str 
        :param simulations: which simulations need to be run. Defaults to None (which means all simulations. I know...)
        :type simulations: list
        :param platform: type of platform on which simulations will be run. Available are 'bash' for running locally, 'myriad' to run on Myriad HPC @ UCL, and 'archer' to run on Archer2 HPC. Defaults to 'bash'. 'myriad' is a good starting point for GridEngine-like submission scripts, 'archer' is a good starting point for SLURM submission scripts.
        :type platform: str, optional
        :param platform_dict: dict with platform specific parameters. Defaults to None.
        :type platform_dict: dict, optional

        """
        if simulations is None:
            simulations=[sim.index for sim in self.simulations] 

        available_platforms=['bash','myriad','archer','custom']
        if platform not in available_platforms:

            print("ERROR: platform to run the simulations is not valid. Available values are {}.".format(available_platforms))
            return
            
        with open('{}/{}'.format(self.path,scriptname),'w') as f:

            if platform=='myriad':

                header="""#!/bin/bash -l

# Batch script to run an MPI parallel job under SGE with Intel MPI.

# Request two hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt={0}

# Request 1 gigabyte of RAM per process (must be an integer followed by M, G, or T)
#$ -l mem=1G

# Request 1 gigabyte of TMPDIR space per node 
# (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=2G

# Set the name of the job.
#$ -N {1}

# Select the MPI parallel environment and 36 processes.
#$ -pe mpi {2}

# Set the working directory to somewhere in your scratch space.
#$ -cwd 

# load gromacs
module unload -f compilers mpi
module load compilers/intel/2018/update3 mpi/intel/2018/update3/intel libmatheval flex plumed/2.5.2/intel-2018 gromacs/2019.3/plumed/intel-2018

export OMP_NUM_THREADS={3}

""".format(platform_dict['wallclock'],platform_dict['job_name'],platform_dict['mpi'],platform_dict['omp'])

                mpirun=platform_dict.get('mpirun','gerun')
                                
                f.write(header)


            elif platform=='archer':

                header="""#!/bin/bash --login

#SBATCH --job-name={0}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={1}
#SBATCH --cpus-per-task=1
#SBATCH --time={2}

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account={4}
#SBATCH --partition=standard
#SBATCH --qos=standard

# Load the gromacs module 
module load gromacs/2022.4+plumed

# Recommended environment settings
export OMP_NUM_THREADS={3}
# Ensure the cpus-per-task option is propagated to srun commands
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# srun launches the parallel program based on the SBATCH options

""".format(platform_dict['job_name'],platform_dict['mpi'],platform_dict['wallclock'],platform_dict['omp'],platform_dict['budget'])

                mpirun=platform_dict.get('mpirun','srun')
                f.write(header)

            elif platform=='custom':

                with open(platform_dict['template_file'],'r') as temp_f:
                    for line in temp_f:
                        f.write(line)

                mpirun=platform_dict.get('mpirun','')
                gmx=platform_dict.get('gmx_bin','gmx')
                
            for sim_index in simulations:

                sim=self.simulations[sim_index]

                f.write('cd {}\n'.format(os.path.basename(sim.path)))

                gmx=platform_dict.get('gmx_bin','gmx_mpi')

                if sim.mdp:
                     f.write('{0} grompp -f {1} -o {2}.tpr -maxwarn {3} -p {4} -c {5} '.format(gmx,
                                                                                              sim.mdp,
                                                                                              sim.name,
                                                                                              sim.maxwarn,
                                                                                              sim.topology,
                                                                                              sim.coordinates))

                if sim.simtype=='posre':
                    f.write('-r {} \n'.format(sim.posre))
                else:
                    f.write('\n')


                f.write('{0} {1} mdrun -deffnm {2} {3} '.format(mpirun, gmx, sim.name,sim.mdrun_options))

                if sim.nsteps:
                    f.write('-nsteps {}'.format(sim.nsteps))
                     
                f.write('\ncd ..\n\n')

            

        print('Written run script file in {}/{}'.format(self.path,scriptname))
 

    def find_simulation_by_name(self,name: str):

        """Find simulations by their name

        :param name: Name of the simulation.
        :type name: str
        :returns simulation: Simulation object
        :rtype simulation: object

        """
        sim_list=[]
        for simulation in self.simulations:
            if simulation.name==name:
                return simulation

    def find_simulations_by_type(self,simtype: str):

        """Find simulations by their type

        :param simtype: Type of simulation (e.g. 'em', 'md').
        :type simtype: str
        :returns: List of simulation indexes.
        :rtype: list.

        """
        sim_list=[]
        for sim in self.simulations:
            if sim.type==simtype:
                sim_list.append(sim.index)

        return sim_list

    @staticmethod
    def _box_matrix(box: list):
        """Given a box vector [a, b, c, alpha, beta, gamma], compute the box matrix.

        :param box: box vector. box vectors (a,b,c) are assumed to be i nm, angles (alpha, beta, gamma) are assumed to be in degrees.
        :type box: list
        :returns M: box matrix
        :rtype M: list

        """

        deg2rad=np.pi/180
        radbox=[0,0,0,0,0,0]
        for i in range(3,6):
            radbox[i]=box[i]*deg2rad
            

        n2=(np.cos(radbox[3])-np.cos(radbox[5])*np.cos(radbox[4]))/np.sin(radbox[5])
        M=np.array([[1, 0, 0],
                    [np.cos(radbox[5]), np.sin(radbox[5]), 0],
                    [np.cos(radbox[4]), n2,             np.sqrt(np.sin(radbox[4])**2-n2*n2)]])

        return M
        
    def delete_molecule(self,delete_list: list):
        """Delete molecules from a system and update the index of the remaining molecules.

        :param delete_list: list of the indexes of the molecules to be deleted.
        :type delete_list: list
        
        """
        
        kept_molecules=[]

        print("Deleting {} molecules from system {}.".format(len(delete_list),self.name))
        for im,m in enumerate(self.molecules):
            if m._index not in delete_list:
                kept_molecules.append(m)

        self.molecules=kept_molecules

        self._update_molecule_indexes()
        print("Now system {} has {} molecules".format(self.name,len(self.molecules)))

        self._update_composition()

        
    def find_molecule_by_resname(self,resname: str):
        """Find the molecules in a system with the given residue name.

        :param resname: residue name of the molecules to be found.
        :type resname: str
        :returns mol_list: list of indexes of the molecules with the given residue name.
        :rtype mol_list: list

        """

        mol_list=[]
        for m in self.molecules:
            if m.resname==resname:
                mol_list.append(m.index)

        return mol_list


    def save_pdb(self,pdb_filename: str):
        """Save a PDB file with the current configuration of the system.

        :param pdb_filename: output PDB file. 
        :type pdb_filename: str

        """
        
        f=open(self.path+'/{}'.format(pdb_filename),'w')

        if self.box:
            f.write('{:6s}{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n'.format("CRYST1",
                                                                             self.box[0]*10,  # sides in A
                                                                             self.box[1]*10,
                                                                             self.box[2]*10,
                                                                             self.box[3],	  # anglse in degrees
                                                                             self.box[4],
                                                                             self.box[5]))

        iatom=1
        for m in self.molecules:
            for a in m.atoms:

                # pdb files have a limit in the representation of number of atoms to 99999 and 9999 residues.
                # to avoid formatting problems here iatom are cut to the last 5 digits and m.index to the last 4 digits

                if iatom>99999:
                    flase_iatom=iatom%100000
                else:
                    false_iatom=iatom

                if m.index+1>9999:
                    false_mindex=(m.index)%10000
                else:
                    false_mindex=m.index+1


                    
                f.write("{:6s}{:5d} {:4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:2s}  \n".format("ATOM",
                                                                                                                         false_iatom,
                                                                                                                         a.name,
                                                                                                                         m.resname,
                                                                                                                         false_mindex,
                                                                                                                         a.coordinates[0]*10,
                                                                                                                         a.coordinates[1]*10,
                                                                                                                         a.coordinates[2]*10,
                                                                                                                         1.00,
                                                                                                                         1.00,
                                                                                                                         a.element))
                iatom+=1

        self.last_saved_structure=f.name
        
        f.close()


        
        

    def add_box(self, box_side: float, shape: str='cubic'):
        """ Create simulation box.
 
        :param box_side: side of the box in nm.
        :type box_side: float
        :param shape: shape of the box. Defaults to 'cubic'.
        :type shape: str, optioal

        """

        import numpy as np
        self.boxshape=shape
        
        if shape=='cubic':
            self.box=[box_side]*3+[90.0]*3
        elif shape=='dodecahedron':
            self.box=[box_side, box_side, 0.5*np.sqrt(2)*box_side] + [60.0]*3 # see Gromacs manual
        elif shape=='octahedron':
            self.box=[box_side, 2/3*np.sqrt(2)*box_side, 1/3*np.sqrt(6)*box_side]+[71.53, 109.47, 71.53] # see Gromacs manual
        else:
            # custom shape
            if len(box_side)!=6:
                print("ERROR: for custom boxes, box_side must be a list of 6 values ([a,b,c,alpha,beta,gamma])")
                exit()
            self.box=box_side

    def center_box(self):

        """Center the atoms in the box with gmx editconf (equivalent to gmx editconf -c)

        """

        import subprocess

        self.save_pdb('TEMP_STRUCT_TO_BE_CENTERED.pdb')
        result=subprocess.run('{0} editconf -f {1}/TEMP_STRUCT_TO_BE_CENTERED.pdb -o {1}/CENTERED_STRUCT.pdb -c'.format(self.gromacs,self.path),stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True)

        self._last_saved_structure=self.path+'/CENTERED_STRUCT.pdb'
        
        self._update_coordinates(self.path+'/CENTERED_STRUCT.pdb')

        
    def _update_coordinates(self,structure_file,start=0,end=np.inf):
        """Assign the coordinates of atoms from a PDB structure file.

        :param structure_file: PDB file from which coordinates are read.
        :type structure_file: str
        :param start: index of the first molecule of the system for which atomic coordinates are updated. Defaults to 0. 
        :type start: int, optional
        :param end: index of the last molecule (included) of the system for which atomic coordinates are updated. Defaults to np.inf.
        :type end: int, optional
        
        """

        iresid=0
        iatom=0

        prev_resid=None
        prev_iresid=None
        with open(structure_file,'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    resid=int(line[22:26])-1


                    if prev_resid==None:
                        prev_resid=resid
                        #print(prev_resid,resid,iatom)
                    if not prev_iresid:
                        prev_iresid=iresid
                        
                    if resid!=prev_resid:
                        iresid+=1
                        prev_resid=resid
                        iatom=0


                    if iresid>=start and iresid<=end:

                        #print(iresid,resid, self.molecules[iresid].index,self.molecules[iresid].resname, len(self.molecules[iresid].atoms), iatom)
                        newcoords=[float(line[30:38])/10, float(line[38:46])/10, float(line[46:54])/10]
                        self.molecules[iresid].atoms[iatom].coordinates=newcoords

                        iatom+=1
                        prev_iresid=iresid

                    


                                            
                    

    def _loadfrommol2(self,name: str, mol2file: str, keep_coordinates: bool=True):
        """Load molecules from a MOL2 structure file.

        :param name: name of the MOL2 structure file.
        :type name: str
        :param mol2file: name of the MOL2 structure file.
        :type mol2file: str
        :param keep_coordinates: keep coordinates read from the structure file. Defaults to True.
        :type keep_coordinates: bool, optional
        
        """
        
        with open(mol2file,'r') as f:

            read_coordinates=False
            read_bonds=False
           
            molecules=[]
            atoms=[]
            
            for line in f:
                if line.strip().endswith('<TRIPOS>ATOM'):
                    read_coordinates=True
                elif line.strip().endswith('<TRIPOS>BOND'):
                    read_bonds=True
                    natoms=len(atoms)
                    #bond_matrix=dia_array((natoms,natoms),dtype='bool')

                if line.startswith('@<TRIPOS>') and not line.strip().endswith('ATOM'):
                    read_coordinates=False
                if line.startswith('@<TRIPOS>') and not line.strip().endswith('BOND'):
                    read_bonds=False

                if read_coordinates and not line.strip().endswith('<TRIPOS>ATOM'):
                    cols=line.split()
                    index=int(cols[0])-1
                    name=cols[1]
                    if keep_coordinates:
                        coordinates=[float(cols[2])/10,float(cols[3])/10,float(cols[4])/10]
                    else:
                        coordinates=[0.0, 0.0, 0.0]
                    atomtype=cols[5]
                    resid=int(cols[6])  # CAVEAT: one residue per molecule
                    resname=cols[7]

                    atoms.append(Atom(name,
                                      index=index,
                                      atomtype=atomtype,
                                      resname=resname,
                                      coordinates=coordinates,
                                      resid=resid))

                elif read_bonds and not line.strip().endswith('<TRIPOS>BOND'):

                    #print(line)

                    cols=line.split()
                    i_index=int(cols[1])-1                  
                    j_index=int(cols[2])-1
                    atoms[i_index].bonds.append(j_index)
                    for ia,a in enumerate(atoms):
                       if j_index in a.bonds and ia not in atoms[j_index].bonds:
                            atoms[j_index].bonds.append(ia)                            
                    #bond_matrix[i_index, j_index]

        self._assign_atoms_to_molecules(name,atoms)
        

    def _assign_atoms_to_molecules(self, name: str, atoms: list):
        """Add atom objects to molecule objects.

        :param name: name of the molecule.
        :type name: str
        :param atoms: list of atom objects.
        :type atoms: list
        
        """

        molecules=[]
        for ia,a in enumerate(atoms):

            if ia==0:
                newmolecule=Molecule(name,resname=a.resname,index=a.resid)
                prev_id=a.resid
            elif a.resid!=prev_id:
                molecules.append(newmolecule)
                newmolecule=Molecule(name,resname=a.resname,index=a.resid)
                prev_id=a.resid

            newmolecule.atoms.append(a)     
        molecules.append(newmolecule)

        for im,m in enumerate(molecules):
            m._renumber_atoms()
            self._molecules.append(m)
          


    def _update_molecule_indexes(self):
        """Renumber the indexes of the molecules (0-based).
        """

        for im,m in enumerate(self.molecules):
            
            m.index=im
            m._update_atom_resid()

        for g in self.groups:
            g._check_duplicates()
                    
        

class Molecule():
    """The molecule class that stores and manage all the information and methods.

    Attributes:\n
       - name : name of the molecule type\n
       - resname : residue name\n
       - index : index of the molecule\n
       - mw : molecular weight in [Da]\n
       - atoms : list of atom types of the molecule type\n
       - contact_matrix : NxN matrix (N=number of atoms in the molecule) of the bonds between the atoms.\n

    Methods:\n
       - help(): print the help for this class.

    """

    def __init__(self, name: str, resname='UNK', index=None):

        """Molecule Class Constructor

        :param name: name of the molecule
        :type name: str
        :param resname: 3-letter residue name of the molecule
        :type resname: 
        :param index: index of the molecule
        :type index: 

        """
        
        self._name=name
        self._resname=resname
        self._index=index
        self._mw=None
        self._atoms=list()
        self._contact_matrix=None
        self._com=None

    @property
    def name(self):
        return self._name

    @property
    def resname(self):
        return self._resname

    @property
    def index(self):
        return self._index

    @property
    def mw(self):
        if self._mw==None:
            self._calc_mass()
        return self._mw

    @property
    def atoms(self):
        return self._atoms

    @property
    def contact_matrix(self):
        if not np.all(self._contact_matrix):
            self._compute_contact_matrix()
        return self._contact_matrix

    @property
    def com(self):
        if np.any(self._com==None):
            self._get_com()
        return self._com

    @resname.setter
    def resname(self,rn):
        self._resname=rn
        return self._resname

    @index.setter
    def index(self,n):
        self._index=n
        return self._index

    def _calc_mass(self):
        """Compute the mass of the molecule and store it in self._mw.
        """

        mw=0
        for a in self.atoms:
            el=a.element
            amw=a.mw
            mw+=amw
            
        self._mw=mw

    def _get_com(self):

        self._com=np.array([0., 0., 0.])
        for a in self.atoms:
            self._com+=a.coordinates

        self._com/=len(self.atoms)
        
    def _renumber_atoms(self):
        """Renumber the atom indexes in the molecule starting from 0.
        """

        base_index=self._atoms[0]._index
        
        for ia,a in enumerate(self._atoms):
            newbonds=[]
            for b in a._bonds:
                newbonds.append(b-base_index)
            a._index-=base_index

            a.bonds=newbonds
            

    def _update_atom_resid(self):
        """Renumber the molecular index of the atoms (i.e. the index of the molecule they belong to).
        """

        for ia,a in enumerate(self._atoms):
            a._resid=self.index


    def _compute_contact_matrix(self):
        """Compute the matrix of bonds of the molecule.
        """

        import numpy as np

        n_atoms=len(self._atoms)
        
        contact_matrix=np.zeros((n_atoms,n_atoms))
        
        for ia,a in enumerate(self._atoms):
            for ja in a._bonds:
                contact_matrix[ia,ja]=1
                contact_matrix[ja,ia]=1
                

        self._contact_matrix=contact_matrix
    
class Atom():
    """
    The atom class that stores and manage all the information and methods.

    Attributes:\n
       - name (str) : name of the atom.\n
       - index (int) : atom index.\n
       - atomtype (str) : atom type.\n
       - resname (str) : name of the residue the atom is part of.\n
       - coordinates (list): coordinates of the atom in nm.\n
       - resid (int) : index of the molecule the atom is part of.\n
       - element (str) : chemical element of the atom.\n
       - mw (float) : atomic mass in Da.\n
       - bonds (list) : list of the index of the atoms in the same molecule that the atom is bound to.\n

    Static Methods:\n
       - _get_element(atomtype) : Assign atomic element of the atom based on its name.\n
       - _get_mw(element) : Assign mass to the atom based on its atomic element.\n
       
    """

    def __init__(self, name: str, index: int=None, atomtype: str=None, resname: str=None, coordinates: list=None, resid: int=None):
        """Atom Class Constructor

        :param name: name of the atom.
        :type name: str
        :param index: atom index. Defaults to None.
        :type index: int, optional
        :param atomtype: atom type. Defaults to None.
        :type atomtype: str, optional
        :param resname: residue name (resname) the atom is part of.  Defaults to None.
        :type resname: str, optional
        :param coordinates: coordinates of the atom in nm.  Defaults to None.
        :type coordinates: list, optional
        :param resid: index of the molecule the atom is part of.  Defaults to None.
        :type resid: int, optional

        """
        
        self._name=name
        self._index=index
        self._absindex=None
        self._atomtype=atomtype
        self._resname=resname
        self._coordinates=coordinates
        self._resid=resid
        self._element=None
        self._mw=None        
        self._bonds=list()

    @property
    def name(self):
        return self._name

    @property
    def index(self):
        return self._index

    @property
    def absindex(self):
        return self._absindex

    @property
    def atomtype(self):
        return self._atomtype
    
    @property
    def element(self):
        if self._element==None:
            #self._name
            self._element=self._get_element(self.atomtype)
        return self._element

    @property
    def mw(self):
        if self._mw==None:
            self._mw=self._get_mw(self.element)
        return self._mw

    @property
    def bonds(self):
        return self._bonds

    @property
    def resname(self):
        return self._resname

    @property
    def resid(self):
        return self._resid

    @property
    def coordinates(self):
        return self._coordinates

    @index.setter
    def index(self,n):
        self._index=n

    @absindex.setter
    def absindex(self,an):
        self._absindex=an

    @resid.setter
    def resid(self,r):
        self._resid=r

    @bonds.setter
    def bonds(self,b):
        self._bonds=b

    @coordinates.setter
    def coordinates(self,x):
        self._coordinates=x

        
    @staticmethod
    def _get_element(atomtype):
        """Assign atomic element of the atom based on its name.
        """

        if atomtype[0].upper() in ['H','C','N','O','F','S']:
            return atomtype[0].upper()
        elif atomtype[0].upper()+atomtype[1].lower() in ['Cl']:
            return atomtype[0].upper()+atomtype[1].lower()
        else:
            print("ERROR: Couldn't recognize the element from the atom type")
            exit()


    @staticmethod
    def _get_mw(element):
        """Assign mass to the atom based on its atomic element.
        """
        
        masses={'H':1.0079,
                'C':12.0107,
                'N':14.0067,
                'O':15.9994,
                'F':18.9984,
                'S':32.065,
                'Cl':35.453,
                }

        if element in masses:
            return masses[element]
        else:
            print('ERROR: element not known')
    

class Group():

    """The group class that manages and stores information about arbitrary groups of atoms
    """ 

    def __init__(self,name):
        """TODO describe function

        :returns: 

        """

        self._name=name
        self._atoms=list()
        self._natoms=0
        self._com=None


    @property
    def name(self):
        return self._name

    @property
    def atoms(self):
        return self._atoms

    @property
    def natoms(self):
        return len(self._atoms)

    @property
    def com(self):
        if self._com==None:
            self._com=np.array([0.,0.,0.])
            totmass=0
            for atom in self._atoms:
                self._com+=np.array(atom.coordinates)*atom.mw
                totmass+=atom.mw
            self._com/=totmass
        return self._com

    @atoms.setter
    def atoms(self,alist):
        self._atoms=alist

        
    def _check_duplicates(self):

        """Check if there are duplicate atoms in the group and return the sorted list of unique atom objects

        :returns: 

        """
        atom_ids=[]

        for atom in self._atoms:
            atom_ids.append(atom.absindex)

        new_atom_ids, sort_index=np.unique(atom_ids,return_index=True)

        new_atoms_list=[]
        for ival in sort_index:
            new_atoms_list.append(self._atoms[ival])

        self.atoms=new_atoms_list
        

    
        

    
        

