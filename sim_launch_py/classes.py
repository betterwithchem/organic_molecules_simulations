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
        
        self._project_path=path
        self._mdp_path=path+'/mdp/'
        self._init_struct_path=path+'/Initial_structures/'
        self._topology_path=path+'/Topologies/'
        self._systems_path=path+'/Systems/'
        self._pickle_path=path+'/.multisim.pkl'
        self._name=name
        self._systems=list()
        self._molecules=list()
        self._job_script_path=None

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

    @job_script_path.setter
    def job_script_path(self,sp):
        self._job_script_path=sp
        
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
    - gromacs : path of the binary for gromacs
    - ambertools : path of the directory where the binaries of the antechamber tools are (namely: antechamber, atomtype, leap...)
    ** these path could be added in an automatic fashion by using shutil.which() on "gmx","gmx_mpi","tleap",... and fail in case they are not found (or give the possibility to add them manually at prompt). I would personally go for the a
utomated fashion (with clear indication in a log file), because, in any case, they need to be in the environment of the system with all the loaded libraries (so, sourcing need to be done BEFORE than running the project). Also, in the cas
es where these programs are NOT needed (e.g. when adding molecules to the project/systems and systems to the project, it is not important that they are defined at runtime).
 
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

    def add_molecule(self,name=None,resname='UNK', structure=None):
        """Add molecule to project

        Args:
            name (str, optional): molecule name. Defaults to None.
            resname (str, optional): residue name. Defaults to 'UNK'.
            structure (str, optional): molecular structure file. Defaults to None.
        """
        import shutil
        for mol in self._molecules:
            if name == mol._name:                
                print("Ignoring molecule {} ({}), it already exists a molecule with the same name in the project.".format(name,resname))
                return

        newmolecule=Molecule(name=name,resname=resname,structure=structure)
        self._molecules.append(newmolecule)        
        shutil.copy(structure, self._init_struct_path)
        

        
    def add_system(self,name=None):
        """Add system to project

        Args:
            name (str, optional): system name. Defaults to None.
        """
        for sys in self._systems:
            if name==sys.name:
                print("System {} already exists! ###### This will need to be changed to an error ##### ".format(name))
                return

        syspath=self.project_path+'/'+name
        if os.path.isdir(syspath) is False:
            os.makedirs(syspath)
        
        newsystem=System(name=name,path=syspath)
        self._systems.append(newsystem)
        


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
        #self._write_output()
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
            print("No PyPol project found in '{}'. Use the 'Project.new_project' module to create a new project."
                  "".format(project_folder))

    def write_sub_command(self,scriptname,system='', template=''):

        import shutil
        recognised_systems=['bash','myriad']

        scriptname=self.project_path+'/'+scriptname

        if template=='':
            template=self.job_script_path+'/myriad.job'
        else:
            template=os.path.abspath(template)

        filename=os.path.basename(template)
        
        if system=='bash':
            for sys in self.systems:
                print("hey!")

        elif system=='myriad':
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
    - new_simulation(self, simtype: str, mdrun_options='', mdp='', print_bash=True, name='',maxwarn=0,start_coord='',
gmxbin=''):
    - print_command(self, bash_file):
   
    """

    def __init__(self,name=None,path=None):
        """System Class Constructor

        Args:
            name (str, optional): Systen name. Defaults to None.
            path (str, optional): System path. Defaults to None.
        """

        self._name=name
        self._path=path
        self._molecules=list()
        self._temperature=0
        self._box=list()
        self._simulations=list()
        self._run_command=None
        

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
    
    @temperature.setter
    def temperature(self,T):
        self._temperature=T

    @box.setter
    def box(self,box):
        if isinstance(box,float) or isinstance(box,int) or len(box)==1:
            self._box=[box, box, box, 90, 90, 90]
        elif len(box)==6:
            self._box=box

    @run_command.setter
    def run_command(self,command):
        self._run_command=command

    
    def add_molecule(self,name=None,moltype=None,knownmolecules=None):

        for mol in knownmolecules:
            if name==mol.name:
                newmolecule=mol
                newmolecule._mol_attributes=moltype
                self._molecules.append(newmolecule)
                return

        print("Error: Couldn't add molecule {} to system {}. Molecule unknown.".format(name,self.name))
        exit()


    def createSolventBox(self,solvent,output_structure="solvent_box.pdb",density=None,nmols=None):

        
        out_path=self.path+'/'+output_structure

        if (density is not None) and (nmols is not None):
            
            volume_box=solvent.mw*10/(density*6.022)*nmol
            self.box[0]=volume_box**(1/3)
            print("Both density and nmol have been defined: changing side of the box to {}.".format(self.box[0]))

        else:
            volume_box=self.box[0]**3

        nmols=util.estimate_n_molecules(volume_box,solvent.mw,density)-1
        
        os.system("gmx -nobackup editconf -f {0} -o {1} -box {2} {2} {2} -angles 90 90 90 -c".format(solvent.structure_path,
                                                                                           out_path,
                                                                                                     self.box[0]*(1.05)))
               
        os.system("gmx -nobackup insert-molecules -f {0} -o {0} -ci {0} -nmol {1} -try 20000".format(out_path,
                                                                        nmols))

        print("Checking number of molecules..")
        inserted_mols=util.check_number_molecules(out_path,solvent)

        if inserted_mols!=(nmols+1):
            print("Error! Inserted number of mol{} molecules required, but after {} trials only {} where inserted!".format(nmols, 20000, inserted_mols))
            exit()

        
    def insertSolute(self,solute,solvent,solvent_box="solvent_box.pdb",concentration=0, output_structure="start.pdb"):

        in_file=os.path.abspath(self.path+'/'+solvent_box)
        out_file=os.path.abspath(self.path+'/'+output_structure)

        if concentration == 0:
            # then single molecule
            nmols=1
        elif concentration > 0:
            from math import ceil
            nmols=util.estimate_n_molecules( self.box[0]**3,solute.mw,concentration )

        #print(self.name,solute.name,solvent.name,concentration,nmols)

        os.system("gmx insert-molecules -f {0} -o {1} -ci {2} -nmol {3} -try 10000 -replace {4} ".format( in_file,
                                                                                                          out_file,
                                                                                                          solute.structure_path,
                                                                                                          nmols,
                                                                                                          solvent.resname))

        print("Checking number of molecules..")
        inserted_mols=util.check_number_molecules(out_file,solute)

        if inserted_mols!=nmols:
            print("Error! Inserted number of mol{} molecules required, but after {} trials only {} where inserted!".format(nmols, 10000, inserted_mols))
            exit()

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
                    

    def new_simulation(self, simtype: str, mdrun_options='', mdp='', print_bash=True, name='',maxwarn=0,start_coord='',gmxbin=''):
       
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

    
    def print_command(self, bash_file):

        bash_file=self.path+'/'+bash_file

        with open(bash_file,'w') as f:
            f.write("#!/bin/bash\n\n")
            for s in self.simulations:
                f.write("# {}\n".format(s.name))
                f.write(s.bash_command+'\n\n')


                
            
        
            

    

        

    

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



    

    
    
