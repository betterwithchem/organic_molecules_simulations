import os
import sim_launch_py.utilities as util

class Project():

    def __init__(self,name=None, path=None):

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
    def name(self):
        return self._name

    @property
    def systems(self):
        return self._systems
    
    @property
    def molecules(self):
        return self._molecules


    @staticmethod
    def help():
        print("""Help! I need somebody
Help! Not just anybody
Help! You know I need someone
Help!""")


    def add_molecule(self,name=None,resname='UNK', structure=None):

        import shutil
        for mol in self._molecules:
            if name == mol._name:                
                print("Ignoring molecule {} ({}), it already exists a molecule with the same name in the project.".format(name,resname))
                return

        newmolecule=Molecule(name=name,resname=resname,structure=structure)
        self._molecules.append(newmolecule)        
        shutil.copy(structure, self._init_struct_path)

        
    def add_system(self,name=None):

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

    

        
class System():

    def __init__(self,name=None,path=None):

        self._name=name
        self._path=path
        self._molecules=list()
        self._temperature=0
        self._box=list()

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

    @temperature.setter
    def temperature(self,T):
        self._temperature=T

    @box.setter
    def box(self,box):
        if isinstance(box,float) or isinstance(box,int) or len(box)==1:
            self._box=[box, box, box, 90, 90, 90]
        elif len(box)==6:
            self._box=box

    def add_molecule(self,name=None,moltype=None,knownmolecules=None):

        for mol in knownmolecules:
            if name==mol.name:
                newmolecule=mol
                newmolecule._mol_attributes=moltype
                self._molecules.append(newmolecule)
                return

        print("Couldn't add molecule {} to system {}. Molecule unknown.".format(name,self.name))
        exit()


    def createSolventBox(self,solvent,output_structure="solvent_box.pdb",density=1000):

        out_path=self.path+'/'+output_structure
        
        os.system("gmx -nobackup editconf -f {0} -o {1} -box {2} {2} {2} -angles 90 90 90 -c".format(solvent.structure_path,
                                                                                           out_path,
                                                                                           self.box[0]))

        nmols=util.estimate_n_molecules(self.box[0]**3,solvent.mw,density)-1
        
        os.system("gmx -nobackup insert-molecules -f {0} -o {0} -ci {0} -nmol {1} -try 10000".format(out_path,
                                                                        nmols))

    def insertSolute(self,solute,solvent,solvent_box="solvent_box.pdb",concentration=0, output_structure="start.pdb"):

        if concentration == 0:
            # then single molecule
            nmols=1
        elif concentration > 0:
            from math import ceil
            nmols=ceil(util.estimate_n_molecules( self.box[0]**3,solute.mw,concentration ))
            

    

        
    #def get_
    

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



    

    
    
