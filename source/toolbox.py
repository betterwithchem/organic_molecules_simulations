import os
from units import unit
from units.predefined import define_units
import scripts.constants

commentchars=['#','@']

define_units()

class Molecule():

    # this class defines molecule types
    
    def __init__(self):

        name=''
        mw=0
        resname='UNK'
        path=''
        structure=''
        top=''
        nmol=0
        
        self._name=name #molparms.species
        self._mw=mw #molparms.MW
        self._resname=resname #molparms.resname
        self._path=path #molparms.path
        self._structure=structure
        self._top=top
        
        #self._smiles=smiles

        #self.addStructPath(self._path+'/'+self._name+'.pdb')


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,val):
        self._name=val

    @property
    def mw(self):
        return self._mw

    @mw.setter
    def mw(self,val):
        self._mw=val

    @property
    def resname(self):
        return self._resname

    @resname.setter
    def resname(self,val):
        self._resname=val

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self,val):
        self._path=val

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self,val):
        self._structure=val

    @property
    def nmol(self):
        return self._nmol

    @nmol.setter
    def nmol(self,nmol):
        self._nmol=nmol

    def addIncludePath(self,path):

        if os.path.isfile(path):
            self._includepath=path
        else:
            exit('ERROR: parameters file not found in {}'.format(path))

    @property
    def includepath(self):
        return self._includepath

    def addStructPath(self,path):

        # todo: add check on extension of the structure file
        
        if os.path.isfile(path):
            self._structurepath=path
        else:
            exit('ERROR: structure file not found in {}'.format(path))

    @property
    def structurepath(self):
        return self._structurepath


class System():

    # this class define the system to be simulated.
    # It is a collection of molecules plus other attributes


    #def __init__(self,solvent='',solute='',systemid=0,path=''):
    def __init__(self):

        """
        self._solvent=solvent
        self._solute=solute
        self._systemid=systemid
        self._solvent_density=0
        self._solute_conc=0
        self._path=path 
        """
        self._solvent=None
        self._solute=None
        self._systemid=0
        self._solvent_density=0
        self._solute_conc=0
        self._path=None
        self._temperature=0

    @property
    def solvent(self):
        return self._solvent

    @solvent.setter
    def solvent(self,species):
        self._solvent=species

    @property
    def solute(self):
        return self._solute

    @solute.setter
    def solute(self,species):
        self._solute=species

    @property
    def systemid(self):
        return self._systemid

    @systemid.setter
    def systemid(self,sid):
        self._systemid=sid

    @property
    def solvent_density(self):
        return self._solvent_density

    @solvent_density.setter
    def solvent_density(self,density):
        self._solvent_density=density

    @property
    def solute_conc(self):
        return self._solute_conc

    @solute_conc.setter
    def solute_conc(self,conc):
        self._solute_conc=conc

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self,path):
        self._path=path

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self,temperature):
        self._temperature=temperature


    def addBox(self,box_vector=[],shape='cubic'):

        if shape in ['cubic']:
            self._box_shape=shape

        else:
            exit("ERROR: Invalid box shape ({})".shape)

        if shape=='cubic':
            if type(box_vector)==float or type(box_vector)==int:
                box_vector=[box_vector, box_vector, box_vector, 90., 90., 90.]
            elif len(box_vector)==1:
                box_vector=[box_vector[0], box_vector[0], box_vector[0], 90., 90., 90.]
            else:
                narg_shape_error(shape)
                
                

        self._box_vector=box_vector

    @property
    def shape(self):
        return self._box_shape

    @property
    def box_vector(self):
        return self._box_vector

def narg_shape_error(shape):
    exit('ERROR: invalid number of arguments for {0} box'.format(shape))
    
class Project():

    # this is a class where different systems are grouped

    def __init__(self,title='new_project',path='./new_project',overwrite=False,continuation=False):

        # check that needed commands are available
        used_commands=['antechamber','parmchk2','tleap', 'gmx']

        print("Checking that needed commands are available...")
        
        for c in used_commands:
            if isthere(c):
                exit("I cannot find command {}, add it to your PATH before running the job!".format(c))


        path=os.path.abspath(path)

        self._title=title
        self._path=path
        self._species_list=[]

        self.setAtomtypes()
        
        if os.path.isdir(path):
            os.chdir(path)
            if overwrite:
                print("W: Project {0} already exists in directory {1} and WILL BE OVERWRITTEN".format(title,path))
                files=os.listdir(path+'/')
                for fname in files:
                    if os.path.isfile(fname):
                        #print("remove file {0}/{1}".format(path,fname))
                        os.unlink(fname)
            else:
                if not continuation:
                    exit("Project {0} already exists and will NOT be overwritten and will NOT continue. Exiting... ".format(title))
                elif continuation:
                    print("W:  Project {0} already exists in directory {1}. Going on with parameters generation and systems building.".format(title,path))                           
        else:
            os.makedirs(path)        
            os.chdir(path)
        
    @property
    def title(self):
        return self._title

    @property
    def path(self):
        return self._path

    @property
    def species(self):
        return self._species_list

    @species.setter
    def species(self,val):
        self.species.append(val)

    def addSpecies(self,species_df):

        self.molecules=species_df

        #print(species_df.loc[0])
        #print(self.molecules)
             
        for i,s in enumerate(self.molecules.molname):
            
            #print(species_df[species_df.molname==s].path[i],species_df[species_df.molname==s].resname[i])
            
            setattr(self,s,Molecule())
            
            name=self.molecules.loc[i].molname
            mw=self.molecules.loc[i].MW
            path=self.molecules.loc[i].path
            resname=self.molecules.loc[i].resname
            
            getattr(self,s).name=name
            getattr(self,s).mw=mw
            getattr(self,s).path=path
            getattr(self,s).resname=resname
            

            if os.path.isdir(name) is False:
                try:
                    os.makedirs(name)
                except:
                    exit('Could not create directory {0} in {1}. Exiting...'.format(name,os.path.abspath('.')))

            for filetype in ['.pdb']:
                if os.path.isfile(path+'/'+name+filetype) :
                    copyFiles( path+'/'+name+filetype, path=name )

                    newpath=os.path.abspath(os.path.curdir)+'/'+s
                    
                    getattr(self,s).path=newpath
                    getattr(self,s).structure=newpath+'/'+name+filetype
                    
                else:
                    print("couldn't find file {}".format(p+'/'+s+filetype))

            #for filetype in ['.top','.itp']:
            #    if os.path.isfile(p+'/'+resname+filetype) :
            #        copyFiles( p+'/'+resname+filetype, path=s )
            #    #else:
            #    #    print("couldn't find file {}".format(p+'/'+s+filetype))

            self.species=name

    def setAtomtypes(self,atomtypes={}):

        self._atomtypes=atomtypes

    @property
    def atomtypes(self):
        return self._atomtypes
    
    @atomtypes.setter
    def atomtypes(self,atype_dict):
        self._atomtypes=atype_dict
               
    def addSystems(self,dataframe):

        setattr(self,'nsystems',0)
        setattr(self,'systems',[])

        wd=self.path

        for i,d in dataframe.iterrows():

            name='s'+str(i)   # need to find a smarter way to name systems

            setattr(self,name,System())
            

            getattr(self,name).solvent=dataframe.at[i,'solvent']
            getattr(self,name).solvent_dens=dataframe.at[i,'dens_solv']
            getattr(self,name).solute=dataframe.at[i,'solute']
            getattr(self,name).solute_conc=dataframe.at[i,'conc_solute']
            getattr(self,name).systemid=i
            getattr(self,name).temperature=dataframe.at[i,'temperature']

            getattr(self,name).path=wd+'/'+name

            print("Creating system {0} with {1},{2}, and {3},{4} in {5} at {6}K".format(name,
                                                                                        getattr(self,name).solvent,
                                                                                        getattr(self,name).solvent_dens,
                                                                                        getattr(self,name).solute,
                                                                                        getattr(self,name).solute_conc,
                                                                                        getattr(self,name).path,
                                                                                        getattr(self,name).temperature))
            
            os.makedirs(name,exist_ok=True)
            
            self.nsystems+=1
            self.systems.append(name)

    def simulationBox(self):

        for i,s in enumerate(self.systems):
            getattr(self,s).addBox(box_vector=[10, 10, 10, 90, 90, 90])

            print(getattr(self,s).box_vector)

    def createIncludeFiles(self):

        for i,s in enumerate(self.species):

            setattr(getattr(self,s),'itp',getattr(self,s).path+'/'+s+'.itp')
            #print(getattr(self,s).top,getattr(self,s).itp)
            extract_molecule_from_gmx_top(getattr(self,s).top,getattr(self,s).itp)

        
        
        
        

        
    



def copyFiles(*args, path='',verbose=False,overwrite=False):

    if os.path.isdir(path) is not True:
        exit('Error while copying files: You need to specify a reachable path to copy files to.')

    for fname in args:
        if not overwrite:
            try:
                os.system('cp -n {0} {1}/'.format( fname, path ) )
            except:
                print('W: failed to copy file {0} in directory {1}'.format( fname, path) )
            else:
                if verbose:
                    print("copied file {0} in {1}".format( fname, path ))
        elif overwrite:
            try:
                os.system('cp -f {0} {1}/'.format( fname, path ) )
            except:
                print('W: failed to copy file {0} in directory {1}'.format( fname, path) )
            else:
                if verbose:
                    print("copied file {0} in {1}".format( fname, path ))

    

    

     


        
def isthere(command):

    # this function simply checks that a command exists.
    # it can be used to check if external programs called by this library
    # are available at run time
    
    return (os.system("which "+command)!=0)

def conc_from_molvol(molecules,volume):

    # this function computes the concentration of a simulation box

    return molecules/volume

def vol_from_molconc(molecules,concentration,mw):

    # this function computes the volume needed for a simulation of a system
    # at a given concentration with a given number of molecules

    # it assumes:
    # - concentration in g/L
    # - molecular weight (mw) in g/mol
    # - Nav in molec/mol
    # - volume in nm**3

    #print(concentration, scripts.constants.Nav)

    if concentration>0:
        volume = molecules/concentration * mw/6.022 * 1e4
        return volume.astype(float)
    elif concentration==0:
        volume=10**3
        return volume

    

def get_side(volume,shape='cubic'):

    if shape=='cubic':
        side=(volume)**(1/3)
    elif shape=='dodecahedron':
        side=(volume/0.707)**(1/3)

    return side

def write_gmx_topology(proj,path_itp,nmol):

    # this function combines multiple itps (in the order they are given!!)
    # and outputs a full topology file in gromacs format
    
    return 0

def extract_molecule_from_gmx_top(topfile,path_itp_file):

    # this function takes a gromacs topology file of a single molecule type
    # and extracts the following sections:
    # [ moleculetype ], [ atoms ], [ bonds ], [ pairs ], [ angles ], [ dihedrals ])

    print_line=False
    
    with open(topfile,'r') as top:
        with open(path_itp_file,'w') as itp:
            for line in top:
                if line.strip() == "[ moleculetype ]":
                    print_line=True
                if line.strip() == "[ system ]":
                    print_line=False
                if print_line:
                    itp.write(line)

    return

def extract_atomtypes_from_gmx_top(topfile,atomtypes):

    # this function takes a gromacs topology file
    # and extracts the atom types in the [ atomtypes ] section
    
    readline=False
    readnext=False
    with open(topfile,'r') as f:
        for line in f:
            if line[0]!=";":
                if line.strip() == "[ atomtypes ]":
                    readnext=True
                if (not readnext) and (line[0] == "["):
                    readline=False
                if readline:
                    cols=line.split()
                    if len(cols)>0:
                        atomtypes=addatomtype(cols,atomtypes)
                if readnext:
                    readline=True
                    readnext=False

    
    return atomtypes

def addatomtype(values,atomtypes):

    # ; name    at.num    mass    charge ptype  sigma      epsilon
    atname=values[0]
    atnum=values[1]
    mass=values[2]
    charge=values[3]
    ptype=values[4]
    sigma=values[5]
    epsilon=values[6]

    #print(values)
    print(atname,atomtypes)
    
    if atname not in atomtypes:
        atomtypes[atname]={}
        atomtypes[atname]["atnum"]=atnum
        atomtypes[atname]["mass"]=mass
        atomtypes[atname]["charge"]=charge
        atomtypes[atname]["ptype"]=ptype
        atomtypes[atname]["sigma"]=sigma
        atomtypes[atname]["epsilon"]=epsilon



    else: # if already exists, check that these parameters are consistent
    
        ndiff=0
        
        if checkValues(atnum,atomtypes[atname]["atnum"]):
            print("atomic number of {} is different from what already in database".format(atname))
            ndiff+=1
        if checkValues(mass,atomtypes[atname]["mass"]):
            print("mass of {} is different from what already in database".format(atname))
            ndiff+=1
        if checkValues(charge,atomtypes[atname]["charge"]):
            print("charge of {} is different from what already in database".format(atname))
            ndiff+=1
        if checkValues(ptype,atomtypes[atname]["ptype"]):
            print("ptype of {} is different from what already in database".format(atname))
            ndiff+=1
        if checkValues(sigma,atomtypes[atname]["sigma"]):
            print("sigma of {} is different from what already in database".format(atname))
            ndiff+=1
        if checkValues(epsilon,atomtypes[atname]["epsilon"]):
            print("epsilon of {} is different from what already in database".format(atname))
            ndiff+=1

        if ndiff>0:
            exit("ERROR: atomtype {0} already exists and has different parameters from those that you tried to add.\n{0} is \n{1}, \nyou tried to add \n{2}".format(atname,atomtypes[atname],values))

    
    return atomtypes

    

def checkValues(val1,val2):

    return val1!=val2
