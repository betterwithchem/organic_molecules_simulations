class Molecule():
    """The molecule class that stores and manage all the information and methods.

    Attributes:
       name : name of the molecule type
       resname : residue name
       structure_path : location of the structure of the molecule type
       topology_path : location of the topology (.top) file of the molecule type
       include_path : location of the include topology (.itp) file of the molecule type 
       mw : molecular weight in [g/mol]
       mol_attributes : arbitrary attributes of the molecule
       nmols : number of molecules of this type
       natoms : number of atoms per molecule 
       atoms : list of atom types of the molecule type

    Methods:
       help(): print the help for this class.
    """

    def __init__(self,name: str, index=None, resname='UNK', structure=None):
        """Molecule Class Constructor
  
        Args: 
           name (str) : name of the molecule
           resname (str, optional) : residue name for the molecule. Defaults to 'UNK'.
           structure_path (str, optional) : position of the structure file. Defaults to None.
        """

        self._name=name
        self._resname=resname
        self._structure_path=structure
        self._topology_path=None
        self._include_path=None
        self._mw=None
        self._mol_attributes=[]
        #self._nmols=0
        self._natoms=0
        self._atoms=list()
        self._index=index
        self._contact_matrix=None

        #to be considered
        #self._rotatable=list()
        #self._com=list()

        import os
        
        if structure is not None:

            # if a template structure is provided, get
            # - atoms
            # - contact matrix
            # - mass
            # - path of the structure (copy the structure locally to the project and use that as template)

            # input structure can be in any format accepted by antechamber.
            # antechamber is used here to convert the input file to mol2 (if not already in mol2)
            # in order to get the bonds
            
            accepted_formats=["ac", "mol2", "pdb", "mpdb", "prepi", "prepc", "gzmat", "gcrt", "mopint",
                            "mopcrt", "gout", "mopout", "alc", "csd", "mdl", "hin", "rst"]

            extension=structure.split('.')[-1]
            
            if extension not in accepted_formats:
                print("ERROR: provided structure is not in a format compatible with Antechamber. Provided file ({}) is a {} file. Accepted formats are {}.".format(structure,extension,accepted_formats))
                exit()

            if extension!="mol2":

                mol2file=os.path.splitext(structure)+'.mol2'
                os.system("{0} -i {1} -fi {2} -o {3} -fo mol2 -j 5 -at sybyl -dr no".format('antechamber',
                                                                                            structure,
                                                                                            extension,
                                                                                            mol2file))
            else:
                mol2file=structure


            
            

                
        
    
