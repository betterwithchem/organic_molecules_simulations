import sys
import os
import os.path
import parmed as pmd
import shutil
import sim_launch_py.utilities as util

'''
def gaff(mol_path: str, output_path: str, molecule_name: str='UNK',  generate_charges: str='bcc', atomtype: str='gaff2', overwrite: bool=False):

    """This function is used to compute gaff parameters for small molecules using tleap, and to transform them in a gromacs friendly form. 

    :param mol_path: Structure file of the molecule.
    :type mol_path: str
    :param path_output: Output topology (.top) file. 
    :type path_output: str
    :param res_name: Name of the molecule, it will be used in the .top file. Defaults to UNK.
    :type res_name: str, optional
    :param generate_charges: Charge assignment method. Defaults to 'bcc'.
    :type generate_charges: str, optional
    :param atomtype: Atomtype scheme definition (gaff, gaff2). Defaults to 'gaff2'.
    :type atomtype: str, optional
    :param overwrite: Overwrite existing files? Defaults to False.
    :type overwrite: bool, optional
    :returns: 

    """

    
    path_parmchk =  "parmchk2"
    path_antechamber = "antechamber"
    path_leap = "tleap"



'''


def gaff(molecule: str, path_output: str, res_name: str='UNK', generate_charges: str='bcc', atomtype: str='gaff2', overwrite: bool=False):
        
    path_parmchk =  "parmchk2"
    path_antechamber = "antechamber"
    path_leap = "tleap"
    
    # for the moment let's assume that the extension of the coordinates is pdb
    available_extensions = ["ac", "mol2", "pdb", "mpdb", "prepi", "prepc", "gzmat", "gcrt", "mopint",
                            "mopcrt", "gout", "mopout", "alc", "csd", "mdl", "hin", "rst"]

    available_charge_methods = ["resp", "cm2", "mul", "rc", "bcc", "esp", "gas", "wc"]

    available_atomtypes = ['gaff', 'gaff2']

    molecule_name=molecule.name
    path_coord=molecule.structure_path

    
    # sanity checks and errors:
    project_dir=os.getcwd()
    
    if not os.path.isdir(path_output):
        exit("ERROR: destination directory {} does not exist.".format(path_output))
    if atomtype not in available_atomtypes:
        exit("ERROR: you selected an invalid atomtype '{0}'. Possible atom types are {1}".format(atomtype,available_atomtypes))
    if generate_charges not in available_charge_methods:
        exit("ERROR: you selected an invalid charge method '{0}'. Possible charge methods are {1}".format(generate_charges,available_charge_methods))

    input_file_type = None
    for ext in available_extensions:
        if path_coord.endswith(ext):
            input_file_type = ext
            break

    if not input_file_type:
        print("Error: File '{}' does not have one of the required extension!      \n"
              "List of the File Formats:                                          \n"
              "     file format type  abbre. index | file format type abbre. index\n"
              "     --------------------------------------------------------------\n"
              "     Antechamber        ac       1  | Sybyl Mol2         mol2    2 \n"
              "     PDB                pdb      3  | Modified PDB       mpdb    4 \n"
              "     AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 \n"
              "     Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 \n"
              "     Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 \n"
              "     Gaussian Output    gout    11  | Mopac Output       mopout 12 \n"
              "     Alchemy            alc     13  | CSD                csd    14 \n"
              "     MDL                mdl     15  | Hyper              hin    16 \n"
              "     AMBER Restart      rst     17 ".format(os.path.basename(path_coord)))
        exit()



    if not overwrite:
        if os.path.isfile(path_output+'/{}.top'.format(molecule_name)):
            print("Parameters file for {0} already exists ({0}.top), skipping...".format(molecule_name))
            return


    # change working directory
    tempdir=path_output + '/.temp'
    if os.path.isdir(tempdir) is False:
        os.makedirs(tempdir)
        
    os.chdir( tempdir )
    
    # create mol2 file from input coordinates
    os.system(
        path_antechamber + " -i {0} -fi {4} -o {1}/{2}.mol2 -fo mol2 -c {3} -rn {2} -pf y -at {5} -nc 0 ".format(
            path_coord,
            '.',
            res_name,
            generate_charges,
            ext,
            atomtype))
    
    # create compound library
    with open("lib.leap",'w') as f:
        f.write("{0}=loadmol2 {0}.mol2\n"
                "saveoff {0} {0}.lib\n"
                "savepdb {0} {0}_leap.pdb\n"
                "quit".format(res_name))


    pmd.load_file('{}.mol2'.format(res_name)).fix_charges(precision=4).save('{}.mol2'.format(res_name))  

        
    os.system(
        path_leap + " -f lib.leap"
        )

    os.system(
        path_parmchk + " -i {0}.mol2 -f mol2 -o {0}.frcmod".format(res_name)
    )

    # Check output (atom types such as SO have no parameters yet)
    skip_structure = False
    if os.path.exists("{}.frcmod".format(res_name)):
        file_frcmod = open("{}.frcmod".format(res_name), "r")
        for line in file_frcmod:
            if "ATTN, need revision" in line:
                skip_structure = True
                break
        file_frcmod.close()
    else:
        skip_structure = True

    if skip_structure:
        print("Error: Some of the atoms involved have not been parametrized yet.")
        return

    # save parameters in amber format
    with open("parm.leap","w") as f:
        f.write("source leaprc.{1}\n"
                "loadamberparams {0}.frcmod\n"
                "loadoff {0}.lib\n"
                "a=loadpdb {0}_leap.pdb\n"
                "saveamberparm a {2}.prmtop {2}.inpcrd\n"
                "quit".format(res_name,
                              atomtype,
                              molecule_name)
                )
        
    os.system(path_leap + " -f parm.leap")

    amb2gmx(molecule_name)

    shutil.copy("{}.top".format(molecule_name), path_output)

    os.chdir(path_output)
    extract_molecule_from_gmx_top("{}/{}.top".format(path_output,molecule_name),
                                  "{}/{}.itp".format(path_output,molecule_name))

    extract_atomtypes_from_gmx_top("{}/{}.top".format(path_output,molecule_name),
                                   "{}/atomtypes.itp".format(path_output))

    molecule.topology_path="{}/{}.top".format(path_output,molecule_name)
    molecule.include_path="{}/{}.itp".format(path_output,molecule_name)
    
    
    molecule.natoms=util.numberOfAtomsFromTop(molecule.topology_path)


    os.chdir(project_dir)
    
    return 0

def getTop(molecule: str,fromPath: str='',toPath: str=''):
    """obtain topology file

    :param molecule: molecule name
    :type molecule: str
    :param fromPath: path where the topology is. Defaults to ''.
    :type fromPath: str, optional
    :param toPath: path where the topology goes. Defaults to ''.
    :type toPath: str, optional

    """
    
    fromPath=os.path.abspath(fromPath)
    toPath=os.path.abspath(toPath)
    
    if os.path.isdir(toPath) is False:
        os.makedirs(toPath)

    if os.path.isfile("{0}/{1.name}.top".format(fromPath,molecule)):
        shutil.copy(fromPath+"/"+molecule.name+".top",toPath)
        molecule.topology_path=toPath+"/"+molecule.name+".top"
    else:
        print("top file for {0.name} not found".format(molecule))
        print("looked into {}/{}.top".format(fromPath,molecule.name))
        exit()

    extract_molecule_from_gmx_top(toPath+"/"+molecule.name+".top",
                                  toPath+"/"+molecule.name+".itp")
    molecule.include_path=toPath+"/"+molecule.name+".itp"

    extract_atomtypes_from_gmx_top(toPath+"/"+molecule.name+".top",
                                   toPath+"/atomtypes.itp")

    
    

def amb2gmx(molecule_name: str):
    """Run amb2gmx for a molecule

    :param molecule_name: molecule name
    :type molecule_name: str

    """
    
    amber = pmd.load_file('{}.prmtop'.format(molecule_name), '{}.inpcrd'.format(molecule_name))
    amber.save('{}.top'.format(molecule_name),overwrite=True)
    amber.save('{}.gro'.format(molecule_name),overwrite=True)
    
def extract_molecule_from_gmx_top(topfile: str):
    """This function takes a gromacs topology file of a single molecule type and extracts the following sections:     [ moleculetype ], [ atoms ], [ bonds ], [ pairs ], [ angles ], [ dihedrals ])

    :param topfile: path of the gromacs topology (.top) file
    :type topfile: str

    """
    print_line=False

    itp=''

    with open(topfile,'r') as top:
        for line in top:
            if line.strip() == "[ moleculetype ]":
                print_line=True
            if line.strip() == "[ system ]":
                print_line=False
            if print_line:
                itp+=line

    return itp


def extract_atomtypes_from_gmx_top(topfile: str):
    """This function takes a gromacs topology file and extracts the atom types in the [ atomtypes ] section

    :param topfile: name of the gromacs topology (.top) file
    :type topfile: str

    """
    existing_atypes=[]

    #if os.path.isfile(atomtypefile) is False:
    #    os.system("touch {}".format(atomtypefile))

    atomtypes={}
        
    readline=False
    readnext=False
    with open(topfile,'r') as tf:
        for line in tf:
            if line[0]!=";":
                if line.strip() == "[ atomtypes ]":
                    readnext=True
                if (not readnext) and (line[0] == "["):
                    readline=False
                if readline:
                    cols=line.split()
                    if len(cols)>0:
                        atomtypes.update({cols[0]:{}})
                        values=''
                        for icol in range(1,len(cols)):
                            values+=cols[icol]+"    "
                        atomtypes[cols[0]]['parms']=values
                        #atomtypes.append(atype)

                if readnext:
                    readline=True
                    readnext=False

    return atomtypes

