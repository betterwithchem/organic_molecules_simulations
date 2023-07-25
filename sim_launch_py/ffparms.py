import sys
import os
import os.path
import parmed as pmd
import shutil
import sim_launch_py.utilities as util

def gaff(molecule, path_output, res_name='UNK', generate_charges='bcc', atomtype='gaff2', overwrite=False):

    """
    This function is used to compute gaff parameters for small molecules.
    Options are:
    - res_name: name of the residue. For the moment it is used 
    """
        
    #path_ambertools='/home/matteo/Source/amber22/bin/'
    #path_parmchk = path_ambertools + "/parmchk2"
    #path_antechamber = path_ambertools + "/antechamber"
    #path_leap = path_ambertools + "/tleap"

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
    extract_molecule_from_gmx_top(path_output + molecule_name + ".top",
                                  path_output + molecule_name + ".itp")

    extract_atomtypes_from_gmx_top(path_output + molecule_name + ".top",
                                   path_output + "atomtypes.itp")

    molecule.topology_path=path_output + molecule_name + ".top"
    molecule.include_path=path_output + molecule_name + ".itp"
    
    molecule.mw=util.molecularWeightFromTop(molecule.topology_path)
    molecule.natoms=util.numberOfAtomsFromTop(molecule.topology_path)


    os.chdir(project_dir)
    
    return 0

def amb2gmx(molecule_name):
    amber = pmd.load_file('{}.prmtop'.format(molecule_name), '{}.inpcrd'.format(molecule_name))
    amber.save('{}.top'.format(molecule_name),overwrite=True)
    amber.save('{}.gro'.format(molecule_name),overwrite=True)



    
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

def extract_atomtypes_from_gmx_top(topfile,atomtypefile):

    # this function takes a gromacs topology file
    # and extracts the atom types in the [ atomtypes ] section

    existing_atypes=[]

    if os.path.isfile(atomtypefile) is False:
        os.system("touch {}".format(atomtypefile))
        
    with open(atomtypefile,'r') as f:
        for line in f:
            cols=line.split()
            existing_atypes.append(cols[0])
    
    readline=False
    readnext=False
    with open(atomtypefile,'a') as af:
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
                            atype=cols[0]
                            if atype in existing_atypes:
                                print("Skipping atom type {}: it already exists in {}".format(atype,atomtypefile))
                            else:
                                af.write(line)


                    if readnext:
                        readline=True
                        readnext=False

    

