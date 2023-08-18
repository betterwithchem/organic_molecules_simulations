from typing import Union
import numpy as np
import copy
import os

# General
def create(path, arg_type, backup=True):
    """
    Generate a new directory or a new file.
    TODO Rewrite in a more compact way. Use os.remove instead of os.system(rm ...)
    :param path: Path of the directory/file to generate
    :param arg_type: Is it a file or directory?
    :param backup: If the directory/file already exists, create a backup directory/file
    :return:
    """
    import os
    if backup:
        if arg_type == 'dir':
            path = os.path.relpath(path)
            if not os.path.exists(path):
                os.makedirs(path)
            else:
                print("Directory '{}' already exists!".format(path))
                i = 0
                while True:
                    if os.path.exists(path + ".bkp." + str(i)):
                        i += 1
                        continue
                    else:
                        os.system("mv {0} {0}.bkp.{1}".format(path, i))
                        os.makedirs(path)
                        break
        elif arg_type == 'file':
            if not os.path.exists(path):
                os.mknod(path)
            else:
                print("File '{}' already exists!".format(path))
                i = 0
                while True:
                    if os.path.exists(path + ".bkp." + str(i)):
                        i += 1
                        continue
                    else:
                        os.system("mv {0} {0}.bkp.{1}".format(path, str(i)))
                        os.mknod(path)
                        break
        else:
            print("Only 'dir' and 'file' are available as arg_type")
    else:
        if arg_type == 'dir':
            if not os.path.exists(path):
                os.makedirs(path)
            else:
                os.system("rm -r " + path)
                os.makedirs(path)
        elif arg_type == 'file':
            if not os.path.exists(path):
                os.mknod(path)
            else:
                os.system("rm " + path)
                os.makedirs(path)
        else:
            print("Only 'dir' and 'file' are available as arg_type")

"""
def molecularWeightFromTop(top_path):

    read=False
    mass=0
    
    with open(os.path.abspath(top_path),'r') as f:
        for line in f:
            if line[0]!=';':
                cols=line.split()
                if len(cols)<1:
                    read=False                    
                elif cols[0]==r"[":
                    if cols[1]=="atoms":
                        read=True
                elif read:
                    mass+=float(cols[7])

    return mass

def numberOfAtomsFromTop(top_path):

    read=False
    natoms=0
    
    with open(os.path.abspath(top_path),'r') as f:
        for line in f:
            if line[0]!=';':
                cols=line.split()
                if len(cols)<1:
                    read=False                    
                elif cols[0]==r"[":
                    if cols[1]=="atoms":
                        read=True
                elif read:
                    natoms+=1

    return natoms
"""                    

def estimateNumMolecules(box_volume,mass,density):
    """ Estimate the number of molecules in a box at a given density

    Args:
       box_volume (float) : volume of the simulation box in [nm**3]
       mass (float) : mass of the molecule type in [g/mol]
       density (float) : density of molecules in [g/L]

    Return:
       n (int) : number of molecules
    """
    
    n = density/mass * 6.022 / 10 * box_volume

    return int(n)


def countMolecules(structure_file,molecule):
    """ Count number of molecules of a given type in a structure file

    Args:
       structure_file (str) : name of the structure file
       molecule (Molecule() object) : molecule

    Return:
       nmols (int) : number of molecules
    """

    ext=os.path.splitext(structure_file)
    
    if ext[1]=='.pdb':
        natoms=0
        with open(structure_file,'r') as f:
            for line in f:
                cols=line.split()
                if line[0:4]=='ATOM' or line[0:6]=='HETATM':
                    if molecule.resname in line:
                        natoms+=1
            nmols=natoms/molecule.natoms
            return int(nmols)
    else:
        print("for the moment only pdb files are supported")
        exit()
            

