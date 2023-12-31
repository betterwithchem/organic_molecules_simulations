
"""
  Functions to define molecular properties that can be used in
  the calculation of collective variables and analysis of simulations 
"""

def findTorsionalAngles(structure: str, structure_format='pdb'):
    """
       Args:
           structure (str): structure file of the molecule.
           structure_format (str, optional): format of the structure file. Defaults to 'pdb'.
       Returns:
           torsion_angle_indices (list): list of lists identifying the dihedral angles of the molecule. Values are 0-indexed.

      Adapted from code by Edgar Olehnovics
    """

    from rdkit import Chem
    from rdkit.Chem.Lipinski import RotatableBondSmarts
    import numpy as np
    from scipy.spatial.distance import cdist

    accepted_formats=['pdb','mol2']
    
    if structure_format.lower() == 'pdb':    
        #Loading pdb file of molecule as mol
        mol = Chem.MolFromPDBFile(structure)
    elif structure_format.lower() == 'mol2':
        mol = Chem.MolFromMol2File(structure)
    else:
        print("Error: structure formats accepted for the evaluation of torsional angles are {}".format(accepted_formats))
        exit()
       
    # mapping atom indices to atoms
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i)

    #Create a list of all bond atom pairs in the molecule
    bond_atom_pairs = [ [b.GetBeginAtomIdx(),b.GetEndAtomIdx()] for b in mol.GetBonds() ]

    #Identify all bond atom pairs which match with the Rotatable bonds SMARTS string, which identifies
    #the atom bond pairs which are rotatable
    rot_atom_pairs = mol.GetSubstructMatches(RotatableBondSmarts)
    rot_atom_pairs = np.array([list(pair) for pair in rot_atom_pairs])

    n_torsional_angles = len(rot_atom_pairs)

    #identify the atoms with rotatable bonds
    ha_rot = np.unique(rot_atom_pairs.flatten()).tolist()

    #convert RDkit mol object to a MolBlock object (a string representation)
    mb = Chem.MolToMolBlock(mol,confId=-1)

    n_atoms = int(mb[26:29]) # <-- or, from checking the string above is all ok.

    #Get atom coordinates
    coor = np.array([[mol.GetConformer().GetAtomPosition(i).x,
                      mol.GetConformer().GetAtomPosition(i).y,
                      mol.GetConformer().GetAtomPosition(i).z] for i in range(n_atoms)])

    #work out distances between all atoms
    D_coor = cdist(coor, coor, metric='euclidean')
    ' distance < 1.862 (A) between pair of atoms is a bond: '
    #Class atom pairs less than 1.862 A apart as bonded#
    D_coor_b = np.where(D_coor<=1.862,1,0)-np.eye(len(coor))
    node_degree = D_coor_b.sum(0)

    #So far, rotatable bonds between two atoms have been identified. To calculate torsional angles, two additional
    #bonds are needed to either side of the rotatable one. This code here works out the indices for two additional atoms
    #to find those bonds
    torsion_angle_indices = []

    #iterating over torsional angles...
    for i in range(n_torsional_angles):
        #...picking out atom indices from listed atom pairs
        b = rot_atom_pairs[i,0]
        c = rot_atom_pairs[i,1]
        #picking out atoms on either side of rotatable bonded atoms(to facilitate torsional angle identification)
        a_ = list(set(np.where(D_coor_b[b]>0)[0]) - {c})#[0]
        d_ = list(set(np.where(D_coor_b[c]>0)[0]) - {b})#[0]

        a = a_[np.argmax(node_degree[a_])]
        d = d_[np.argmax(node_degree[d_])]

        torsion_angle_indices.append([a,b,c,d])

    torsion_angle_indices = np.array(torsion_angle_indices)

    return torsion_angle_indices
