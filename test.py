from sim_launch_py.classes import Project

p=Project.new_project(name='test_rot',path='./test_rot',overwrite=True)

p.add_system('s0')
s=p.systems[0]

p.add_molecule(name='sulfamerazine',resname='SUR', structure='Structures/sulfamerazine/sulfamerazine.pdb')
p.add_molecule(name='water',resname='WAT', structure='Structures/water/water.pdb')

s.add_molecule(name='water',knownmolecules=p.molecules)
s.add_molecule(name='sulfamerazine',knownmolecules=p.molecules)

wat=s.molecules[0]
sur=s.molecules[1]

s.box=5

s.createSolventBox(wat,output_structure="{}/box.pdb".format(s.path),density=200)
s.insertSolute(sur,wat,solvent_box="{}/box.pdb".format(s.path),concentration=200,output_structure="{}/start.pdb".format(s.path))




#print(len(s.atoms))
#print([a.atomID for a in s.atoms])
