import sim_launch_py.classes as mp
import sys
import numpy as np

from sim_launch_py.utilities import get_distance

pp=mp.Project.new_project('load_pdb_proj','./load_pdb_proj',overwrite=True)

s0=pp.add_system('s0')

s0.add_molecule('olanzapine',pdbfile='./UNOGIN_eq406_bulk.pdb')

print('added...')

comfile=open('com.pdb','w')

delete=[]

#print(box_center)
#print(s0.box[0:3])

cutoff=1.5  # nm

com=[]

print(s0.box)

for im,m in enumerate(s0.molecules):
    
    com.append(m.get_molecule_com(s0.box))

refmol=int(len(s0.molecules)/2)

deg2rad=np.pi/180
n2=(np.cos(s0.box[3]*deg2rad)-np.cos(s0.box[5]*deg2rad)*np.cos(s0.box[4]*deg2rad))/np.sin(s0.box[5]*deg2rad)
M=M=np.array([[s0.box[0], 0, 0],
              [s0.box[0]*np.cos(s0.box[5]*deg2rad), s0.box[1]*np.sin(s0.box[5]*deg2rad), 0],
              [s0.box[0]*np.cos(s0.box[4]*deg2rad), s0.box[1]*n2, s0.box[2]*np.sqrt(np.sin(s0.box[4]*deg2rad)**2-n2*n2)]])

refcoords=M.dot(np.array([0.7,0.5,0.3]))

print('ref coords mol {}: {}'.format(refmol,refcoords))
    
for im,m in enumerate(s0.molecules):
        
    comfile.write("ATOM  {:5d} {:4s} UNK     1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:1s}\n".format(im,'COM',com[im][0]*10,com[im][1]*10,com[im][2]*10,1.00,3.95,'X'))

    #print('....{}...'.format(im))
    dist=get_distance(com[im],ref=refcoords,box=s0.box)
   
    if dist>cutoff:
        #print(im, com[im],com[257], dist, cutoff)
        delete.append(m.index)
    #else:
    #    print(im,dist,cutoff, com[im],com[257])
    
print(len(s0.molecules))
#print(delete)
s0.delete_molecules(delete)
print(len(s0.molecules))

fout=open('filtered_molecules.pdb','w')

fout.write("CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}{:10s}".format(s0.box[0]*10,
                                                                           s0.box[1]*10,
                                                                           s0.box[2]*10,
                                                                           s0.box[3],
                                                                           s0.box[4],
                                                                           s0.box[5],
                                                                           'P-1'))
for im,m in enumerate(s0.molecules):
    for ia,a in enumerate(m.atoms):
           aa=len(m.atoms)*im+ia
           fout.write("ATOM  {:5d} {:4s} UNK     1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:1s}\n".format(aa,a.name,a.coordinates[0]*10,a.coordinates[1]*10,a.coordinates[2]*10,1.00,3.95,a.element))
        
