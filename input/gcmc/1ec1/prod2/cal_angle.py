import sys, os
import mdtraj as md
import numpy as np


ind = [[778,779,782,783]]
ind2 = [[2344,2345,2348,2349]]
t = md.load('1ec1-gcmc-md2.xtc',top='first.gro')
angle = md.compute_dihedrals(t,ind)
angle2 = md.compute_dihedrals(t,ind2)
np.save('torsion3.npy',angle)
np.save('torsion4.npy',angle2)



