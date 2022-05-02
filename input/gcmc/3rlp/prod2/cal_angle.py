import sys, os
import mdtraj as md
import numpy as np

ind = [[680,681,684,685]]
t = md.load('3rlp-gcmc-md2.xtc',top='first.gro')
angle = md.compute_dihedrals(t,ind)
np.save('torsion.npy',angle)
