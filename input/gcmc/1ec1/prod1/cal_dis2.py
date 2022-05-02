import sys, os
import mdtraj as md
import numpy as np

#dis0 = []
#dis1 = []

t = md.load('1ec1-gcmc-md.xtc',top='first.gro')
ind = [[37101,784],[37101,2350]]
dis = md.compute_distances(t,ind)
dis0 = (dis[:,0])
dis1 = (dis[:,1])
np.save('dis3.npy',dis0)
np.save('dis4.npy',dis1)

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.plot(range(len(dis0)),dis0,'o-')
#plt.plot([0,len(dis0)],[0.30943495,0.30943495])
#plt.plot([0,len(dis0)],[0.975,0.975])
plt.plot([0,len(dis0)],[0.55714,0.55714])
plt.ylim(0.2,0.8)
plt.xlabel('frame #')
plt.ylabel('distance (nm)')
plt.savefig('dis3.png')
plt.close()


plt.plot(range(len(dis1)),dis1,'o-')
#plt.plot([0,len(dis1)],[0.29198116,0.29198116])
#plt.plot([0,len(dis1)],[0.975,0.975])
plt.plot([0,len(dis1)],[0.59268,0.59268])
plt.ylim(0.2,0.8)
plt.xlabel('frame #')
plt.ylabel('distance (nm)')
plt.savefig('dis4.png')
plt.close()

