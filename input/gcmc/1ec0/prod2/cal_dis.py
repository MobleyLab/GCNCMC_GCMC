import sys, os
import mdtraj as md
import numpy as np

#dis0 = []
#dis1 = []

t = md.load('1ec0-gcmc-md2.xtc',top='first.gro')
ind = [[443,779],[2007,2343],[780,2349]]
dis = md.compute_distances(t,ind)
dis0 = (dis[:,0])
dis1 = (dis[:,1])
dis2 = dis[:,2]
np.save('dis.npy',dis)
np.save('dis0.npy',dis0)
np.save('dis1.npy',dis1)
np.save('dis2.npy',dis2)

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.plot(range(len(dis0)),dis0,'o-')
#plt.plot([0,len(dis0)],[0.30943495,0.30943495])
plt.plot([0,len(dis0)],[0.975,0.975])
plt.ylim(0.75,1.15)
plt.xlabel('frame #')
plt.ylabel('distance (nm)')
plt.savefig('dis0.png')
plt.close()


plt.plot(range(len(dis1)),dis1,'o-')
#plt.plot([0,len(dis1)],[0.29198116,0.29198116])
plt.plot([0,len(dis1)],[0.975,0.975])
plt.ylim(0.75,1.15)
plt.xlabel('frame #')
plt.ylabel('distance (nm)')
plt.savefig('dis1.png')
plt.close()

plt.plot(range(len(dis2)),dis2,'o-')
#plt.plot([0,len(dis1)],[0.29198116,0.29198116])
plt.plot([0,len(dis2)],[0.557,0.557])
plt.ylim(0.3,0.7)
plt.xlabel('frame #')
plt.ylabel('distance (nm)')
plt.savefig('dis2.png')
plt.close()

