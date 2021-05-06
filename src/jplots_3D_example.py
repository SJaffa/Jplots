"""
Example script for a simple J3D analysis

"""
import numpy as np
import matplotlib.pyplot as plt
import momentsmod as mm
from astrodendro import Dendrogram


#Read in data and log
print("=============Read in data")

infile='../data/RUNI032_griddump.out'
xind,yind,zind,val=np.loadtxt(infile,unpack=True)

#Make into 3D array
xarr=np.unique(xind)
yarr=np.unique(yind)
zarr=np.unique(zind)

data=np.zeros((len(xarr), len(yarr), len(zarr)))
n=0
for i in range(len(xarr)):
    for j in range(len(yarr)):
        for k in range(len(zarr)):
            data[i,j,k] = val[n]
            n+=1

data[data!=0]=np.log(data[data!=0])
#%%

#Extract the structures of interest
print('============Computing dendrogram')
d = Dendrogram.compute(data, 
                       min_value=2.0,
                       min_npix=9,
                       min_delta=0.1)
d.viewer()

#%%

# Find J values of the leaves (top-level structures)
print('==============Calculating J values')

ids=[]
J1s=[]
J2s=[]
J3s=[]

bigax=mm.plot3dax()
plot_text=False
lim=0.15 #accuracy limit for checking shape type (symmetrical/prolate/oblate/triaxial)

for s in d.leaves:
    ids.append(s.idx)
    g=np.zeros_like(data)
    g[s.get_mask()]=data[s.get_mask()]

    J1, J2, J3 = mm.moments_3d(g)
    J1s.append(J1)
    J2s.append(J2)
    J3s.append(J3)
    
    x=J1
    y=J2
    z=J3
    c=s.idx
    
    if abs(J1-J2)<lim:
        if abs(J2-J3)<lim:
            bigax.plot([x], [y], [z],marker='s', color='b')
            #print '======sym'
            if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='b', fontsize=12)
        else:
            bigax.plot([x], [y], [z],marker='o', color='r')
            #print '======oblate'
            if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='r', fontsize=12)
    elif abs(J2-J3)<lim:
        bigax.plot([x], [y], [z],marker='^', color='g')
        #print '=======prolate'
        if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='g', fontsize=12)
    else:
        bigax.plot([x], [y], [z],marker='v', color='k')
        #print '======tri'
        if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='k', fontsize=12)
    
#%%
   