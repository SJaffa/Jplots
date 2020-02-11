#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 17:09:25 2017

@author: sjaffa
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.ndimage.measurements import center_of_mass as spcom
import momentsmod as mm
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
from astrodendro import Dendrogram


#Read in data file to 3D numpy array
print('....reading data')


data,header= fits.getdata('/home/sarah/Dropbox/PhDwork/170117Moments/SeamusFilaments/Data/PPV_Cubes/BinnedS2.fits',header=True)
#data, header1 = fits.getdata('/home/sarah/Dropbox/PhDwork/170117Moments/SeamusFilaments/Data/PPP_Cubes/S2.fits',header=True)

#data=np.loadtxt('/home/sarah/Dropbox/Herts/GandalfStuff/smallTurb/turbSmall.sf.00201_density_grid.dat',skiprows=19)
#data=data.reshape((512,452,611))


#infile='./Data files and outputs/RUNI032_griddump.out'
#xind,yind,zind,val=np.loadtxt(infile,unpack=True)
#
##Make into 3D array
#xarr=np.unique(xind)
#yarr=np.unique(yind)
#zarr=np.unique(zind)
#
#data=np.zeros((len(xarr), len(yarr), len(zarr)))
#n=0
#for i in range(len(xarr)):
#    for j in range(len(yarr)):
#        for k in range(len(zarr)):
#            data[i,j,k] = val[n]
#            n+=1

#%%
#Compute dendrogram to segment image

print('computing dend')
d = Dendrogram.compute(data)#, 
                       #   min_value=10**0.65,
                       #   min_delta=5,
                       #   min_npix=10)

d.viewer() #Dendrogram viewer (see astrodendro package for detalis)

#%%
#set up 3D plot
fig = plt.figure()
gridspec.GridSpec(5,5)
fsize=8
fig.set_size_inches(w=fsize,h=fsize)
bigax=plt.subplot2grid((5,5), (0,0), colspan=4, rowspan=4, projection='3d')

#Analyse structures (leaves) from dendrogram
print('analysing')

J1s=[]
J2s=[]
J3s=[]
mass=[]
vols=[]
ids=[]
c=0
for shape in d.leaves:
    print(shape)
    ids.append(shape.idx)
    
    #get pixel mask for structure
    mask = np.zeros(data.shape, dtype=bool)
    mask = mask | shape.get_mask()
    gt=np.where(mask,data,0)
    
    #print gt.shape
    xt,yt,zt=shape.indices()
    g=gt[min(xt):max(xt)+1,min(yt):max(yt)+1,min(zt):max(zt)+1]
    #print min(xt),max(xt),min(yt),max(yt),min(zt),max(zt)
    print(g.shape)

    com=spcom(g)
    Vtot=np.count_nonzero(g)
    Mtot=np.sum(g)
    #print Mtot
    mass.append(Mtot)
    vols.append(Vtot)
    
    #Calculate moment of inertia
    yind=np.arange(g.shape[0])
    xind=np.arange(g.shape[1])
    zind=np.arange(g.shape[2])
    
    dx, dy, dz = np.meshgrid(xind,yind,zind)

        
    dx=dx-com[0]
    dy=dy-com[1]
    dz=dz-com[2]
    
    Ixx=np.sum(((dy**2)+(dz**2))*g)
    Iyy=np.sum(((dx**2)+(dz**2))*g)
    Izz=np.sum(((dx**2)+(dy**2))*g)
    
    Ixy=-np.sum(dx*dy*g)
    Iyz=-np.sum(dy*dz*g)
    Ixz=-np.sum(dx*dz*g)
    
    Iyx=Ixy
    Izy=Iyz
    Izx=Ixz
    
    Inertia_tensor=[[Ixx, Ixy, Ixz],
                    [Ixy, Iyy, Iyz],
                    [Ixz, Iyz, Izz]]

    ow,ov = np.linalg.eig(Inertia_tensor)
    eig_ord = np.argsort(ow)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.

    w = ow[eig_ord]
    v = ov[:, eig_ord].T

    I0=(2./5.)*Mtot*(((3.*Vtot)/(4.*np.pi))**(2./3.))
    J1, J2, J3 = (I0-w)/(I0+w)
    
    print(J1, J2, J3)

    J1s.append(J1)
    J2s.append(J2)
    J3s.append(J3)
    c=shape.idx
    
    #Calculate things to plot (J1,2,3 or I1,2,3, or ratios etc.)
    x=J1
    y=J2
    z=J3
    
    plot_text=False
    lim=0.15 #accuracy limit for checking shape type (symmetrical/prolate/oblate/triaxial)
    
    if abs(J1-J2)<lim:
        if abs(J2-J3)<lim:
            bigax.plot([x], [y], [z],marker='s', color='b')
            print('======sym')
            if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='b', fontsize=12)
        else:
            bigax.plot([x], [y], [z],marker='o', color='r')
            print('======oblate')
            if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='r', fontsize=12)
    elif abs(J2-J3)<lim:
        bigax.plot([x], [y], [z],marker='^', color='g')
        print('=======prolate')
        if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='g', fontsize=12)
    else:
        bigax.plot([x], [y], [z],marker='v', color='k')
        print('======tri')
        if plot_text:
                bigax.text(x, y+0.05, z+0.05,"%i"%c,color='k', fontsize=12)
    
    print()
    c+=1
    

#%%
    
J1s=np.array(J1s)
J2s=np.array(J2s)
J3s=np.array(J3s)

#Save structure id, J values, mass nad volume to file
outfile='.'+infile.split('.')[1]+'.jout'
np.savetxt(outfile, 
           np.vstack((ids, J1s, J2s, J3s, mass, vols)).T, 
           header='ID\tJ1\tJ2\tJ3\tMass\tnpix',
           fmt='%i\t%7.5f\t%7.5f\t%7.5f\t%6f\t%i')

#Tidy up 3D plot
bigax.axis('equal')
bigax.set_xlim((-1.1,1.1))
bigax.set_ylim((-1.1,1.1))
bigax.xaxis.tick_top()
bigax.grid()


#Plot axes
bigax.plot([0,0],[0,0],[-1.3,1.3],'k-')
bigax.plot([0,0],[-1.3,1.3],[0,0],'k-')
bigax.plot([-1.3,1.3],[0,0],[0,0],'k-')
#
bigax.set_xlabel(r'$J_{1}$')
bigax.set_ylabel(r'$J_{2}$')
bigax.set_zlabel(r'$J_{3}$')
bigax.xaxis.set_label_position('top') 

plt.savefig('.'+infile.split('.')[1]+'_jplot.pdf')