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

#3D
testlist=['sphere',
          'flat-disk','thick-disk',
          'shell-thin','shell-thick',#'shell-half', 
          #'shell-thin-noise','shell-thick-noise',
          #'ring-thin', 'ring-half',
          'plane-thick','plane-thin','rect-thin',
          'cc',
          'fil-thick','fil-thin','fil-asym',
          'ellipsoid-pro', 'ellipsoid-ob',]#'ellipsoid-ob2','ellipsoid-ob4' ]#'ellipsoid-tri']
        
fig = plt.figure()
gridspec.GridSpec(5,5)
fsize=8
fig.set_size_inches(w=fsize,h=fsize)

bigax=plt.subplot2grid((5,5), (0,0), colspan=4, rowspan=4, projection='3d')
c=0    
loclist=[(0,4),(1,4),(2,4),(3,4),
         (4,0),(4,1),(4,2),(4,3),(4,4)]
"""
"""
print('reading data')


#data,header= fits.getdata('/home/sarah/Dropbox/PhDwork/170117Moments/SeamusFilaments/Data/PPV_Cubes/BinnedS2.fits',header=True)
#data, header1 = fits.getdata('/home/sarah/Dropbox/PhDwork/170117Moments/SeamusFilaments/Data/PPP_Cubes/S2.fits',header=True)

#data=np.loadtxt('/home/sarah/Dropbox/Herts/GandalfStuff/smallTurb/turbSmall.sf.00201_density_grid.dat',skiprows=19)
#data=data.reshape((512,452,611))
xind,yind,zind,val=np.loadtxt('/home/sarah/Dropbox/Herts/GandalfStuff/RUNI032_griddump.out',unpack=True)

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

#%%
print('computing dend')
d = Dendrogram.compute(data, 
                          min_value=10**0.65,
                          min_delta=5,
                          min_npix=10)

d.viewer()

#%%
print('analysing')

J1s=[]
J2s=[]
J3s=[]
mass=[]
vols=[]
ids=[]
c=0
for shape in d.leaves:#testlist:#
    print shape,c
    ids.append(shape.idx)
    
    mask = np.zeros(data.shape, dtype=bool)

    mask = mask | shape.get_mask()
    gt=np.where(mask,data,0)
    #print gt.shape
    xt,yt,zt=shape.indices()
    g=gt[min(xt):max(xt)+1,min(yt):max(yt)+1,min(zt):max(zt)+1]
    #print min(xt),max(xt),min(yt),max(yt),min(zt),max(zt)
    print g.shape
    """
    g=mm.testdata3d(20,shape=shape)
    """
    #scale to have always 500 mass
    #g*=500./np.sum(g)
    
    #g[0,0]=0
    #smallax=plt.subplot2grid((5,5), loclist[c],projection='3d')
    #mm.small_3d_plot(smallax,g)
#    #plt.figure()
#    #smallax=plt.gca()
#    im=smallax.imshow(np.sum(np.sqrt(g.T), axis=0),interpolation='None',
#                      origin='lower',cmap='Greys',vmin=0)
#    plt.setp( smallax.get_xticklabels(), visible=False)
#    plt.setp( smallax.get_yticklabels(), visible=False)
#    smallax.text(5,85,"%i"%c,color='k', fontsize=14)
    #plt.colorbar(im)
#    

    com=mm.find_com_3d(g)
    #smallax.plot([com[0]],[com[1]],'y^')

    com=spcom(g)
    Vtot=np.count_nonzero(g)
    Mtot=np.sum(g)
    #print Mtot
    mass.append(Mtot)
    vols.append(Vtot)
    
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
    #print
    #print(Inertia_tensor)
    ow,ov = np.linalg.eig(Inertia_tensor)
    eig_ord = np.argsort(ow)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.

    w = ow[eig_ord]
    v = ov[:, eig_ord].T
    #print w,v

    #print()
    #print("Inertia tensor z, y, x eigenvalues:")
    #print()
    #print(w)
    #print()
    #print("Inertia tensor z, y, x eigenvectors:")
    #print()
    #print(v)
    ord_vals=w
    ord_vects=v

    normalizer = np.sqrt(max(ord_vals) / Mtot)
    norm_vals = np.sqrt(ord_vals / Mtot)
    
    start = com
    ends = [[norm_vals[0] * ord_vects[0][0], norm_vals[0] * ord_vects[0][1], norm_vals[0] * ord_vects[0][2]],
               [norm_vals[1] * ord_vects[1][0], norm_vals[1] * ord_vects[1][1], norm_vals[1] * ord_vects[1][2]],
               [norm_vals[2] * ord_vects[2][0], norm_vals[2] * ord_vects[2][1], norm_vals[2] * ord_vects[2][2]]]
    #ends = [[(norm_vals[0] - 1) * ord_vects[0][0], (norm_vals[0] - 1) * ord_vects[0][1], (norm_vals[0] - 1) * ord_vects[0][2]],
    #      [(norm_vals[1] - 1) * ord_vects[1][0], (norm_vals[1] - 1) * ord_vects[1][1], (norm_vals[1] - 1) * ord_vects[1][2]],
    #      [(norm_vals[2] - 1) * ord_vects[2][0], (norm_vals[2] - 1) * ord_vects[2][1], (norm_vals[2] - 1) * ord_vects[2][2]]]
    #for e in ends:
    #    smallax.plot([start[0],e[0]],[start[1],e[1]],[start[2],e[2]],'g-')
    
    #I1,I2,t1=mm.i12(g,com)
 
    #this is the 3d way
    I0=(2./5.)*Mtot*(((3.*Vtot)/(4.*np.pi))**(2./3.))
    J1, J2, J3 = (I0-w)/(I0+w)
    
    #this is me playing with 2d way
    #Atot=np.array([np.sum(np.maximum.reduce(g,axis=i)) for i in [0,1,2]])
    #ma=Atot*Mtot
    #J1,J2,J3=(ma - 4*np.pi*w)/(ma + 4*np.pi*w)

    print (J1, J2, J3)

    J1s.append(J1)
    J2s.append(J2)
    J3s.append(J3)
    c=shape.idx
    lim=0.15
    if abs(J1-J2)<lim:
        if abs(J2-J3)<lim:
            #bigax.plot([J1], [J1/J2], [J1/J3],marker='s', color='b')
            print '======sym'
            #bigax.text(J1, (J1/J2)+0.05, (J1/J3)+0.05,"%i"%c,color='b', fontsize=12)
        else:
            #bigax.plot([J1], [J1/J2], [J1/J3],marker='o', color='r')
            print '======oblate'
            #bigax.text(J1, (J1/J2)+0.05, (J1/J3)+0.05,"%i"%c,color='r', fontsize=12)
    elif abs(J2-J3)<lim:
        #bigax.plot([J1], [J1/J2], [J1/J3],marker='*', color='g')
        print '=======prolate'
        #bigax.text(J1, (J1/J2)+0.05, (J1/J3)+0.05,"%i"%c,color='g', fontsize=12)
    else:
        #bigax.plot([J1], [J1/J2], [J1/J3],marker='p', color='k')
        print '======tri'
        #bigax.text(J1, (J1/J2)+0.05, (J1/J3)+0.05,"%i"%c,color='k', fontsize=12)
    
    print
    c+=1
    

#%%
    
J1s=np.array(J1s)
J2s=np.array(J2s)
J3s=np.array(J3s)

np.savetxt('RUNI032_griddump.jplots', 
           np.vstack((ids, J1s, J2s, J3s, mass, vols)).T, 
           header='ID\tJ1\tJ2\tJ3\tMass\tnpix',
           fmt='%i\t%7.5f\t%7.5f\t%7.5f\t%6f\t%i')

bigax.axis('equal')
bigax.set_xlim((-1.1,1.1))
bigax.set_ylim((-1.1,1.1))
bigax.xaxis.tick_top()

bigax.grid()
p=bigax.scatter(J1s,J2s,J3s,c=abs(J3s/J2s))
fig.colorbar(p)

#fs=18
#
##bigax.text(0.1,0.8,'Centrally concentrated',color='y',fontsize=fs)
##bigax.text(-0.9,-0.9,'Shell',color='b',fontsize=fs)
##bigax.text(0.5,-0.75,'Filament',color='r',fontsize=fs)
#
#bigax.add_patch(patches.Rectangle((-1.0, -1.0),1,1,alpha=0.2,color='b'))
#bigax.add_patch(patches.Rectangle((0.0, 0.0),1,1,alpha=0.2,color='y'))
#bigax.add_patch(patches.Rectangle((0.0, -1.0),1,1,alpha=0.2,color='r'))
#
bigax.plot([0,0],[0,0],[-1.3,1.3],'k-')
bigax.plot([0,0],[-1.3,1.3],[0,0],'k-')
bigax.plot([-1.3,1.3],[0,0],[0,0],'k-')
#
bigax.set_xlabel(r'$J_{1}$')
bigax.set_ylabel(r'$J_{2}$')
bigax.set_zlabel(r'$J_{3}$')
bigax.xaxis.set_label_position('top') 
#
#fig.subplots_adjust(hspace=0,wspace=0)
#plt.setp([a.get_xticklabels() for a in fig.axes[1:]], visible=False)
#plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible=False)
#
#plt.savefig('proof.pdf')