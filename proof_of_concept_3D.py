#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:44:53 2020

@author: sarah
"""

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
          'ring-thin', 'ring-half',
          'plane-thick','plane-thin','rect-thin',
          'cc',
          'fil-thick','fil-thin','fil-asym',
          'ellipsoid-pro', 'ellipsoid-ob',]#'ellipsoid-ob2','ellipsoid-ob4' 'ellipsoid-tri'
        
fig = plt.figure()
gridspec.GridSpec(5,5)
fsize=8
fig.set_size_inches(w=fsize,h=fsize)

bigax=plt.subplot2grid((5,5), (0,0), colspan=4, rowspan=4, projection='3d')
c=0    
loclist=[(0,4),(1,4),(2,4),(3,4),
         (4,0),(4,1),(4,2),(4,3),(4,4)]


#%%
print('analysing')

J1s=[]
J2s=[]
J3s=[]
mass=[]
vols=[]
ids=[]
c=0
for shape in testlist:
    #print shape,c
    ids.append(c)
    
    #get test shape
    g=mm.testdata3d(20,shape=shape)
    
    #scale to have always 500 mass
    g*=500./np.sum(g)
    
    #calculate centre of mass
    com=spcom(g)
    
    
    Vtot=np.count_nonzero(g) #volume
    Mtot=np.sum(g) #mass
    #print Mtot
    mass.append(Mtot)
    vols.append(Vtot)
    
    #get pixel indices
    yind=np.arange(g.shape[0])
    xind=np.arange(g.shape[1])
    zind=np.arange(g.shape[2])
    
    #calculate moment of inertia
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

    #Eigenvectors of inetia tensor are principle moments
    ow,ov = np.linalg.eig(Inertia_tensor)
    #Sort by size
    eig_ord = np.argsort(ow)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.

    w = ow[eig_ord] #Principle moments
    v = ov[:, eig_ord].T #Vectors

    #Calulate J values
    I0=(2./5.)*Mtot*(((3.*Vtot)/(4.*np.pi))**(2./3.))
    J1, J2, J3 = (I0-w)/(I0+w)
    
    #print (J1, J2, J3)

    J1s.append(J1)
    J2s.append(J2)
    J3s.append(J3)

    #Calculate things to plot (J1,2,3 or I1,2,3, or ratios etc.)
    x=J1
    y=J2
    z=J3
    
    plot_text=False
    lim=0.15 #accuracy limit for checking shape type (symmetrical/prolate/oblate/triaxial)
    
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
    
    #print
    c+=1
    

#%%
    
#Save shape id, J values, mass and volume to file
J1s=np.array(J1s)
J2s=np.array(J2s)
J3s=np.array(J3s)

np.savetxt('3d_testshapes.out', 
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

bigax.set_xlabel(r'$J_{1}$')
bigax.set_ylabel(r'$J_{2}$')
bigax.set_zlabel(r'$J_{3}$')
bigax.xaxis.set_label_position('top') 

plt.savefig('proof3D.pdf')
