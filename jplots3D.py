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

testlist=['sphere','shell-thin','shell-thick','cc']
          
fig = plt.figure(1)
gridspec.GridSpec(5,5)
fig.set_size_inches(w=10,h=10)

bigax=plt.subplot2grid((5,5), (0,0), colspan=4, rowspan=4, projection='3d')
c=0    
loclist=[(0,4),(1,4),(2,4),(3,4),
         (4,0),(4,1),(4,2),(4,3),(4,4)]

for n in testlist:
    #print n    
    g=mm.testdata3d(100,shape=n)
    #g[0,0]=0
    smallax=plt.subplot2grid((5,5), loclist[c])
    #plt.figure()
    #smallax=plt.gca()
    im=smallax.imshow(np.sum(np.sqrt(g.T), axis=0),interpolation='None',
                      origin='lower',cmap='Greys',vmin=0)
    plt.setp( smallax.get_xticklabels(), visible=False)
    plt.setp( smallax.get_yticklabels(), visible=False)
    smallax.text(5,85,"%i"%c,color='k', fontsize=14)
    #plt.colorbar(im)
    
    if n[:2]=='cc':
        mk='s'
    else:
        mk='o'
    com=mm.find_com_3d(g)
    #smallax.plot([com[0]],[com[1]],'y^')

    com=spcom(g)
    Vtot=np.count_nonzero(g)
    Mtot=np.sum(g)
    
    xind=np.arange(g.shape[0])
    yind=np.arange(g.shape[1])
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
    
    w,v = np.linalg.eigh(Inertia_tensor)
    
    print w,v
    
    #I1,I2,t1=mm.i12(g,com)
    first=w[0]
    second=w[1]
    third=w[2]
    
    I0=(2./5.)*Mtot*(((3.*Vtot)/(4.*np.pi))**(2./3.))
    
    J1, J2, J3 = (I0-w)/(I0+w)
    #print II1,II2
    
    
    bigax.plot([J1], [J2], [J3],marker=mk,label=n,)# color='k')
    #bigax.text(J1, J2+0.05, J3+0.05,"%i"%c,color='k', fontsize=14)
    print
    c+=1
    
#    print "%i\t%i\t%i\t%3.2e\t%3.2e\t%3.2f\t%3.2f"%(c, Atot, Mtot, first, 
#                                                    second, II1, II2)
    
bigax.axis('equal')
bigax.set_xlim((-1.1,1.1))
bigax.set_ylim((-1.1,1.1))
bigax.xaxis.tick_top()

bigax.grid()

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