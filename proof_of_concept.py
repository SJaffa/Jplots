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
import momentsmod as mm

testlist=[#'diag-fil','diag-fil2',
          'disk','disk-cc',#'peak',
          #'noisy',#'dot',
          'fil-thin','fil-thick2',#'fil-thin2',
          'ellipse','ellipse-cc',#'disk-cc2',
          #'half-ring-thick','half-ring-thin',
          'ring-thick','ring-thin',
          'ring-lumpy',#'asym_cc','fil-thin_dot',
          ]
          
fig = plt.figure(1)
gridspec.GridSpec(5,5)
fig.set_size_inches(w=10,h=10)

bigax=plt.subplot2grid((5,5), (0,0), colspan=4, rowspan=4)
c=0    
loclist=[(0,4),(1,4),(2,4),(3,4),
         (4,0),(4,1),(4,2),(4,3),(4,4)]

for n in testlist:
    #print n    
    g=mm.testdata(100,shape=n)
    #g[0,0]=0
    smallax=plt.subplot2grid((5,5), loclist[c])
    #plt.figure()
    #smallax=plt.gca()
    im=smallax.imshow(np.sqrt(g.T),interpolation='None',
                      origin='lower',cmap='Greys',vmin=0)
    plt.setp( smallax.get_xticklabels(), visible=False)
    plt.setp( smallax.get_yticklabels(), visible=False)
    smallax.text(5,85,"%i"%c,color='k', fontsize=14)
    #plt.colorbar(im)
    
    if n[:2]=='cc':
        mk='s'
    else:
        mk='o'
    com=mm.find_com(g)
    #smallax.plot([com[0]],[com[1]],'y^')

    Atot=np.count_nonzero(g)
    Mtot=np.sum(g)
    I1,t1=mm.i1(g,com)
    I2,_=mm.i2(g,com)
    #print "Area=%.2e,M=%.2e"%(Atot,Mtot)
    #print I1,I2
    
    #print t1
    
    first=min([I1,I2]) #first moment must always be smallest
    second=max([I1,I2])
    
    
    ma=Atot*Mtot
    mI= ma/(4*np.pi)
    
    #print mI
    
    II1=(ma - 4*np.pi*first)/(ma + 4*np.pi*first)
    II2=(ma - 4*np.pi*second)/(ma + 4*np.pi*second)
    
    #print II1,II2
    
    
    bigax.plot(II1,II2,marker=mk,label=n, color='k')
    bigax.text(II1,II2+0.05,"%i"%c,color='k', fontsize=14)
    print
    c+=1
    
    print "%i\t%i\t%i\t%3.2e\t%3.2e\t%3.2f\t%3.2f"%(c, Atot, Mtot, first, 
                                                    second, II1, II2)
    
bigax.axis('equal')
bigax.set_xlim((-1.1,1.1))
bigax.set_ylim((-1.1,1.1))
bigax.xaxis.tick_top()

bigax.grid()

fs=18

#bigax.text(0.1,0.8,'Centrally concentrated',color='y',fontsize=fs)
#bigax.text(-0.9,-0.9,'Shell',color='b',fontsize=fs)
#bigax.text(0.5,-0.75,'Filament',color='r',fontsize=fs)

bigax.add_patch(patches.Rectangle((-1.0, -1.0),1,1,alpha=0.2,color='b'))
bigax.add_patch(patches.Rectangle((0.0, 0.0),1,1,alpha=0.2,color='y'))
bigax.add_patch(patches.Rectangle((0.0, -1.0),1,1,alpha=0.2,color='r'))

bigax.vlines(0,-1.3,1.3)
bigax.hlines(0,-1.3,1.3)

bigax.set_xlabel(r'$J_{1}$')
bigax.set_ylabel(r'$J_{2}$')
bigax.xaxis.set_label_position('top') 

fig.subplots_adjust(hspace=0,wspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[1:]], visible=False)
plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible=False)

plt.savefig('proof.pdf')