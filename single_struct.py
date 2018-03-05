#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 14:54:55 2017

@author: sjaffa

Testing a single structure

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import momentsmod as mm
from astropy.io import fits
from astrodendro import Dendrogram
import matplotlib.gridspec as gridspec
import glob as gb
import sys
import matplotlib.lines as mlines

#================================================fits file
thing=mm.ReadParameters("./SeamusNatSub2D_Params.txt")

print thing['filename']

#if "*" in thing['filename']:
#    #choose file
#    thing['filename']=thing['filename'][:-1]+'2'

if "*" in thing['filename']:
    filelist=gb.glob("%s%s.fits"%(thing['root'],thing['filename']))
    nfiles=len(filelist)
else:
    nfiles=1
    filelist=["%s%s.fits"%(thing['root'],thing['filename'])]
              
#check
if nfiles==0:
    print "Error: files not found, check parameter file."
    sys.exit()

#filelist=["%s%s.fits"%(thing['root'],thing['filename'])]
#nfiles=1
nloop=0

while nloop<nfiles:
    
    
    #%% READ IN DATA, MASK SELECTED AREA AND CHECK DIMENSIONS
    data, header1 = fits.getdata(filelist[nloop],header=True)
     
    data=data[thing['mask_xmin']:thing['mask_xmax'], thing['mask_ymin']:thing['mask_ymax']]
    
    vslice=6
    if 0:#slice data cube along velocity axis
            #threedee=False
            data=data[vslice,:,:]
    
    ##### Ok so here I check if the input was actually a 3D cube rather than a 2D image
    if(len(np.shape(data))==3):
        #data = data.T
        threedee=True
        print "3D"
    
    else:
        threedee=False
        print "2D"
    
    #%% BUILD OR READ IN DENDROGRAM
    #===========================================================================#
    #THIS SEEMS TO BREAK INTERACTIVE PLOT FOR UNKNOWN REASONS. BOLLOCKS         #
    #params=[thing['min_value'],thing['min_delta'],thing['min_npix']]           #
    #d=mm.compute_dend(params, filename=thing['filename'],                      #
    #                  root=thing['root'], data=data)                           #
    #===========================================================================#
    d = Dendrogram.compute(data, 
                           min_value=thing['min_value'],
                           min_delta=thing['min_delta'],
                           min_npix =thing['min_npix'])
    
    #%% SET UP PLOT
    
        
    nstruct=0#input("Enter ID of structure to test:")
    nslice=10
    
    l=d[nstruct]
    
    idx=l.idx           #ID of structure
    if threedee:
        x=l.indices()[1]    #x coordinates of structure
        y=l.indices()[2]    #y coordinates of structure
    else:
        x=l.indices()[0]    #x coordinates of structure
        y=l.indices()[1]    #y coordinates of structure    
    v=l.values()        #pixel values of structure
    n=l.get_npix()      #number of pixels
                        
    xmax=np.amax(x)+1   #coordinates of bounding box
    ymax=np.amax(y)+1
    xmin=np.amin(x)-1
    ymin=np.amin(y)-1
    
    vmin=thing['min_value']#np.amin(v) #bounding instensity values
    vmax=10#np.amax(v)
    
    grid=np.zeros((xmax-xmin,ymax-ymin))
        
    for j in range(n):
        grid[x[j]-xmin,y[j]-ymin]+=v[j]
    #%%
    
    grid=data
    j1s=[]
    j2s=[]
    intensities=np.logspace(np.log10(vmin),np.log10(vmax),nslice)
    
    #plt.figure(figsize=(4,12))
    #im=plt.imshow(data.T,extent=(-0.5,0.5,0,3),cmap='jet',origin='lower')#,vmin=vmin,vmax=vmax)
    #v = plt.axis()
    #plt.contour(data.T,levels=intensities,colors='k',extent=(-0.5,0.5,0,3))#,vmin=vmin,vmax=vmax)
    #plt.axis(v)
    #plt.colorbar(im)
    ##plt.axes('equal')
    
    #    fig=plt.figure(figsize=(20,10))
    #    gridspec.GridSpec(5,10)
    #    bigax=plt.subplot2grid((5,10), (0,0), colspan=10, rowspan=3)
    #    fig.subplots_adjust(hspace=0,wspace=0)
    
    #Analyse to nslice intensity slices:
    i=-1
    I=vmin
    jj=0
    fig=plt.figure(figsize=(8,6))
    #while I<vmax:
    for I in intensities:
        try:
            i+=1
            print I,i
            
            #if np.sum(grid>I)<thing['min_npix']:
            #    break
            #intensities.append(I)
            levelgrid=np.zeros_like(grid)
            levelgrid[grid>I]=grid[grid>I]
                        
            #calculate moments and plot
            J1,J2,_,_ =mm.moments_2(levelgrid)
            j1s.append(J1)
            j2s.append(J2)
            
            if i in [0,1,2,3]: #save plot of individual structures
                ax=plt.subplot2grid((4,1), (jj,0))
                #plt.figure(figsize=(1,1))
                #ax=plt.gca()
                im=ax.imshow(levelgrid,interpolation='none',
                             origin='lower', cmap='hot',norm=LogNorm(),
                             vmin=1,vmax=10**3)
                             #extent=(-0.5,0.5,0,3))
                ax.text(0.8,0.05,"%1.1e"%I,color='k',transform=ax.transAxes,)
                #ax.plot(com[1],com[0], 'ko',)
                #plt.savefig("%sStructures_%s/%i.png"%(thing['root'],thing['filename'],idx))
                #plt.close()
                ax.set_ylim(15,65)
                ax.set_xticks([])
                ax.set_yticks([])
                jj+=1
            if 0: #save data of individual structures
                filename="%sStructures_%s/%i.dat"%(thing['root'],thing['filename'],idx)
                np.savetxt(filename, levelgrid)
            #I+=thing['min_delta']
        except IndexError:
            break

    
    fig.subplots_adjust(top=0.95,bottom=0.22)
    cbar_ax = fig.add_axes([0.2, 0.13, 0.6, 0.05])
    cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label(r'Column density (g cm$^{-2}$)')
        
    fig=plt.figure(7)
    bigax=plt.gca()
    
    bigax.set_xlabel(r'Column density (g cm$^{-2}$)')
    bigax.set_ylabel("$J$")
    #bigax.xaxis.set_label_position('top')
    #bigax.xaxis.tick_top()
    bigax.set_title('Subsonic')
     
    
    bigax.semilogx(intensities, j1s, 'b-s', label=r'$J_{1}$')
    bigax.semilogx(intensities, j2s, 'g-^', label=r'$J_{2}$')
  
    bigax.hlines(0,vmin,vmax,'k')
    #bigax.legend(loc='center right')
    
    bigax.set_xlim((vmin,I))
    bigax.set_ylim((-1.1,1.1))
    
    #ax.colorbar(im)
    #f,jax=mm.plot_moments_axes()
    #jax.plot(j1s,j2s,'^-')#,alpha=0.3)

    nloop+=1
fig.subplots_adjust(bottom=0.2)
line_j1 = mlines.Line2D([], [], color='blue', marker='s',
                          markersize=15, label=r'$J_{1}$')
line_j2 = mlines.Line2D([], [], color='green', marker='^',
                          markersize=15, label=r'$J_{2}$')
plt.legend(handles=[line_j1,line_j2],loc='center right')
