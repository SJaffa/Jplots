#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 14:08:47 2017

@author: sjaffa
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astrodendro import Dendrogram
from momentsmod import *
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LogNorm


rcw120={'filename':'nh2C_mosaic341352',
        'root':'/home/sjaffa/Dropbox/PhDwork/170117Moments/HiGAL_tiles/RCW120/',
        'mask_xmin':0,
        'mask_xmax':-1,
        'mask_ymin':1600,
        'mask_ymax':2200,
        'min_value':5e21,
        'min_delta':1e22,
        'min_npix':200}
        
liz_ridge={'filename':'l316_cdens2',
           'root':'/home/sjaffa/Dropbox/PhDwork/170117Moments/Liz_ridge/',
            'mask_xmin':640,
            'mask_xmax':797,
            'mask_ymin':642,
            'mask_ymax':803,
            'min_value':3e2,
            'min_delta':3e2,
            'min_npix':5}
        
        
thing=rcw120

print "Reading data"
        
data = fits.getdata(thing['root']+thing['filename']+'.fits')

if (1 in data.shape):
    if data.shape[0]==1:
        data=data[0,:,:]
    elif data.shape[1]==1:
        data=data[:,0,:]
    elif data.shape[2]==1:
        data=data[:,:,0]

data=data[thing['mask_xmin']:thing['mask_xmax'], thing['mask_ymin']:thing['mask_ymax']]

print "Computing dendrogram"
d = Dendrogram.compute(data,
                       min_value=thing['min_value'],
                       min_delta=thing['min_delta'],
                       min_npix=thing['min_npix'])
#d.save_to(thing['root']+thing['filename']+'_dend.fits')

#%%

#d = Dendrogram.load_from(thing['root']+thing['filename']+'_dend.fits')
def onpick(event):
    print "name: ",event.name
    print "canvas: ",event.canvas
    print "guievent: ",event.guiEvent
    print "mouseevent: ", event.mouseevent
    print "artist: ", event.artist
#    datax,datay = event.artist.get_data()
#    print datax, len(datay)
#    msx, msy = event.mouseevent.xdata, event.mouseevent.ydata
#    dist = np.sqrt((np.array(datax)-msx)**2+(np.array(datay)-msy)**2)
#    ind = np.argmin(dist)
#    print ind
#    print dist[ind]
    
    ind = event.ind
    s = event.artist.get_gid()
    print s, ind
    #print len(ind)
#    if len(ind)>1:
#        print "lots"
#        datax,datay = event.artist.get_data()
#        datax,datay = [datax[i] for i in ind],[datay[i] for i in ind]
#        msx, msy = event.mouseevent.xdata, event.mouseevent.ydata
#        dist = np.sqrt((np.array(datax)-msx)**2+(np.array(datay)-msy)**2)
#        ind = [ind[np.argmin(dist)]]
#        print ind
    #ind=int(ind)
    #s = int(event.artist.get_gid())
    #print "Structure %i, level %i:%.2e"%(s, ind, levels[s,ind])
    #print s, ind
    p.plot_contour(ax, structure=int(s), lw=3, label='Structure %s'%s)#, colors='red')
    plt.legend()
    
    
n_cuts=2
nl=sum([len(d.trunk[i].descendants) for i in range(len(d.trunk))])+len(d.trunk)
print 'Testing %i structures with %i cuts'%(nl, n_cuts)

I1s=np.zeros((nl,n_cuts))
I2s=np.zeros_like(I1s)
levels=np.zeros_like(I1s)

ids=np.zeros(nl)
leafs=np.zeros_like(ids)
#n_levs=[]
struct=[]


#lc=-1
for l in d.all_structures:
    #lc+=1
    #if lc>2:
    #    break
    s=l.idx
    struct.append(s)
    print "Analysing structure %i"%s
    ids[s]=s
    leafs[s]=l.is_leaf
    #read in dimensions
    xx=l.indices()[0]
    yy=l.indices()[1]
    #zz=l.indices()[2]
    v=l.values()
    
    #rearrange dimensions
    x=xx#zz
    y=yy
    #z=xx
    
    n=l.get_npix()
    xmax=np.amax(x)+1
    ymax=np.amax(y)+1
    xmin=np.amin(x)-1
    ymin=np.amin(y)-1
    #zmax=np.amax(z)+1
    #zmin=np.amin(z)-1
    
    if xmax==0:
        xmax=zmax
        x=z
    grid=np.zeros((xmax-xmin,ymax-ymin))

    for j in range(n):
        grid[x[j]-xmin,y[j]-ymin]=v[j]
    
    vmax=l.vmax
    vmin=l.vmin
    #vrange=vmax-vmin
    #n_levels=max((1,np.floor(vrange/thing['min_delta'])))
    #print n_levels
    #n_levs.append(n_levels)
 
    if leafs[s]:
        #skip top 10% of leaves to avoid "dot"
        vmax=0.9*vmax
        
    levels[s,:]=np.logspace(np.log10(vmin),np.log10(vmax),n_cuts)
    
    if 0:
        plt.figure()
        im=plt.imshow(grid, origin='lower',interpolation='None')
        plt.contour(grid,levels[s,:],colors='white')
        plt.colorbar(im)
        plt.title("Structure %i"%s)
        plt.savefig('%s%s_struct%i.png'%(thing['root'],thing['filename'],s))
        plt.close()
    
    
    vcut=levels[s,0]
    gridcut=np.zeros_like(grid)
    gridcut[grid>vcut]=grid[grid>vcut]
    I1s[s,0],I2s[s,0]=moments(gridcut)

    #if leafs[s]:
    #    bigax.plot(I1o,I2o,'^',color=colors(s),picker=5)5
    #else:
    #    bigax.plot(I1o,I2o,'o',color=colors(s),picker=5)
    for c in range(1,n_cuts):#vcut in levels[s,1:]:
        vcut=levels[s,c]
        gridcut=np.zeros_like(grid)
        gridcut[grid>vcut]=1
        
        I1s[s,c],I2s[s,c]=moments(gridcut)
        #if leafs[s]:
        #    line,=bigax.plot(I1n,I2n,'^',color=colors(l.idx),picker=5)
        #else:
        #    line,=bigax.plot(I1n,I2n,'o',color=colors(l.idx),picker=5)
        #bigax.plot([I1o,I1n],[I2o,I2n],'-',
        #           color=colors(l.idx),alpha=0.3)
        #I1o,I2o=I1n,I2n

dI=np.sqrt([((I1s[:,i]-I1s[:,i+1])**2 + (I2s[:,i]-I2s[:,i+1])**2) for i in np.arange(n_cuts-1)]).T
 
#%%           
bigax,fig=plot_moments_axes(d)
colors = cm.get_cmap('terrain', nl)

# make space for colorbar:
fig.tight_layout(rect=[0, 0, 0.5, 1])

# adding colorbar:
ax_dend = fig.add_axes([0.5, 0.05, 0.45, 0.85])
#norm = mpl.colors.Normalize(vmin=0, vmax=nl)
#cb = mpl.colorbar.ColorbarBase(ax_cb, cmap='terrain', 
#                               norm=norm, orientation='vertical')

for s in range(nl):
    bigax.plot(I1s[s,:],I2s[s,:],'-',color=colors(s),alpha=0.7,picker=4,gid=str(s))
    for i in range(n_cuts-1):
        if dI[s,i]<0.05:
            if  leafs[s]:
                bigax.plot(I1s[s,i:i+1],I2s[s,i:i+1],'^',color=colors(s))#,picker=4,gid=str(s))
            else:
                bigax.plot(I1s[s,i:i+1],I2s[s,i:i+1],'o',color=colors(s))#,picker=4,gid=str(s))
        
    
p = d.plotter()
ax = ax_dend
im=ax.imshow(data, origin='lower', interpolation='nearest',norm=LogNorm(), 
             cmap=plt.cm.get_cmap('Greys'))#,vmin=thing['min_value'])
# Show contour for ``min_value``
#p.plot_contour(ax, color='black')

# Highlight two branches
#p.plot_contour(ax, structure=8, lw=3, colors='red')
#p.plot_contour(ax, structure=34, lw=3, colors='blue')
#p.plot_contour(ax, structure=5, lw=3, colors='green')
#p.plot_contour(ax, structure=47, lw=3, colors='yellow')
fig.canvas.mpl_connect('pick_event', onpick)
