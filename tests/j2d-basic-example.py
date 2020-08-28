#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:56:24 2020

@author: Sarah Jaffa
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../src/')

import momentsmod as mm

from astropy.io import fits #for reading data
from astrodendro import Dendrogram #for segmenting image



# Read in data
print('====Importing data')
filename='nh2C_mosaic341352'
data, header = fits.getdata(filename+'.fits',header=True)



# Mask interesting region of data
mask_xmin=0
mask_xmax=-1
mask_ymin=1600
mask_ymax=2200


data=data[mask_xmin:mask_xmax, 
              mask_ymin:mask_ymax]


# Compute dendrogram to segment image into regions 
# You can use whatever method you like to define interesting bits of the image
# e.g. thresholding, SCIMES, clumpfind, fellwalker or just masking.
#  Dendrograms is just one that I find interesting

print('====Building dendrogram')
min_value=5e21
min_delta=1e22
min_npix=200

d = Dendrogram.compute(data, 
                       min_value=min_value,
                       min_delta=min_delta,
                       min_npix =min_npix)

#%%

#Plot image
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,6))
im=ax2.imshow(data,interpolation='none',
              origin='lower',
              vmin=min_value)

plt.colorbar(im,ax=ax2)
plt.title(filename)

ax1=mm.plot_moments_axes(ax1,text=True)

#%%

################################################## neat to here
        
nl=sum([len(d.trunk[i].descendants) for i in range(len(d.trunk))])+len(d.trunk)
print('====Analysing %i structures'%nl)
trunks=[t.idx for t in d.trunk]        

structs=np.zeros((nl,3)) #IDs and J values of each structure

for l in d.all_structures:
    idx=l.idx #ID of structure

    grid_2d,npx=mm.make_grid_2d(l) #mask this structure from full image
    J1,J2,com,_ =mm.moments_2d(grid_2d) #calculate moments

    # plot J values
    ax1.plot(J1,J2,'^',label=l.idx,mec='k',mfc='g')
    
    # save structure ID and J values
    structs[idx,:]=[idx,J1,J2]

    if 0: #optionally add structure ID to plot
        ax1.text(J1,J2,"%i"%idx,color='k')

#Save plot
plt.savefig(filename+'_jplot.pdf')

# Save J values to file
outfile=filename+'.jvals'
np.savetxt(outfile, structs, fmt=['%i', '%5.4f', '%5.4f'],
           header=filename+'\nID\tJ1\tJ2')