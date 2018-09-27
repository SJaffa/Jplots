#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:24:29 2017

@author: sjaffa
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import momentsmod as mm
import os
from astropy.io import fits
from astrodendro import Dendrogram


import glob as gb
import sys


#Read in parameters file
thing=mm.ReadParameters("./example_Params.txt")


if "*" in thing['filename']: #Many fits files to run
    filelist=gb.glob("%s%s.fits"%(thing['root'],thing['filename']))
    nfiles=len(filelist)
else: #Just one fits file
    nfiles=1
    filelist=["%s%s.fits"%(thing['root'],thing['filename'])]
              
#check
if nfiles==0:
    print "Error: files not found, check parameter file."
    sys.exit()
        
nloop=0 #loop counter for multiple files
all_structs=np.zeros((1,6)) #info on j-plot coordinates of each structure and if the structure is a branch, trunk or leaf
all_npix=[]
clump_mass=[]

vslice=0

while nloop<nfiles:

#while vslice<datao.shape[0]:
    print filelist[nloop]
    # READ IN DATA, MASK SELECTED AREA AND CHECK DIMENSIONS
    datao, header1 = fits.getdata(filelist[nloop],header=True)
           
    #Check if the input was a 3D cube or 2D image
    data,threedee=mm.check_dimensions(datao)

    if 0:#slice data cube along velocity axis
        threedee=False
        data=datao[vslice,:,:]
    
    data=data[thing['mask_xmin']:thing['mask_xmax'], 
              thing['mask_ymin']:thing['mask_ymax']]
              
    # BUILD DENDROGRAM
    print "====Building dend"

    d = Dendrogram.compute(data, 
                           min_value=thing['min_value'],
                           min_delta=thing['min_delta'],
                           min_npix =thing['min_npix'])
    if len(d.trunk)>0:
        # SET UP PLOT
    
        bigax,dendax,f=mm.plot_interactive_axes(d,threedee=threedee)
        
        if threedee:
            im=dendax.imshow(np.sum(data,axis=0),interpolation='none',
                             origin='lower',#norm=LogNorm(),cmap='jet',
                             vmin=thing['min_value'],)# vmax=1e3)    
        else:                      
            im=dendax.imshow(data,interpolation='none',
                             origin='lower', #norm=LogNorm(),cmap='jet',
                             vmin=thing['min_value'],)# vmax=1e3)

        
        plt.colorbar(im,ax=dendax)
        plt.title(filelist[nloop])
        
            
        # CHECK HOW MANY STRUCTURES   
        nl=sum([len(d.trunk[i].descendants) for i in range(len(d.trunk))])+len(d.trunk)
        print '====Testing %i structures'%nl
        trunks=[t.idx for t in d.trunk]        

        
        desc=np.zeros((nl,2))    #info on descendants of each structure
        structs=np.zeros((nl,6)) #info on j-plot coordinates of each structure and if the structure is a branch, trunk or leaf
        
        for l in d.all_structures:
            idx=l.idx           #ID of structure
            if threedee:
                grid_2d,grid_3d,npx=mm.make_grid_3d(l)
                mass_kg, mass_msun = mm.calculate_mass(grid_3d,
                                                       thing['pixel_pc'],
                                                       thing['pixel_vel'])
                clump_mass.append(mass_msun)
                I1,I2,I3,com,_ =mm.moments_3d(grid_2d)
                #print idx,mass_msun
            else:
                grid_2d,npx=mm.make_grid_2d(l)
                clump_mass.append(0)
                I1,I2,com,_ =mm.moments_2d(grid_2d)
            
            all_npix.append(npx)
            
            #calculate moments and plot
            
        

            if 0: #save plot of individual structures
    
                ##### Just a check to see if the Structures folder exists, if it isn't we make it and if it is we delete it then make it again :D 
                checkname = filelist[nloop]+"_Structures"#thing['root'] + "Structures_"+thing['filename']
                
                if(os.path.isdir(checkname)):
                    cmd = "rm -r " + checkname
                    os.system(cmd)
                
                cmd = "mkdir " + checkname
                os.system(cmd)
                #Fill with new figures
                plt.figure(figsize=(4,4))
                ax=plt.gca()
                ax.imshow(grid_2d,interpolation='none',origin='lower', cmap='rainbow',norm=LogNorm())
                ax.text(0.05,0.90,"%i"%idx,color='k',transform=ax.transAxes,)
                ax.plot(com[1],com[0], 'ko',)
                plt.savefig(filelist[nloop]+"_Structures/%i.png"%idx)
                plt.close()
    
                
            if 0: #save data of individual structures
                filename=filelist[nloop]+"_Structures/%i.dat"%idx
                np.savetxt(filename, grid_2d)
        
       
            #save id and location for plotting merger paths later
    
            #check trunk/leaf
            if l.is_leaf:
                if l.idx in trunks:
                    bigax.plot(I1,I2,'^',label=l.idx,mec='k',mfc='g',picker=4,gid=str(idx))
                    #print idx, l.get_peak()
                else:
                    bigax.plot(I1,I2,'^',label=l.idx,mec='k',mfc='white',picker=4,gid=str(idx))
                structs[idx,:]=[idx,I1,I2,1,l.get_npix(),clump_mass[-1]]
            elif l.idx in trunks:
                bigax.plot(I1,I2,'s',label=l.idx,mec='k',mfc='k',picker=4,gid=str(idx))
                desc[idx,:]=[tl.idx for tl in l.children]
                structs[idx,:]=[idx,I1,I2,2,l.get_npix(),clump_mass[-1]]
            else: #branch
                bigax.plot(I1,I2,'o',label=l.idx,mec='k',mfc='grey',picker=4,gid=str(idx))
                desc[idx,:]=[tl.idx for tl in l.children]
                structs[idx,:]=[idx,I1,I2,3,l.get_npix(),clump_mass[-1]]
            if 0: #add structure ID to plot?
                bigax.text(I1,I2,"%i"%idx,color='k')


            
        if 0: #save results file with data on structures in this fits file
            resfile="%s_RES_%s_%s_%s"%(filelist[nloop],thing['min_value'],thing['min_npix'],thing['min_delta'])

            if not (os.path.isdir(thing['root']+resfile)):
                print"===========SAVING"
                head="""%s
    Dendrogram parameters: 
                min_value =%s
                min_npix  =%s
                min_delta =%s\n\n
    ID\tJ1\tJ2\ttype\tnpix\n"""%(filelist[nloop],
                                 thing['min_value'],
                                 thing['min_npix'],
                                 thing['min_delta'])
                np.savetxt(resfile, structs, 
                           fmt="%i\t%4.3f\t%4.3f\t%i\t%i",
                           header=head)
        
        if 0: #save leaves
            resfile="%s_RES_%s_%s_%s"%(filelist[nloop],thing['min_value'],thing['min_npix'],thing['min_delta'])
            leaves=structs[structs[:,3]==1]
            leaves=np.concatenate((leaves[:,:3],leaves[:,4:]),axis=1)
            leaves=np.pad(leaves, ((0,0),(0,4)),'constant',constant_values=0)
            for n in range(leaves.shape[0]):
                n_id=leaves[n,0]
                peak=d[n_id].get_peak()
                leaves[n,5]=peak[1] #peak val
                leaves[n,6:]=peak[0] #peak v,x,y
                
            print"===========SAVING"
            head="""%s
Dendrogram parameters: 
    min_value =%s
    min_npix  =%s
    min_delta =%s\n\n
ID\tJ1\tJ2\tnpix\tMass(m_sun)\tTpeak\tpeak_v\tpeak_x\tpeak_y\n"""%(thing['filename'],
                             thing['min_value'],
                             thing['min_npix'],
                             thing['min_delta'])
            np.savetxt(resfile, leaves, 
                       fmt="%i\t%4.3f\t%4.3f\t%i\t%.2f\t%.2f\t%i\t%i\t%i",
                       header=head)
        
        all_structs=np.concatenate((all_structs,structs), axis=0)
    nloop+=1 
    vslice+=1
    
if 0: #save results file with data on structures in all fits files

    resfile="%s_RES_%s_%s_%s"%(thing['filename'],thing['min_value'],thing['min_npix'],thing['min_delta'])
    
    if not (os.path.isdir(thing['root']+resfile)):
        np.savetxt(thing['root']+resfile, structs)
        
if 0:  #plot tree
    p = d.plotter()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    p.plot_tree(ax, color='blue')
    ax.set_ylabel(r'Column density (H$_{2}$ cm$^{-2}$)')
    ax.set_xticks([])

if 1: #plots whole population
    print "====Building kde"
    all_structs=all_structs[1:,:]
    kdefig,kdeax=mm.plot_moments_axes()
    kdefig.set_size_inches(7,6)
    ms=8
    maxmass=np.amax(all_structs[:,4])
    plt.text(-0.8,0.8,'(a)')
    
    if 0:#plot all, coloured by npix
        cmap = plt.get_cmap('RdYlBu')
        scales=all_structs[:,4]/maxmass
        mapcol=cmap(scales)
        im=kdeax.scatter(all_structs[:,1],all_structs[:,2], marker='^',facecolor='white',s=20+(200*scales))
        #plt.colorbar()
    if 1:#plot leaves
        leaves=all_structs[all_structs[:,3]==1] 
        kdeax.plot(leaves[:,1],leaves[:,2], '^',mec='k',fillstyle='none',markersize=ms)
    if 1:#plot branches
        branch=all_structs[all_structs[:,3]==3] 
        kdeax.plot(branch[:,1],branch[:,2], 'o',mec='k',mfc='grey',markersize=ms)
    if 1:#plot trunks
        trunks=all_structs[all_structs[:,3]==2] 
        kdeax.plot(trunks[:,1],trunks[:,2], 's',mec='k',mfc='k',markersize=ms)
    if 0: #plot KDE
        k1=mm.KDE_2D(all_structs,1)
        im=kdeax.contour(k1, origin='lower', extent=(-1,1,-1,1), cmap='gray', 
                         levels=np.logspace(-1,1,5))
        #plt.clabel(im, inline=1, fontsize=10)
    if 0: #plot merger paths
        for ii in range(nl):
            ds=desc[ii]
            if ds[0]!=0:
                for dj in ds:
                    j=int(dj)
                    xi=structs[ii,1]
                    yi=structs[ii,2]
                    xj=structs[j,1]
                    yj=structs[j,2]
                    kdeax.plot([xi,xj],[yi,yj],'k-',alpha=0.3)
    kdeax.set_xlim((-1.1,1.1))
    kdeax.set_ylim((-1.1,1.1))

    #kdefig.savefig(filelist[nloop-1][:-5]+"_kde.pdf")
    
if 0: #analyse masses
    leaves=all_structs[all_structs[:,3]==1]
    print "N leaves = %i"%leaves.shape[0]
    leaf_mass=leaves[:,5]
    print "Mean mass (M_sun): %f\nMax: %f\nMin: %f\n"%(np.mean(leaf_mass),
                                               np.amax(leaf_mass),
                                               np.amin(leaf_mass))




