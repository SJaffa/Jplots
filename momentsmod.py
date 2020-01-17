#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 10:45:04 2017

@author: sjaffa
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass as spcom

def ColDens_13co(grid,dv):
    """
    Convert grid of brightness temp (PPV) to column density of 13CO (PP)
    
    Parameters
    ----------
    grid : array_like
        3D array of brightness temperature
    dv :
        velocity resolution

    Returns
    -------
    N_CO : array_like
        2D array of projected 13CO column density
    """
    
    T_ex = 25 #excitation temperature, K
    bracket = (1./(np.exp(10.6/T_ex)-1))-0.02
    bracket2 = (T_ex*np.exp(5.3/T_ex))/(1-np.exp(-10.6/T_ex))
    #print bracket, bracket2
    N_CO=np.zeros((grid.shape[1],grid.shape[2]))
    
    for i in range(grid.shape[1]):
        for j in range(grid.shape[2]):
            tau_int=0
            for k in range(grid.shape[0]):
                #if grid[i,j,k]>0:
                tau_int += (-np.log(1 - (grid[k,i,j]/10.6)*(bracket**-1)))*dv
            #print tau_int
            N_CO[i,j] = 1.51E14*bracket2*tau_int
    return N_CO

def calculate_mass(grid,pix_pc,pix_vel):
    X=0.735 #Conversion factor from H2 mass to total mass
    mu_H2 =  3.35E-27 #Mass of hydrogen molecule in kg
    px_cm = pix_pc*3.1E18 #pixel width in cm
    A_px_cm = px_cm**2 #pixel area, in cm^2
    
    #print px_cm,A_px_cm
    N_13co = ColDens_13co(grid,pix_vel)
    #print N_13co
                
    N_H2 = 6.5E5 * N_13co
    m_H2 = mu_H2*np.sum(N_H2) * A_px_cm
    
    mass = m_H2/X #mass in kg
    solar_mass=mass/(1.99E30)
    return mass,solar_mass
    
def make_grid_3d(l):

    vel=l.indices()[0]  #velovity coordinates of structure
    x=l.indices()[1]    #x coordinates of structure
    y=l.indices()[2]    #y coordinates of structure

    v=l.values()        #pixel values of structure
    n=l.get_npix()      #number of pixels
    

    xmax=np.amax(x)+1   #coordinates of bounding box
    ymax=np.amax(y)+1
    xmin=np.amin(x)-1
    ymin=np.amin(y)-1
    vmin=np.amin(vel)-1
    vmax=np.amax(vel)+1

    #make empty grid
    grid_3d=np.zeros((vmax-vmin,xmax-xmin,ymax-ymin))

    #populate 3d grid
    for j in range(n):
        grid_3d[vel[j]-vmin,x[j]-xmin,y[j]-ymin]=v[j]
        
    #integrate velocity axis for 2d grid
    grid_2d=np.sum(grid_3d,axis=0)
        
    return grid_2d,grid_3d,n
    
def make_grid_2d(l):

    x=l.indices()[0]    #x coordinates of structure
    y=l.indices()[1]    #y coordinates of structure

    v=l.values()        #pixel values of structure
    n=l.get_npix()      #number of pixels
    
    xmax=np.amax(x)+1   #coordinates of bounding box
    ymax=np.amax(y)+1
    xmin=np.amin(x)-1
    ymin=np.amin(y)-1

    #make empty grid
    grid=np.zeros((xmax-xmin,ymax-ymin))

    #populate 2d grid
    for j in range(n):
        grid[x[j]-xmin,y[j]-ymin]=v[j]
        
    return grid,n

def justMST(x,y,z=[8]):
    from scipy.sparse.csgraph import minimum_spanning_tree
    if len(z)==1: #2D data
        grid=np.zeros((len(x),len(y)))
        all_edges=[]
        for i in range(len(x)):
          for j in range(len(x)):
            if i!=j:
              if grid[j][i]==0:
                grid[i][j]=np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
                all_edges.append(grid[i][j])   
        Tcsr = minimum_spanning_tree(grid,overwrite=True)
        mst=Tcsr.toarray().astype(float)
    
    elif len(z)==len(x):
        length=len(x)
        grid=np.zeros((length,length))
        all_edges=[]
        for i in range(length):
          for j in range(length):
            if i!=j:
              if grid[j][i]==0:
                grid[i][j]=np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
                all_edges.append(grid[i][j])   
        #make mst
        Tcsr = minimum_spanning_tree(grid,overwrite=True)
        mst=Tcsr.toarray().astype(float)
    else:
        print "Error: x,y, and z must have the same length"
        quit()
    return mst, all_edges

def plot_flat_mask(ax,d,idx,threedee,col='w'):
    if threedee:
        mask_3d=d[idx].get_mask()
        mask_2d=np.sum(mask_3d,axis=0)
    else:
        mask_2d=d[idx].get_mask()
    ax.contour(mask_2d,colors=col,levels=[0])

def single_structure_evolution(d,struct_id,nslices=20):
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LogNorm
    
    l=d[struct_id]      #get structure from dendrogram

    x=l.indices()[0]    #x coordinates of structure
    y=l.indices()[1]    #y coordinates of structure
    v=l.values()        #pixel values of structure
    n=l.get_npix()      #number of pixels
                        
    xmax=np.amax(x)+1   #coordinates of bounding box
    ymax=np.amax(y)+1
    xmin=np.amin(x)-1
    ymin=np.amin(y)-1
    
    vmin=np.amin(v) #bounding instensity values
    vmax=np.amax(v)
    
    grid=np.zeros((xmax-xmin,ymax-ymin)) #blank grid
        
    for j in range(n): #fill in structure
        grid[x[j]-xmin,y[j]-ymin]+=v[j]

    j1s=[] 
    j2s=[]
    
    intensities=np.logspace(np.log10(vmin),np.log10(vmax),nslices)

    #set up plot
    fig=plt.figure(figsize=(20,10))
    gridspec.GridSpec(5,10)
    bigax=plt.subplot2grid((5,10), (0,0), colspan=10, rowspan=3)
    fig.subplots_adjust(hspace=0,wspace=0)
    
    #Analyse to nslices intensity slices:
    i=-1
    for I in intensities:
        i+=1
        #print I
        
        #select above intensity cutoff
        levelgrid=np.zeros_like(grid)
        levelgrid[grid>I]=grid[grid>I]
        
        #plot thresholded structure
        ax=plt.subplot2grid((5,10), (3+(i/10),i%10))
        ax.imshow(levelgrid,interpolation='none',origin='lower', cmap='rainbow',norm=LogNorm())
        ax.set_xticks([])
        ax.set_yticks([])
                    
        #calculate moments
        J1,J2,_,_ =moments_2d(levelgrid)
        j1s.append(J1)
        j2s.append(J2)
        
    #plot J values vs intensity
    bigax.set_xlabel("Intensity")
    bigax.set_ylabel("J")
    bigax.xaxis.set_label_position('top')
    bigax.xaxis.tick_top()
    bigax.set_xlim((vmin,vmax))
    
    bigax.semilogx(intensities, j1s, 'b:x', label="J1")
    bigax.semilogx(intensities, j2s, 'g:x', label="J2")

    bigax.hlines(0,vmin,vmax,'k')
    bigax.legend(loc='lower right')
    
    return j1s,j2s,I,fig
    
def find_com_2d(grid):
    com=np.array([0,0])
    nx=grid.shape[0]
    ny=grid.shape[1]
    for i in range(nx):
        for j in range(ny):
            com =com + grid[i,j]*np.array([i,j])
    com=com/np.sum(grid)

    return com

def find_com_3d(grid):
    com=np.array([0,0,0])
    nx, ny, nz=grid.shape[:]

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                com =com + grid[i,j, k]*np.array([i,j, k])
    com=com/np.sum(grid)

    return com

def check_dimensions(data):
    print data.shape
    if 0 in data.shape:
        print "zero dimension"
        i=0
        while data.shape[i]!=0:
            i+=1
        data=data.sum(axis=i)
        print data.shape
    if 1 in data.shape:
        print "one dimension"
        i=0
        while data.shape[i]!=1:
            i+=1
        data=data.sum(axis=i)
        print data.shape
    if(len(np.shape(data))==3):
        threedee=True
        print "3D"
    else:
        threedee=False
        print "2D"
    return data,threedee

def KDE_2D(structs,sig):

    import scipy.stats

    ns = len(structs[:,0])

    g = np.zeros((2,ns))

    g[0,:] = structs[:,1]
    g[1,:] = structs[:,2]

    kde = scipy.stats.gaussian_kde(g,bw_method="scott")
    #kde2 = scipy.stats.gaussian_kde(g,bw_method="silverman")
    #kde3 = scipy.stats.gaussian_kde(g,bw_method=sig)

    x,y = np.mgrid[-1.0:1.0:1000j,-1.0:1.0:1000j]
    p = np.vstack([x.ravel(),y.ravel()])

    k1 = np.reshape(kde(p).T,x.shape)
    #k2 = np.reshape(kde2(p).T,x.shape)
    #k3 = np.reshape(kde3(p).T,x.shape)

    k1 = k1.T
    #k2 = k2.T
    #k3 = k3.T

    return k1#,k2,k3


def ReadParameters(param_file):

    type_of_var = {"filename":"str",
                   "root":"str",
                   "dims":"int",
                   "mask_xmin":"int",
                   "mask_xmax":"int",
                   "mask_ymin":"int",
                   "mask_ymax":"int",
                   "min_value":"float",
                   "min_delta":"float",
                   "min_npix":"int",
                   "pixel_pc":"float",
                   "pixel_vel":"float"}

    param = {"filename":"holder",
             "root":"holder",
             "dims":2,
             "mask_xmin":1,
             "mask_xmax":1,
             "mask_ymin":1,
             "mask_ymax":1,
             "min_value":1.0,
             "min_delta":1.0,
             "min_npix":1,
             "pixel_pc":1.0,
             "pixel_vel":1.0,}
    

    with open(param_file) as f:

        for line in f:

	    words = line.split()

            try:
            
                var = type_of_var[words[0]]
                
                if(var=="str"):
                    param[words[0]]=words[1]
                elif(var=="int"):
                    param[words[0]]=int(words[1])                    
                elif(var=="float"):
                    param[words[0]]=float(words[1])
                else:

                    print "The variable is neither a string, float or integer. I don't know how to deal with this"

            except KeyError:

                print "There is no such parameter. Add it to the type_of_var and param dictionaries"
	    

    f.close()

    return param


    

    
def filament_grid(xpix, ypix, nfils, fil_length, fil_width,seed=0):
    """
    nfils randomly oriented filaments of length fil_length and 
    aspect ratio fil_aspecton an xpix by ypix grid. 
    Filaments have uniform density.
    """
    if seed!=0:
        np.random.seed(seed)
    grid=np.zeros((xpix,ypix))
    angles = np.random.rand(nfils)*np.pi
    xpos = np.random.rand(nfils)*xpix
    ypos = np.random.rand(nfils)*ypix
    
    for i in range(nfils):
        grad=np.tan(angles[i])
        for y in range(ypix):
            for x in range(xpix):
                if ((x-xpos[i])**2 + (y-ypos[i])**2) < (0.5*fil_length)**2:
                    #within distance
                    if (y-ypos[i])-(grad*(x-xpos[i]))<0.0001:
                        grid[x,y]+=1
    return grid

    


def m11(grid,com):
    ny=grid.shape[1]
    dr11=0

    ay = np.linspace(0,ny-1,ny)

    dr11 = np.sum(np.sum(grid,axis=0)*(ay-com[1])**2)

    return dr11
    


def m22(grid,com):
    nx=grid.shape[0]
    dr22=0

    ax = np.linspace(0,nx-1,nx)

    dr22 = np.sum(np.sum(grid,axis=1)*(ax-com[0])**2)

    return dr22
    



def m12(grid,com):
    nx=grid.shape[0]
    ny=grid.shape[1]
    dr12=0

    ax = np.linspace(0,nx-1,nx) - com[0]
    ay = np.linspace(0,ny-1,ny) - com[1]

    ax2 = np.array([ax,]*ny).T
    ay2 = np.array([ay,]*nx)

    dr12 = np.sum(grid*ay2*ax2)

    return dr12
    
def theta(M11,M22,M12):

    th=0.5*np.arctan2((2*M12),(M11-M22))

    return th
    

def i12(grid,com):
    M11=m11(grid,com)
    M22=m22(grid,com)
    M12=m12(grid,com)
    tt=theta(M11,M22,M12)

    i1=(M11*np.cos(tt)*np.cos(tt)) + \
        (M22*np.sin(tt)*np.sin(tt)) + \
        (M12*np.sin(2*tt))

    i2=(M11*np.sin(tt)*np.sin(tt)) + \
        (M22*np.cos(tt)*np.cos(tt)) - \
        (M12*np.sin(2*tt))

    return i1,i2,tt
    

def moments_2d(g):
    com=spcom(g)
    Atot=np.count_nonzero(g)
    Mtot=np.sum(g)
    I1,I2,t1=i12(g,com)
    first=min([I1,I2])
    second=max([I1,I2])
    ma=Atot*Mtot
    II1=(ma - 4*np.pi*first)/(ma + 4*np.pi*first)
    II2=(ma - 4*np.pi*second)/(ma + 4*np.pi*second)

    return II1,II2,com,t1

def moments_3d(g):
    com=spcom(g)
    Atot=np.count_nonzero(g)
    Mtot=np.sum(g)
    I1,I2,t1=i12(g,com)
    first=min([I1,I2])
    second=max([I1,I2])
    ma=Atot*Mtot
    II1=(ma - 4*np.pi*first)/(ma + 4*np.pi*first)
    II2=(ma - 4*np.pi*second)/(ma + 4*np.pi*second)

    return II1,II2,II3,com,t1
    
def plot_interactive_axes(d,text=False,threedee=False):
    
    f,bigax=plot_moments_axes(text=text)
    
    f.set_size_inches(12,6)
    
    f.tight_layout(rect=[0, 0, 0.5, 1])
    dendax = f.add_axes([0.5, 0.05, 0.45, 0.85])
    
    #p = d.plotter()
    
    def onpick(event):
        s = event.artist.get_gid()
        print s
        #p.plot_contour(dendax, structure=int(s), lw=3)
        plot_flat_mask(dendax,d,int(s),threedee)
    
    f.canvas.mpl_connect('pick_event', onpick)
    
    return bigax,dendax,f
    
def plot_moments_axes(text=False):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
       
    plt.rcParams.update({'font.size': 18, 'font.family':'serif','text.usetex':False})
    
    f=plt.figure()
    
    bigax=plt.gca()
    bigax.axis('equal')
    bigax.set_xlim((-1.1,1.1))
    bigax.set_ylim((-1.1,1.1))
    bigax.xaxis.tick_top()
    bigax.grid()
    
    bigax.vlines(0,-1.3,1.3)
    bigax.hlines(0,-1.3,1.3)
    
   
    bigax.add_patch(patches.Rectangle((-1.0, -1.0),1,1,alpha=0.2,color='b'))
    bigax.add_patch(patches.Rectangle((0.0, 0.0),1,1,alpha=0.2,color='y'))
    bigax.add_patch(patches.Rectangle((0.0, -1.0),1,1,alpha=0.2,color='r'))
    #bigax.add_patch(patches.Rectangle((0.0, -0.5),0.5,0.5,alpha=0.1,color='r'))
    
    if text:
        fs=16
        bigax.text(0.1,0.8,'Centrally concentrated',color='y',fontsize=fs)
        bigax.text(-0.9,-0.9,'Shell',color='b',fontsize=fs)
        bigax.text(0.5,-0.9,'Filament',color='r',fontsize=fs)
    
    bigax.set_xlabel(r'$J_{1}$')
    bigax.set_ylabel(r'$J_{2}$')
    bigax.xaxis.set_label_position('top')
    return f, bigax
    

def figure_wsc_GC(figure_number, fits_header):
    """
    From Elizabeth Watkins.
    
    This function creates a figure with galactic co-ordinates
    To add annotations or to limit the size, you use pixel co-ordinates
    and normal matplotlib functions.
    To add or change stuff see: 
    http://wcsaxes.readthedocs.io/en/latest/ticks_labels_grid.html#
    """
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt
    
    wcs = WCS(fits_header)

    fig = plt.figure(figure_number)
    ax = fig.add_subplot(111, projection=wcs)

    lon = ax.coords[0]
    lon = ax.coords['glon']
    lon.set_axislabel( 'Galactic Longitude', fontsize = 15 )
    lon.set_major_formatter('dd:mm') # so only shows degs and min
    lon.set_ticks(color='black', exclude_overlapping=True, size = 6)
    lon.set_minor_frequency(5)
    lon.display_minor_ticks(True)

    lat = ax.coords[1]
    lat = ax.coords['glat']
    lat.set_axislabel( 'Galactic Latitude', fontsize = 15 )
    lat.set_major_formatter('dd:mm') # so only shows degs and min
    lat.set_ticks(color='black', exclude_overlapping=True, size = 6)
    lat.set_minor_frequency(5)
    lat.display_minor_ticks(True)
    
def find_coords(gal_lat,gal_lon):
    """
    Convert galactic latitude and longitude in degrees
    to FK5 coordinates.
    gal_lat = 'XXdXXmXXs'
    gal_lon = 'XXdXXmXXs'
    """
    from astropy.coordinates import SkyCoord
    rcw = SkyCoord(gal_lat,gal_lon, frame='galactic')
    return rcw.fk5
    
def compute_dend(params, filename='',root='', data=[]):
    """
    If dend file already exists, read in.
    Else, compute and save dend.
    """
    from astropy.io import fits
    from astrodendro import Dendrogram
    import os
    
    dendfolder='%s/%s_dend_%1.1e_%1.1e_%1.1e'%(root, 
                                                   filename,
                                                   params[0],
                                                   params[1],
                                                   params[2])
    
    if (os.path.isdir(dendfolder) and os.path.isfile("%s/dendrogram.fits"%dendfolder)):
        print "Reading dend"
        d = Dendrogram.load_from("%s/dendrogram.fits"%dendfolder)
    else:
        print "Building dend"
        if len(data)>0:
            array=data
        elif len(filename)>0:
            array = fits.getdata('%s/%s.fits'%(root,filename))
        else:
            print "Please specify a filename and root or data for 'compute_dend' function."
            return 0
        d = Dendrogram.compute(array, 
                               min_value=params[0],
                               min_delta=params[1],
                               min_npix =params[2])
        d.save_to("%s/dendrogram.fits"%dendfolder)

    return d

def distance_to_fil(x,y):
    """Calculate desitance of point (x,y)
    from diagonal line y=-x"""
    return (x+y)/np.sqrt(2)
    
def testdata2d(n,shape='noisy'):
    # n is gridsize
    grid=np.zeros((n,n))
    if shape=='ellipse-cc':
        a=n/2.
        b=n/4.
        
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(((dx**2)/(a**2)) + ((dy**2)/(b**2)))
                if dr<1:
                    grid[i,j]=1
                if dr<0.5:
                    grid[i,j]=3
                if dr<0.2:
                    grid[i,j]=5
    elif shape=='ellipse':
        a=n/2.
        b=n/4.
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(((dx**2)/(a**2)) + ((dy**2)/(b**2)))
                if dr<1.:
                    grid[i,j]=1
    elif shape=='fil-thin':
        grid[n*0.4:n*0.45,:]=1
    elif shape=='fil-thick':
        grid[n*0.3:n*0.6,:]=1
    elif shape=='fil-thin2':
        grid[:,n*0.4:n*0.45]=1
    elif shape=='fil-thick2':
        grid[:,n*0.3:n*0.6]=1
    elif shape=='disk-cc':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/3.:
                    grid[i,j]=3
                if dr<n/5.:
                    grid[i,j]=5
    elif shape=='disk-cc2':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/3.:
                    grid[i,j]=7
                if dr<n/4.:
                    grid[i,j]=15                   
    elif shape=='ring-thick':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/4:
                    grid[i,j]=0
    elif shape=='ring-thin':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/2.5:
                    grid[i,j]=0
    elif shape=='ring-lumpy':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/4:
                    grid[i,j]=0
        noise=np.random.rand(0.01*n**2).reshape(0.1*n,0.1*n)
        for i in range(n):
            for j in range(n):
                grid[i,j]=np.floor(grid[i,j]*2*noise[i/10,j/10])

    elif shape=='disk':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/3.:
                    grid[i,j]=1 
                
    elif shape=='peak':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/20.:
                    grid[i,j]=1
    elif shape=='dot':
        grid[n/2,n/2]=1
    elif shape=='noisy':
        grid=np.random.rand(n**2).reshape((n,n))
    elif shape=='diag-fil':
        for i in range(n):
            grid[i,i]=1   
    elif shape=='diag-fil2':
        for i in range(n):
            grid[i,-i]=1
    elif shape=='asym_cc':
        grid[n*0.3:n*0.6,:]=1
        grid[n*0.4:n*0.5,:]=10
        grid[n*0.4:n*0.5,n*0.2:n*0.8]=30
        #grid[n*0.3:n*0.6,:]=1
        
    elif shape=='fil-thin_dot':
        grid[n*0.4:n*0.45,:]=1
        grid[n*0.6:n*0.62,n*0.6:n*0.62]=1
        
    elif shape=='half-ring-thick':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/4.:
                    grid[i,j]=0
        grid[:,:int(n/2.)]=0
    elif shape=='half-ring-thin':
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j]=1
                if dr<n/2.5:
                    grid[i,j]=0
        grid[:,:int(n/2.)]=0
        
    return grid

def ellipsoid(n,a,b,c):
    grid=np.zeros((n,n,n))
    cen=int(n/2.)
    for i in range(-cen,cen):
        for j in range(-cen,cen):
            for k in range(-cen,cen):
                if (((i**2)/(a**2)) + ((j**2)/(b**2)) + ((k**2)/(c**2)))<1.:
                    grid[i+cen,j+cen,k+cen]=1
    return grid

def testdata3d(n,shape='sphere'):
    # n is gridsize
    grid=np.zeros((n,n,n))
    if shape=='sphere':
        c=int(n/2.)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        grid[i,j,k]=1
                        
    elif shape=='ellipsoid-pro':
        grid=ellipsoid(n,(n/3.),(n/3.),(n/2.))
    
    elif shape=='ellipsoid-ob':
        grid=ellipsoid(n,(n/2.),(n/2.),(n/3.))
        
    elif shape=='ellipsoid-ob2':
        grid=ellipsoid(n,(n/2.),(n/2.),(n/4.))
        
    elif shape=='ellipsoid-ob3':
        grid=ellipsoid(n,(n/2.),(n/2.),(n/5.))
        
    elif shape=='ellipsoid-ob4':
        grid=ellipsoid(n,(n/2.),(n/2.),(n/10.))
        
    elif shape=='ellipsoid-tri':
        grid=ellipsoid(n,(n/2.),(n/4.),(n/7.))
        
    elif shape=='cc':
        c=int(n/2.)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        grid[i,j,k]=1./(0.01+((i-c)**2 + (j-c)**2 + (k-c)**2))
                        
    elif shape=='shell-thick':
        c=int(n/2.)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        if ((i-c)**2 + (j-c)**2 + (k-c)**2)>((n/4.)**2):
                            grid[i,j,k]=1
                            
    elif shape=='shell-thin':
        c=int(n/2.)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        if ((i-c)**2 + (j-c)**2 + (k-c)**2)>((n/2.1)**2):
                            grid[i,j,k]=1

    elif shape=='shell-thick-noise':
        c=int(n/2.)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        if ((i-c)**2 + (j-c)**2 + (k-c)**2)>((n/4.)**2):
                            grid[i,j,k]=1
        noise=np.random.rand(int(0.001*n**3)).reshape(int(0.1*n),int(0.1*n),int(0.1*n))
        noise=np.random.rand(int(n**3)).reshape(int(n),int(n),int(n))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    #grid[i,j,k]=np.floor(grid[i,j,k]*2*noise[i/10,j/10,k/10])
                    grid[i,j,k]=np.floor(grid[i,j,k]*2*noise[i,j,k])

                            
    elif shape=='shell-thin-noise':
        c=int(n/2.)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        if ((i-c)**2 + (j-c)**2 + (k-c)**2)>((n/2.1)**2):
                            grid[i,j,k]=1
        noise=np.random.rand(int(0.001*n**3)).reshape(int(0.1*n),int(0.1*n),int(0.1*n))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    grid[i,j,k]=np.floor(grid[i,j,k]*2*noise[i/10,j/10,k/10])
        
    elif shape=='shell-half':
        c=int(n/2.)
        for i in range(n/2):
            for j in range(n):
                for k in range(n):
                    if ((i-c)**2 + (j-c)**2 + (k-c)**2)<((n/2.)**2):
                        if ((i-c)**2 + (j-c)**2 + (k-c)**2)>((n/2.1)**2):
                            grid[i,j,k]=1
    
    elif shape=='ring-thin':
        k=0
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j,k]=1
                if dr<n/2.5:
                    grid[i,j,k]=0
                    
    elif shape=='ring-half':
        k=0
        for i in range(n/2):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/2.:
                    grid[i,j,k]=1
                if dr<n/2.5:
                    grid[i,j,k]=0
                            
    elif shape=='fil-thick':
                grid[int(n*0.4):int(n*0.6),int(n*0.4):int(n*0.6),:]=1
                            
    elif shape=='fil-thin':
                grid[int(n*0.4):int(n*0.45),int(n*0.4):int(n*0.45),:]=1
                
    elif shape=='fil-asym':
                grid[int(n*0.4):int(n*0.45),int(n*0.4):int(n*0.6),:]=1
                
    elif shape=='plane-thin':
                grid[int(n/2.),:,:]=1
                
    elif shape=='rect-thin':
                grid[int(n/2.),int(n/4.):int(3*n/4.),:]=1
                
    elif shape=='plane-thick':
                grid[int(n*0.4):int(n*0.6),:,:]=1
                
    elif shape=='flat-disk':
        k=0
        for i in range(n):
            dx=i-(n/2.)
            for j in range(n):
                dy=j-(n/2.)
                dr=np.sqrt(dx**2 + dy**2)
                if dr<n/3.:
                    grid[i,j,k]=1 

    elif shape=='thick-disk':
        for k in range(int(0.4*n),int(0.6*n)):
            for i in range(n):
                dx=i-(n/2.)
                for j in range(n):
                    dy=j-(n/2.)
                    dr=np.sqrt(dx**2 + dy**2)
                    if dr<n/3.:
                        grid[i,j,k]=1 
                    
    return grid

def small_3d_plot(ax,g):
    ax.set_xlim3d(0,g.shape[0])
    ax.set_ylim3d(0,g.shape[1])
    ax.set_zlim3d(0,g.shape[2])

    for i in range(g.shape[0]):
        for j in range(g.shape[1]):
            for k in range(g.shape[2]):
                if (g[i,j,k]>0):
                    ax.scatter([i],[j],[k],c='k',marker='.',alpha=0.1)
                    
def image_moment(image, p, q, r=-1, dims=2):
    if dims==2:
        m = im_2d(image, p, q)
    elif dims==3:
        if r==-1:
            raise ValueError("Need order for 3rd dimension moment")
        m = im_3d(image, p, q, r)
    else:
        raise ValueError("Number of dimensions is invalid")
    return m
        
def im_2d(image,p,q):
    nx, ny = image.shape
    x_ind, y_ind = np.mgrid[:nx, :ny]
    moment=(image * x_ind**p * y_ind**q).sum()
    return moment
        
def im_3d(image,p,q,r):
    nx, ny, nz = image.shape
    x_ind, y_ind, z_ind = np.mgrid[:nx, :ny, :nz]
    moment=(image * x_ind**p * y_ind**q * z_ind**r).sum()
    return moment

def im_com(image):
    #0th moment is summ of pixel values
    m00=image_moment(image,0,0)
    m10=image_moment(image,1,0)
    m01=image_moment(image,0,1)
    
    x_com=m10/m00
    y_com=m01/m00
    
    return [x_com,y_com]
    
def plot3dax():    
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure()
    fig.set_size_inches(w=10,h=10)
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlim((-1.1,1.1))
    ax.set_ylim((-1.1,1.1))
    ax.set_zlim((-1.1,1.1))
    
    ax.set_xlabel('$J_{1}$')
    ax.set_ylabel('$J_{2}$')
    ax.set_zlabel('$J_{3}$')
    
    ax.plot([0,0],[0,0],[-1.3,1.3],'k-')
    ax.plot([0,0],[-1.3,1.3],[0,0],'k-')
    ax.plot([-1.3,1.3],[0,0],[0,0],'k-')
    
    ax.axis('equal')
    
    return ax