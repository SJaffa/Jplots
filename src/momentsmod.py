# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 10:45:04 2017

@author: sjaffa
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass as spcom
    
def make_grid_2d(l):
    """
    Make 2D projected density grid from astrodendro structure.

    Parameters
    ----------
    l : astrodendro Structure instance

    Returns
    ----------
    grid : array_like
        2D grid of projected density along axis 0
    n : int
        number of pixels in the structure
    """
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



def plot_flat_mask(ax,d,idx,threedee,col='w'):
    """
    Plot mask of Structure over current axes.

    Parameters
    ----------
    ax : matplotlib axes instance
        Axes on wihch to plot
    d : astrodendro Dendrogram instance
        Dendrogram
    idx : int
        Index of structure to plot
    threedee : boolean
        True indicates 3D, False indicates 2D
    col : string, optional
        Matplotlib colour code to use when plotting mask contour.
        Default is white: 'w'
    
    Returns
    ----------
    
    """
    if threedee:
        mask_3d=d[idx].get_mask()
        mask_2d=np.sum(mask_3d,axis=0)
    else:
        mask_2d=d[idx].get_mask()
    ax.contour(mask_2d,colors=col,levels=[0])

def single_structure_evolution(d,struct_id,nslices=20):
    """
    Take a single structure form the dendrogram and analyse it at several
    logarithmically spaced intensity levels, plotting the resultant J values.

    Parameters
    ----------
    d : astrodendro Dendrogram instance
        Dendrogram from which toextract structure
    struct_id : int
        Index (id) od structure
    nslices : int, optional
        Number of logarithmically spaced intervals to analyse.
        Default = 20            

    Returns
    ----------
    j1s : array
        J1 values
    j2s : array
        J2 values
    I : float
        maximum intensity analysed
    fig : matplotlib Figure instance
        Plot of changing J values
    """
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
        #print(I)
        
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

def check_dimensions(data):
    """
    Check dimensionality of an array and collapse if 1 or 0 dimensioned
    Parameters
    ----------
    data :array
        Array to check
    
    Returns
    ----------
    data: array
        Array collapsed along any 0 or 1 dimensional axes
    threedee: boolean
        Flag is True if 3D, False if 2D
    """
    print(data.shape)
    if 0 in data.shape:
        print("zero dimension")
        i=0
        while data.shape[i]!=0:
            i+=1
        data=data.sum(axis=i)
        print(data.shape)
    if 1 in data.shape:
        print("one dimension")
        i=0
        while data.shape[i]!=1:
            i+=1
        data=data.sum(axis=i)
        print(data.shape)
    if(len(np.shape(data))==3):
        threedee=True
        print("3D")
    else:
        threedee=False
        print("2D")
    return data,threedee



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
    """
    Calculate 2D J values
    
    """
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
    """
    Compute 3D J values    
    """
    #Find mass, volume, centre of mass
    com=spcom(g)
    Vtot=np.count_nonzero(g) #volume
    Mtot=np.sum(g) #mass

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

    return J1,J2,J3
    
def plot_interactive_axes(d,text=False,threedee=False):
    """
     Create interactive axes   
    """
    f=plt.figure()
    jax=plt.gca()
    
    jax=plot_moments_axes(jax,text=text)
    
    f.set_size_inches(12,6)
    
    f.tight_layout(rect=[0, 0, 0.5, 1])
    dendax = f.add_axes([0.5, 0.05, 0.45, 0.85])
    
    #p = d.plotter()
    
    def onpick(event):
        s = event.artist.get_gid()
        print(s)
        #p.plot_contour(dendax, structure=int(s), lw=3)
        plot_flat_mask(dendax,d,int(s),threedee)
    
    f.canvas.mpl_connect('pick_event', onpick)
    
    return jax,dendax,f
    
def plot_moments_axes(ax,text=False):
    """
    Create Jplots blank axes with labels    
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
       
    plt.rcParams.update({'font.size': 18, 'font.family':'serif','text.usetex':False})

    ax.axis('equal')
    ax.set_xlim((-1.1,1.1))
    ax.set_ylim((-1.1,1.1))
    ax.xaxis.tick_top()
    ax.grid()
    
    ax.vlines(0,-1.3,1.3)
    ax.hlines(0,-1.3,1.3)
    
   
    ax.add_patch(patches.Rectangle((-1.0, -1.0),1,1,alpha=0.2,color='b'))
    ax.add_patch(patches.Rectangle((0.0, 0.0),1,1,alpha=0.2,color='y'))
    ax.add_patch(patches.Rectangle((0.0, -1.0),1,1,alpha=0.2,color='r'))
    #bigax.add_patch(patches.Rectangle((0.0, -0.5),0.5,0.5,alpha=0.1,color='r'))
    
    if text:
        fs=16
        ax.text(0.1,0.8,'Centrally concentrated',color='y',fontsize=fs)
        ax.text(-0.9,-0.9,'Bubble',color='b',fontsize=fs)
        ax.text(0.5,-0.9,'Filament',color='r',fontsize=fs)
    
    ax.set_xlabel(r'$J_{1}$')
    ax.set_ylabel(r'$J_{2}$')
    ax.xaxis.set_label_position('top')
    return ax

    
def testdata2d(n,shape='noisy'):
    """
    Create various test shapes in 3D    
    """
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
        grid[int(n*0.4):int(n*0.45),:]=1
    elif shape=='fil-thick':
        grid[int(n*0.3):int(n*0.6),:]=1
    elif shape=='fil-thin2':
        grid[:,int(n*0.4):int(n*0.45)]=1
    elif shape=='fil-thick2':
        grid[:,int(n*0.3):int(n*0.6)]=1
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
        noise=np.random.rand(int(0.01*n**2)).reshape(int(0.1*n),int(0.1*n))
        for i in range(n):
            for j in range(n):
                grid[i,j]=np.floor(grid[i,j]*2*noise[int(i/10),int(j/10)])

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
    """
    Make 3D pixelated ellipsoid
    """
    grid=np.zeros((n,n,n))
    cen=int(n/2.)
    for i in range(-cen,cen):
        for j in range(-cen,cen):
            for k in range(-cen,cen):
                if (((i**2)/(a**2)) + ((j**2)/(b**2)) + ((k**2)/(c**2)))<1.:
                    grid[i+cen,j+cen,k+cen]=1
    return grid

def testdata3d(n, shape='sphere'):
    """
    Make simple 3D shapes for testing J3D    
    """
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
                    grid[i,j,k]=np.floor(grid[i,j,k]*2*noise[int(i/10),int(j/10),int(k/10)])
        
    elif shape=='shell-half':
        c=int(n/2.)
        for i in range(int(n/2)):
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
        for i in range(int(n/2)):
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
       
    
def plot3dax():
    """
    Make 3D axes and set limits and labels
    """
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
    
    #ax.axis('equal')
    
    return ax
