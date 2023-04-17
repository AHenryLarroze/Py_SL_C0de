from spharm import GaussQuad
import numpy as np
import math
from spharm import sphericalobject
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import stripy
import cartopy
import cartopy.crs as ccrs

from scipy import io

def sph2cart(az, el, r):
    rsin_theta = r * np.sin(el)
    x = rsin_theta * np.cos(az)
    y = rsin_theta * np.sin(az)
    z = r * np.cos(el)
    return x, y, z

class GRID(object):
        """
        A class used to represent the Gaussian Grid

        ...

        Attributes
        ----------
            nlons : int
                number of longitude in the grid
            nlats : int
                number of latitude in the grid
            lats : np.array (maxdeg, maxdeg x2)
                latitude as a meshgrid
            elons : np.array (maxdeg, maxdeg x2)
                longitude as a meshgrid
            colats : np.array (maxdeg, maxdeg x2)
                colatitude as a meshgrid
            time_step : np.array(nb_time_step)
                value of the different time step used in the model.

        Methods
        -------
            add_time(time_step)
            interp_on(grd,lon,lat)
            create_ICE()
            create_SED()
            create_TOPO()
            create_LOVE()
            create_SL()
                     

        """
        def __init__(self,maxdeg):
            """
            Parameters
            ----------
                maxdeg : int
                    The maximum degree of the spherical harmonics used for the creation of the grid
            """
            self.nlons = (maxdeg) * 2 
            self.nlats = (maxdeg)
            x,w=GaussQuad(maxdeg)
            x_GL = np.arccos(x)*180/math.pi - 90
            lon_GL = np.linspace(0,360,2*maxdeg+1)
            lon_GL = lon_GL[:-1]
            self.lats=x_GL
            self.elons=lon_GL
            self.colats = 90 - self.lats
            self.elons,self.lats=np.meshgrid(self.elons,self.lats)
        
        # def add_time(self,time_step):
        #     '''
            
        #     This function add time field to the object.
        
        # Parameters :  
        #     time_step (np.array): a vector array of length number of time step.
             
        # See the documentation of the cited class object for more information on different parameters used in the function.
        
        # Returns :
            
        # Added fields : 
        #     time_step (np.array): an array of each time step of length number of time step.
        
        #     '''
        #     self.time_step=time_step # à ajouter dans ice grid definition
            
        def interp_on(self,grd,lon,lat,smoothing=False):
            '''
            This function interpolate a grid of data on the grid of the model.
        
            Parameters :  
                grd (np.array) : The grid data
                lon (np.array) : The longitude of the data
                lat (np.array) : The latitude of the data          

            See the documentation of the cited class object for more information on different parameters used in the function.
                
            Added fields : 
                grd (np.array) : The interpolated grid with the shape of the gausse legendre grid.
        
            '''
            lon,lat=np.meshgrid(lon,lat)
            vertices_lat=np.radians(lat.flatten())
            vertices_lon=np.radians(lon.flatten())
            spherical_triangulation = stripy.sTriangulation(lons=vertices_lon, lats=vertices_lat,refinement_levels=0)
            if smoothing :
                grd,dds,err=spherical_triangulation.smoothing(grd.flatten(),np.ones_like(grd.flatten()),0.1,0.1,0.01)
            vertices_lats=np.radians(self.lats.ravel())
            vertices_elons=np.radians(self.elons.ravel())
            grd,err=spherical_triangulation.interpolate_nearest(vertices_elons,vertices_lats,data=grd.flatten())
            grd[np.isnan(grd)]=0
            grd=grd.reshape(self.elons.shape)
            return grd
        
        def smooth_on(self,grd,lon,lat):
            lon,lat=np.meshgrid(lon,lat)
            vertices_lat=np.radians(lat.ravel())
            vertices_lon=np.radians(lon.ravel())
            spherical_triangulation = stripy.sTriangulation(lons=vertices_lon, lats=vertices_lat,refinement_levels=0)
            grd,dds,err=spherical_triangulation.smoothing(grd.flatten(),np.ones_like(grd.flatten()),10,0.1,0.01)
            vertices_lats=np.radians(self.lats.ravel())
            vertices_elons=np.radians(self.elons.ravel())
            grd,err=spherical_triangulation.interpolate_nearest(vertices_elons,vertices_lats,data=grd.flatten())
            grd[np.isnan(grd)]=0
            grd=grd.reshape(self.elons.shape)
            return grd

        def disk(self,lat,lon,radius,high,tx=1):
            '''
            A method used to generate a disk. This function have been used for benchmarking of the code. 

            Attribute
            ---------
                lat : nparray
                    Array contaning the latitudinal coordinate of the center of the disk.
                lon : nparray
                    Array containing the longitudinal coordinate of the center of the disk.
                radius : double
                    The radius of the disk in degree
                high : double
                    the thickness of the disk over the considered area.
        
            '''

            grd=np.zeros((tx,self.lats.size,self.elons.size))
            lon_g,lat_g=np.meshgrid(self.elons,self.lats)
            grd[:,((lon-lon_g)**2+(lat-lat_g)**2)<radius]=high
            grd[:,:2,:]=grd[:,:2,:]*0
            return grd

        def zeros(self,tx=1):
            '''
            A method used to generate a zero array with the caracteristics of the grid. 
            '''
            return np.zeros((tx,self.lats.size,self.elons.size))


class TIME_GRID(GRID,sphericalobject):
    """
    This class is used to manage the mass grids. These grids have a time dimenssion this way wa can manage the time variation of the mass. We can define the mass grid by it's mass directly or by coupling a height with a density. If needed you can load spherical harmonics coefficient.

    Attributes
    ----------
        grid_time_step : np.array
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        height_time_grid : np.array
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coedd : np.array
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step

    """

    def __init__(self,time_step_number,maxdeg,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=0):
        super().__init__(maxdeg)
        self.isgrd=False
        self.iscoeff=False

        self.saved=np.array([]) # initialize the save of the spherical harmonic object

        self.rho=rho

        self.height_time_grid=np.zeros((time_step_number,self.nlats,self.nlons))
        self.mass_time_grid=np.zeros((time_step_number,self.nlats,self.nlons))
        self.height_time_coeff=np.zeros((time_step_number,int(maxdeg*(maxdeg+1)/2)))
        self.mass_time_coeff=np.zeros((time_step_number,int(maxdeg*(maxdeg+1)/2)))

        if height_time_grid!=None:
            self.height_time_grid=height_time_grid
            self.mass_time_grid=height_time_grid*rho
            self.grd_0=self.mass_time_grid[0,:,:]
        elif mass_time_grid!=None:
            self.mass_time_grid=mass_time_grid
            self.grd_0=self.mass_time_grid[0,:,:]
        elif height_time_coeff!=None:
            self.height_time_coeff=height_time_coeff
            self.mass_time_coeff=height_time_coeff*rho
            self.coeff_0=self.mass_time_coeff[0,:,:]
        elif mass_time_coeff!=None:
            self.mass_time_coeff=mass_time_coeff
            self.coeff_0=self.mass_time_coeff[0,:,:]

    def interp_on_time(self,grid_to_interp,grid_time_step,model_time_step,interp_type='Thickness_divide',backend='False'):
        """
        This function is used for interpolation upon time and space it call the interpolation function of the GRID parameter. 

        Attributes
        ----------
            time_step : np.array
                this array define the time step on wich the time grid must be interpolated.
            interp_type : str
                Select the type of interpolation you want. 
                    - Thickness_divide : Simply divide the thickness into the different time step. This methods is usefull to preserve the total thickness.
        """

        if interp_type=="Thickness_divide":
            grid_time_step=grid_time_step[::-1]
            model_time_step=model_time_step[::-1]
            grid_to_interp=grid_to_interp[::-1,:,:]
            d_grid_time_step=np.diff(grid_time_step)
            d_grid_to_interp=grid_to_interp/d_grid_time_step[:,np.newaxis,np.newaxis]
            Merge_time_step=np.unique(np.concatenate((grid_time_step,model_time_step)), return_index=True)
            Merge_time_step=(Merge_time_step[0],Merge_time_step[1]+1)
            Merge_time_step[1][Merge_time_step[1]>len(grid_time_step)]=0
            count=np.diff(np.nonzero(Merge_time_step[1])).squeeze()
            to_merge=np.searchsorted(Merge_time_step[0],grid_time_step[np.isin(grid_time_step,model_time_step,invert=True)].squeeze())
            self.height_time_grid=np.array(d_grid_to_interp.repeat(count,axis=0)*np.diff(Merge_time_step[0])[:,np.newaxis,np.newaxis])
            if not(to_merge.dtype==np.int64):
                for i in range(len(to_merge)):
                    if backend :
                        print('time slicing : ' + str(i))
                    self.height_time_grid[to_merge[i]-1,:,:]+=self.height_time_grid[to_merge[i],:,:]
                    self.height_time_grid=np.delete(self.height_time_grid,to_merge[i],0)
                    to_merge-=1
            self.height_time_grid=self.height_time_grid[::-1,:,:]

    def interp_on_time_and_space(self,grid_to_interp,grid_time_step,model_time_step,time_step_number,grid_lon,grid_lat,interp_type='Thickness_divide',backend=False):
        """
        This function is used for interpolation upon time and space it call the interpolation function of the GRID parameter. 

        Attributes
        ----------
            time_step : np.array
                this array define the time step on wich the time grid must be interpolated.
            interp_type : str
                Select the type of interpolation you want. 
                    - Thickness_divide : Simply divide the thickness into the different time step. This methods is usefull to preserve the total thickness.
        """

        #if len(grid_time_step)<len(model_time_step) :
        time_grid_pre_interp=np.zeros((len(grid_time_step)-1,self.elons.shape[0],self.lats.shape[1]))
        for i in range(len(grid_time_step)-1):
            if backend :
                print('interpolation number : ' + str(i))
            time_grid_pre_interp[i,:,:]=self.interp_on(grid_to_interp[i,:,:],grid_lon,grid_lat)
        self.interp_on_time(time_grid_pre_interp,grid_time_step,model_time_step,interp_type,backend)
    
    def grid_from_step(self,t_it):
        """
        This method is used to get the value of the grid at the defined time step
        
        Attributes
        ----------
            t_it : double
                This is the value of the time step iteration on wich you are trying to retreave the grid. It must inside the time_step interpolation you have used during the initialisation of the time grid. 
        """
        if self.mass_time_grid!=None :
            self.grd=self.mass_time_grid[t_it,:,:]
        elif self.height_time_grid!=None:
            self.grd=self.height_time_grid[t_it,:,:]
        self.isgrd=True
        return self

    def coeff_from_step(self,t_it):
        """
        This method is used to get the value of the coefficient at the requested time iteration.
        
        Attributes
        ----------
            t_it : double
                This is the value of the time step iteration on wich you are trying to retreave the coefficient. It must be inside the time_step interpolation you have used during the initialisation of the time grid. 
        """
        if self.mass_time_coeff!=None :
            self.coeff=self.mass_time_coeff[t_it,:,:]
        elif self.height_time_coeff!=None:
            self.coeff=self.height_time_coeff[t_it,:,:]
        self.iscoeff=True
        return self
    
    def zeros_time(self,time_step_number):
        """
        This method is used to define a grid over time with only 0 value
        """
        self.height_time_grid=self.zeros(time_step_number)
    
    def disk_time(self,time_step_number,lat,lon,radius,high):
        """
        This method is used to define a grid over time with a disk defined with it's center coordinate and the height
        """
        self.height_time_grid=self.disk(lat,lon,radius,high,time_step_number)
    
    def update_0(self,time_step_number,type='grd'):
        """
        This function is used to save the first iterration of the model.
        """
        if type=='grd':
            self.grd_0=self.mass_time_grid[time_step_number+1,:,:]
        elif type=='coeff':
            self.coeff_0=self.mass_time_coeff[time_step_number+1,:,:]

    def plot_step_on_sphere(self,time_step,cmap=cm.inferno,vmin=None,vmax=None,clip=False,inverse=False):
        colormap=cmap
        if vmin==None and vmax==None :
            normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
        elif vmax==None and not(vmin==None):
            normaliser = mpl.colors.Normalize(vmin=vmin, vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
        elif vmin==None and not(vmax==None):
            normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=vmax,clip=clip)
        if inverse:
            normaliser.inverse()
        
        u=self.elons/360*2*math.pi
        v=(self.lats+90)/180*math.pi
        x,y,z=sph2cart(u.flatten(),v.flatten(),np.ones((u.flatten().shape)))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(np.reshape(x,(self.lats.shape[0],self.elons.shape[1])),np.reshape(y,(self.lats.shape[0],self.elons.shape[1])),np.reshape(z,(self.lats.shape[0],self.elons.shape[1])),facecolors=colormap(normaliser(self.height_time_grid[time_step,:,:])),cmap=colormap)
        ax.set_aspect('equal')

    def scatter_step_on_sphere(self,time_step,cmap=cm.inferno,vmin=None,vmax=None,clip=False,inverse=False,marker='.',s=0.2):
        colormap=cmap
        if vmin==None and vmax==None :
            normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
        elif vmax==None and not(vmin==None):
            normaliser = mpl.colors.Normalize(vmin=vmin, vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
        elif vmin==None and not(vmax==None):
            normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=vmax,clip=clip)
        if inverse:
            normaliser.inverse()

        u=self.elons/360*2*math.pi
        v=(self.lats+90)/180*math.pi
        x,y,z=sph2cart(u.flatten(),v.flatten(),np.ones((u.flatten().shape)))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x,y,z,marker=marker,c=colormap(normaliser(self.height_time_grid[time_step,:,:].flatten())),s=s)
        ax.set_aspect('equal')

    def plot_step_on_Mollweide_projection(self,time_step,cmap=cm.get_cmap('RdYlBu_r', 10),alpha_ocean=0.5,coast_line_width=0.5):
        fig = plt.figure(figsize=(12, 12), facecolor="none")
        ax  = plt.subplot(111, projection=ccrs.Mollweide())
        ax.set_global()
        colormap = cmap
        m = ax.imshow(self.height_time_grid[time_step,:,:], origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, cmap=colormap, interpolation="gaussian")
        plt.colorbar(mappable=m, orientation="horizontal", shrink=0.5)
        ax.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
        ax.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)

    def plot_step(self,time_step):
        projection2 = ccrs.Mollweide()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.pcolor(self.elons,self.lats,self.height_time_grid[time_step,:,:])


class SEDIMENT_TIME_GRID(TIME_GRID):
    """
    A class used to represent the sediment grid

    Attributes
    ----------
        rho : float
           Density value of the sediment as a constant
        sed : np.array (maxdeg, maxdeg x2)
           Gaussian grid with the sediment thickness
    """
    def __init__(self) : 
        """
        Parameters
        ----------
        """
        # preset the grid and coefficient to none. The sediment density is set to 2300 
        super().__init__()
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
    
    def load(self,way,time_step,type="height"):
        '''
        function to load the sediment grid. The actual version simply create a matrix of 0.
    
        Parameters : 
            grid (object): The output of the class GRID. 
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            sed (np.array): the sediment thickness threw time the shape is number of time step, maxdeg, maxdegx2
                            The definition here is the one of Mitrovica code. 
                            So the thickness of sediment at time is the thickness deposited from that time till present.

        '''
        #self.sed=np.zeros((len(model_p.grid.time_step),model_p.grid.nlats,model_p.grid.nlons))
        data=io.loadmat(way)
        lon=data['lon']
        lat=data['lat']
        grd=data['Sed_map_fin']
        if type=="height":
            self.height_time_grid=self.interp_on_time_and_space(lat,lon,grd,time_step) # NOT ITME step enterd
            self.height_time_grid=self.height_time_grid*self.rho
            self.grd_0=self.mass_time_grid[0,:,:]
        elif type=="mass":
            self.mass_time_grid=self.interp_on_time_and_space(lat,lon,grd,time_step) # NOT ITME step enterd
            self.grd_0=self.mass_time_grid[0,:,:]
        return self

class ICE_TIME_GRID(TIME_GRID):
    """
    A class used to represent the ice grid

    ...

    Attributes
    ----------
        rho : float
           Density value of the ice
        ice : np.array (nb_time_step, maxdeg, maxdeg x2)
           Ice thickness Gaussian grid
        ice_corrected : np.array (np_time_step, maxdeg, maxdeg x2)
           Ice thickness Gaussian grid corrected of the crusted ice
        
    Methods
    -------
        load(grid)
           load the ice from ice_6g model  !!! need to be modified !!!
        quickload(grid)
           load the ice from ice_6g model and is parallelized  !!! need to be modified !!!
        calc_del_ice(t_it)
           update the grd field with the variation of ice thickness
        
    """
    def __init__(self) : 
        """
        Parameters
        ----------
        """
        # initialize to false the grid and coefficient (no grid or coefficient pre loaded). The ice volumetric masse is set to 916.7 kg/m3.
        super().__init__()
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        
    def load(self,time_step,time_step_number):
        '''
        
        Function to load the ice data from the ice6g map. 
        I'll have to create a function to simply interpolate a time grid for ice.
        This way anyone can load a data without having to modify or create a function. 
        
        Parameters :  
            grid (object): Output of the class GRID.
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            ice (np.array): The ice thickness threw time matrix interpolated on the Gaussian grid of GRID class of shape time step, maxdeg, maxdeg x 2.
            ice_corrected (np.array): an array full of zero with the same shape then the ice field. 
        ''' 
        data=io.loadmat('ice6g_data') #load the file, can be modified to load an other file.
        ice_in=data['ice6g']
        ice_time=data['ice_time']
        ice_time=ice_time.squeeze()
        # precreate the ice matrix. 
        # same for the ice_corrected matrix. 
        #create an enlarged ice_lat and ice_long matrix
        ice_lat =data['ice_lat']
        ice_lon =data['ice_long']
        # do the same thing with ice_extended
        # for each time step interpolate the ice grid on the Gaussian grid defined in GRID class
        self.interp_on_time_and_space(ice_lat,ice_lon,ice_in,time_step)
        return self 

    def ice_correction(self,time_step_number,topo,oc):

        for t_it in range(len(time_step_number)):
            topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]-self.ice_corrected[t_it,:,:]+self.ice[t_it,:,:]
        for t_it in range(time_step_number): # to parallelize
            check1 = OCEAN_TIME_GRID().evaluate_ocean(-topo.height_time_grid[t_it,:,:]+self.ice[t_it,:,:]) # generate the ocean function for ice-topo
            check2 = OCEAN_TIME_GRID().evaluate_ocean(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:]).grd*(OCEAN_TIME_GRID().evaluate_ocean(-self.ice[t_it,:,:]*self.rho-(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:])*oc.rho).grd)
            self.ice_corrected[t_it,:,:] = check1.grd*self.ice[t_it,:,:]+check2*self.ice[t_it,:,:] # add the two part of ice over the check1 nd check2 positive area.
        for t_it in range(time_step_number):    
            topo.height_time_grid[t_it,:,:]=self.height_time_grid[t_it,:,:]+self.ice_corrected[t_it,:,:]-self.ice[t_it,:,:]



class OCEAN_TIME_GRID(TIME_GRID):
    """
    A class used to represent the ocean grid

    ...

    Attributes
    ----------
        rho : float
           Density value of the water
        oc_0_grd : np.array (maxdeg, maxdeg x2)
           Gaussian grid with 0 or 1 copy of the ocean function at time 0 normaly. 
        oc_0_coeff : np.array (maxdeg, maxdeg) !!!! A vérifier !!!! 
           spherical harmonic coefficient copy of the ocean function at time 0 normaly
           
    Methods
    -------
        update_oc_0()
            set the oc_0_grd and oc_0_coeff attribute as a copy of the grd and coeff
        evaluate_ocean(topo)
            evaluate the ocean function based on the Gaussian grid of the topography topo

    """
    def __init__(self) :
        """
        Parameters
        ----------
        """
        # initialize the ocean with no grid and no coefficient. The saved is also initialized. The volumetric mass of water is set to 1000 
        super().__init__()
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        
        
    def evaluate_ocean(self,topo) :
        '''
        
        Evaluate the ocean function using the topography. It create a 0-1 matrix wich is 1 where topo<0 and 0 where topo>0.  
    
        Parameters : 
            topo (np.array): topo gaussian grid of size maxdeg, maxdegx2
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            self.grd (np.array): a matrix of 0 and 1 with shape maxdeg, maxeg x 2.

        '''
        # use sign function to optimize the conversion from positive negative value to boolean value.
        out = -0.5*np.sign(topo)+0.5
        out = 0.5*np.sign(out-0.6)+0.5
        # set the grd to the output of the previous computation.
        self.grd=out
        self.isgrd=True
        return self

class TOPOGRAPHIC_TIME_GRID(TIME_GRID):
    """
    A class used to represent the topographic grid

    Attributes
    ----------
        topo_pres : np.array (maxdeg, maxdeg x2)
            The present topography gaussian grid
        topo_initial : np.array (maxdeg, maxdeg x2)
            The topography at the begining of the modelization (maxdeg, maxdeg x2)
        topo : np.array (n_time_step, maxdeg, maxdeg x2)
            The topography threw time
        topo_0 : np.array (maxdeg, maxdeg x2)
            Gaussian grid of the topography at the time of run update_topo_0. Used to save the initial topography after it's modification by the load. 

    Methods
    -------
        load(grid,,ice,sed)
            Load the topographic data from the file topo_SL. Create the topography variation by including the ice thickness and the sediment thickness
        update_topo_0()
            update the topo_0 attribute.
        ice_correction(ice,grid,oc)
            correct the ice and the topography from the land attached ice variation.          

    """
    def __init__(self) : 
        """
        Parameters
        ----------
        """
        super().__init__()
        # initialize coefficient and grid to 0 because no grid or coefficient has been created
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        
    def load(self,SL_model,way):
        '''
        Load the topography and adapt it with the ice and sediment thickness.
        
        Parameters :  
            grid (object): output of the class GRID
            ice (object): output of the class spherical_ice
            sed (object): output of the class spherical_sed
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            topo_pres (np.array): matrix of the present topography with shape maxdeg, maxdeg x 2. This grid is 
                                  defined to compare with the topography at the end of the simulation.
            topo_initial (np.array): matrix of the topography at the begining of the modelisation, with shape maxdeg, maxdeg x 2.
            topo (np.array):matrix of the topography threw time with shape time step, maxdeg, maxdeg x 2.
            grd (np.array): graussian grid set to the topography at the beginning of the simulation, maxdeg, maxdeg x 2.
        
        '''
        data=io.loadmat(way) # load the topographic data from topo_SL
        self.topo_pres=SL_model.grid.interp_on(data['topo_bed'],np.transpose(data['lon_topo'][0]),np.transpose(data['lat_topo'][0]))+SL_model.ice.ice[-1,:,:] 
        self.topo_initial=self.topo_pres - SL_model.ice.ice[-1,:,:]-SL_model.sed.sed[0,:,:]+SL_model.ice.ice[0,:,:]
        for i in range(1,len(SL_model.grid.time_step)):
            self.height_time_grid[i,:,:]=self.topo_pres - SL_model.ice.ice[0,:,:]-SL_model.sed.sed[i,:,:]+SL_model.ice.ice[i,:,:]
        return self
    
    def update_topo_0(self):
        '''
        update the topography 0 wich is the topography at the beginning of the simulation. This function is used to update
        the initial topography after each simulation to restart a loop over it. This way there is a convergence and we take in account the GIA.
        
        Parameters :  
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            topo_0 (np.array): A matrix of the topography at the bigining of the simulation of shape maxdeg, maxdeg x 2.  
        '''
        self.topo_0=self.topo[0,:,:].copy() # update the actual topography 
        return self