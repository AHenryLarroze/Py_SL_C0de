import numpy as np
import math
from .spharm import sphericalobject
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import stripy
from netCDF4 import Dataset
import pyshtools as pysh
from .Load import LOAD

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
        def __init__(self):
            """
            Parameters
            ----------
                maxdeg : int
                    The maximum degree of the spherical harmonics used for the creation of the grid
            """
            self.nlons = (self.maxdeg) * 2 
            self.nlats = (self.maxdeg)
            x,w=pysh.expand.SHGLQ(self.maxdeg-1)
            x_GL = np.arccos(x[::-1])*180/math.pi - 90
            lon_GL = np.linspace(0,360,2*self.maxdeg+1)
            lon_GL = lon_GL[:-1]
            self.lats=x_GL
            self.elons=lon_GL
            self.colats = 90 - self.lats
            #self.elons,self.lats=np.meshgrid(self.elons,self.lats)
        
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
            
        def interp_on(self,grd,lon,lat,smoothing=False,grid_type='global',error=False):
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
            if grid_type=='global':
                lon,lat=np.meshgrid(lon,lat)
                vertices_lat=np.radians(lat.flatten())
                vertices_lon=np.radians(lon.flatten())
                spherical_triangulation = stripy.sTriangulation(lons=vertices_lon, lats=vertices_lat,refinement_levels=0,permute=True)
                if smoothing :
                    grd,dds,err=spherical_triangulation.smoothing(grd.flatten(),np.ones_like(grd.flatten()),0.1,0.1,0.01)
                elons,lats=np.meshgrid(self.elons,self.lats)
                vertices_lats=np.radians(lats.flatten())
                vertices_elons=np.radians(elons.flatten())
                grd,err=spherical_triangulation.interpolate_nearest(vertices_elons,vertices_lats,data=grd.flatten())
                grd[np.isnan(grd)]=0
            elif grid_type=='local':
                # Select the points inside the local studied zone.
                spherical_triangulation = stripy.Triangulation(x=lon, y=lat,permute=True)
                if smoothing :
                    grd,dds,err=spherical_triangulation.smoothing(grd,np.ones_like(grd),0.1,0.1,0.01)
                elons,lats=np.meshgrid(self.elons,self.lats)
                lats=lats.flatten()
                elons=elons.flatten()
                point_in=((elons>lon.min())*(elons<lon.max()))*((lats>lat.min())*(lats<lat.max()))
                elons_in=elons[point_in]
                lats_in=lats[point_in]
                grd_in,err=spherical_triangulation.interpolate_nearest(elons_in,lats_in,grd)
                # I need to complete the grid where the is no data
                grd=np.zeros(elons.shape)
                grd[point_in]=grd_in
            else:
                print('No such grid type, please select one between : regular or samples')
            if error : 
                return grd.reshape((self.nlats,self.nlons)), err.reshape((self.nlats,self.nlons))
            return grd.reshape((self.nlats,self.nlons))
        
        def smooth_on(self,grd,lon,lat):
            lon,lat=np.meshgrid(lon,lat)
            vertices_lat=np.radians(lat.ravel())
            vertices_lon=np.radians(lon.ravel())
            spherical_triangulation = stripy.sTriangulation(lons=vertices_lon, lats=vertices_lat,refinement_levels=0)
            grd,dds,err=spherical_triangulation.smoothing(grd.flatten(),np.ones_like(grd.flatten()),10,0.1,0.01)
            elons,lats=np.meshgrid(self.elons,self.lats)
            vertices_lats=np.radians(lats.ravel())
            vertices_elons=np.radians(elons.ravel())
            grd,err=spherical_triangulation.interpolate_nearest(vertices_elons,vertices_lats,data=grd.flatten())
            grd[np.isnan(grd)]=0
            grd=grd.reshape(elons.shape)
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

    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=0,grid_name='time_grid',from_file=(False,)):
        if not(from_file[0]) :
            self.maxdeg=maxdeg
            super().__init__()
            self.isgrd=False
            self.iscoeff=False


            self.time_grid_name=grid_name

            self.saved=np.array([]) # initialize the save of the spherical harmonic object

            self.rho=rho

            self.time_step=time_step
            self.time_step_number=len(time_step)

            self.height_time_grid=np.zeros((self.time_step_number-1,self.nlats,self.nlons))
            self.mass_time_grid=np.zeros((self.time_step_number-1,self.nlats,self.nlons))
            self.height_time_coeff=np.zeros((self.time_step_number-1,int(maxdeg*(maxdeg+1)/2)))+0j
            self.mass_time_coeff=np.zeros((self.time_step_number-1,int(maxdeg*(maxdeg+1)/2)))+0j

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

        elif from_file[0] :
            ncgrid = Dataset(from_file[1]+'.nc',mode='r',format='NETCDF4_CLASSIC') 
            self.maxdeg=len(ncgrid['lat'][:].data)
            super().__init__()

            self.time_step=ncgrid['time'][:].data
            self.time_step_number=len(self.time_step)

            self.rho=rho

            self.height_time_grid=ncgrid['thickness'][:].data
            self.mass_time_grid=np.zeros((self.time_step_number,self.nlats,self.nlons))
            self.height_time_coeff=np.zeros((self.time_step_number,int(self.maxdeg*(self.maxdeg+1)/2)))+0j
            self.mass_time_coeff=np.zeros((self.time_step_number,int(self.maxdeg*(self.maxdeg+1)/2)))+0j

            self.mass_time_grid=self.height_time_grid*rho

            ncgrid.close()

    def interp_on_time(self,grid_to_interp,grid_time_step,model_time_step,interp_type='Thickness_divide',backend='False',grid_type='regular'):
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
            if grid_type=='global':
                grid_to_interp=grid_to_interp[::-1,:,:]
                d_grid_time_step=np.diff(grid_time_step)
                d_grid_to_interp=grid_to_interp/d_grid_time_step[:,np.newaxis,np.newaxis]

                Merge_time_step=np.unique(np.concatenate((grid_time_step,model_time_step)), return_index=True)
                Merge_time_step=(Merge_time_step[0],Merge_time_step[1]+1)
                Merge_time_step[1][Merge_time_step[1]>len(grid_time_step)]=0

                out=len(model_time_step[model_time_step>grid_time_step.max()])

                count=np.diff(np.nonzero(Merge_time_step[1])).squeeze()
                to_merge=np.searchsorted(Merge_time_step[0],grid_time_step[np.isin(grid_time_step,model_time_step,invert=True)].squeeze())-1
                grd_interpolated=np.zeros((len(Merge_time_step[0])-1,d_grid_to_interp.shape[1],d_grid_to_interp.shape[2]))
                grd_interpolated=np.array(np.concatenate((d_grid_to_interp.repeat(count,axis=0),np.zeros((out,d_grid_to_interp.shape[1],d_grid_to_interp.shape[2]))))*np.diff(Merge_time_step[0])[:,np.newaxis,np.newaxis])

                for i in range(len(to_merge)):
                    if backend :
                        print('time slicing : ' + str(i))
                    grd_interpolated[to_merge[i]-1,:,:]+=grd_interpolated[to_merge[i],:,:]
                    grd_interpolated=np.delete(grd_interpolated,to_merge[i],0)
                    to_merge-=1
                return grd_interpolated[::-1,:,:]
            elif grid_type=='local':
                grid_to_interp=grid_to_interp[::-1,:]
                d_grid_time_step=np.diff(grid_time_step)
                d_grid_to_interp=grid_to_interp/d_grid_time_step[:,np.newaxis]

                Merge_time_step=np.unique(np.concatenate((grid_time_step,model_time_step)), return_index=True)
                Merge_time_step=(Merge_time_step[0],Merge_time_step[1]+1)
                Merge_time_step[1][Merge_time_step[1]>len(grid_time_step)]=0

                out=len(model_time_step[model_time_step>grid_time_step.max()])


                count=np.diff(np.nonzero(Merge_time_step[1])).squeeze()
                to_merge=np.searchsorted(Merge_time_step[0],grid_time_step[np.isin(grid_time_step,model_time_step,invert=True)].squeeze())-1
                grd_interpolated=np.zeros((len(Merge_time_step[0])-1,d_grid_to_interp.shape[1]))
                grd_interpolated=np.array(np.concatenate((d_grid_to_interp.repeat(count,axis=0),np.zeros((out,d_grid_to_interp.shape[1]))))*np.diff(Merge_time_step[0])[:,np.newaxis])

                for i in range(len(to_merge)):
                    if backend :
                        print('time slicing : ' + str(i))
                    grd_interpolated[to_merge[i]-1,:]+=grd_interpolated[to_merge[i],:]
                    grd_interpolated=np.delete(grd_interpolated,to_merge[i],0)
                    to_merge-=1
                return grd_interpolated[::-1,:]
            else :
                print('no such grid type avaiable. try samples or structured')

    def interp_on_time_and_space(self,grid_to_interp,grid_time_step,grid_lon,grid_lat,interp_type='Thickness_divide',backend=False,grid_type='global'):
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
        if len(grid_time_step)<self.time_step_number:
            time_grid_pre_interp=np.zeros((len(grid_time_step)-1,self.nlats,self.nlons))
            for i in range(len(grid_time_step)-1):
                if backend :
                    print('interpolation number : ' + str(i))

                time_grid_pre_interp[i,:,:]=self.interp_on(grid_to_interp[i],grid_lon,grid_lat,grid_type=grid_type)
            grid_type='global'
            self.height_time_grid=self.interp_on_time(time_grid_pre_interp,grid_time_step,self.time_step,interp_type,backend=backend,grid_type=grid_type)
        else :
            grd_to_interpolate=self.interp_on_time(grid_to_interp,grid_time_step,self.time_step,interp_type,backend=backend,grid_type=grid_type)
            for i in range(self.time_step_number-1):
                if backend :
                    print('interpolation number : ' + str(i))
                self.height_time_grid[i,:,:]=self.interp_on(grd_to_interpolate[i,:],grid_lon,grid_lat,grid_type=grid_type)

    
    def grid_from_step(self,t_it):
        """
        This method is used to get the value of the grid at the defined time step
        
        Attributes
        ----------
            t_it : double
                This is the value of the time step iteration on wich you are trying to retreave the grid. It must inside the time_step interpolation you have used during the initialisation of the time grid. 
        """
        # if not(self.mass_time_grid is None) :
        #     self.grd=self.mass_time_grid[t_it,:,:]
        # elif not(self.height_time_grid is None):
        self.grd=self.height_time_grid[t_it,:,:].copy()
        return self

    def coeff_from_step(self,t_it):
        """
        This method is used to get the value of the coefficient at the requested time iteration.
        
        Attributes
        ----------
            t_it : double
                This is the value of the time step iteration on wich you are trying to retreave the coefficient. It must be inside the time_step interpolation you have used during the initialisation of the time grid. 
        """
        # if not(self.mass_time_coeff is None) :
        #     self.coeff=self.mass_time_coeff[t_it,:]
        # elif not(self.height_time_coeff is None):
        self.coeff=self.height_time_coeff[t_it,:].copy()
        return self
    
    def timegrdtotimecoeff(self):
        for i in range(self.time_step_number):
            self.grd=self.height_time_grid[i,:,:]
            self.height_time_coeff[i,:]=self.grdtocoeff().coeff
        return self
    
    def timecoefftotimegrd(self):
        for i in range(self.time_step_number):
            self.coeff=self.height_time_coeff[i,:]
            self.height_time_grid[i,:,:]=self.coefftogrd().grd
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
    
    def update_0(self):
        """
        This function is used to save the first iterration of the model.
        """
        if not(self.height_time_grid is None) :
            self.grd_0=self.height_time_grid[0,:,:].copy()
        if not(self.height_time_coeff is None) :
            self.coeff_0=self.height_time_coeff[0,:].copy()

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
        
        elons,lats=np.meshgrid(self.elons,self.lats)
        
        u=elons/360*2*math.pi
        v=(lats+90)/180*math.pi
        x,y,z=sph2cart(u.flatten(),v.flatten(),np.ones((u.flatten().shape)))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(np.reshape(x,(self.lats.shape[0],self.elons.shape[0])),np.reshape(y,(self.lats.shape[0],self.elons.shape[0])),np.reshape(z,(self.lats.shape[0],self.elons.shape[0])),facecolors=colormap(normaliser(self.height_time_grid[time_step,:,:])),cmap=colormap)
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

        elons,lats=np.meshgrid(self.elons,self.lats)
        
        u=elons/360*2*math.pi
        v=(lats+90)/180*math.pi
        x,y,z=sph2cart(u.flatten(),v.flatten(),np.ones((u.flatten().shape)))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x,y,z,marker=marker,c=colormap(normaliser(self.height_time_grid[time_step,:,:].flatten())),s=s)
        ax.set_aspect('equal')

    def plot_step(self,time_step):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        elons,lats=np.meshgrid(self.elons,self.lats)
        ax.pcolor(elons,lats,self.height_time_grid[time_step,:,:])

    def save(self,save_way=''):
        ncgrid=Dataset(save_way+self.time_grid_name+'.nc','w','NETCDF4_CLASSIC')
        ncgrid.title=self.time_grid_name
        ncgrid.createDimension('lon',self.nlons)
        ncgrid.createDimension('lat',self.nlats)
        ncgrid.createDimension('time',self.time_step_number-1)
        lat=ncgrid.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:]=self.lats
        lon=ncgrid.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        lon[:]=self.elons
        time=ncgrid.createVariable('time', np.float32, ('time',))
        time.units = 'kyr'
        time.long_name = 'time'
        time[:]=self.time_step[:-1]
        thickness=ncgrid.createVariable('thickness',np.float32,('time','lat','lon'))
        thickness.units='m'
        thickness.long_name='layer_thickness'
        thickness[:,:,:]=self.height_time_grid
        ncgrid.close()

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
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=2600,grid_name='time_grid',from_file=(False,)) : 
        """
        Parameters
        ----------
        """
        # preset the grid and coefficient to none. The sediment density is set to 2300 
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file)
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        self.rho=rho

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
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=916.7,grid_name='time_grid',from_file=(False,)) : 
        """
        Parameters
        ----------
        """
        # initialize to false the grid and coefficient (no grid or coefficient pre loaded). The ice volumetric masse is set to 916.7 kg/m3.
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file)
        self.ice=self.height_time_grid
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        self.sdeli_00=0
        self.deli_00_prev=0

    def ice_correction(self,topo,oc):

        for t_it in range(self.time_step_number):
            topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]-self.height_time_grid[t_it,:,:]+self.ice[t_it,:,:]
        for t_it in range(self.time_step_number): # to parallelize
            check1 = OCEAN_TIME_GRID().evaluate_ocean(-topo.height_time_grid[t_it,:,:]+self.ice[t_it,:,:]) # generate the ocean function for ice-topo
            check2 = OCEAN_TIME_GRID().evaluate_ocean(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:]).grd*(OCEAN_TIME_GRID().evaluate_ocean(-self.ice[t_it,:,:]*self.rho-(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:])*oc.rho).grd)
            self.height_time_grid[t_it,:,:] = check1.grd*self.ice[t_it,:,:]+check2*self.ice[t_it,:,:] # add the two part of ice over the check1 nd check2 positive area.
        for t_it in range(self.time_step_number):    
            topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]+self.height_time_grid[t_it,:,:]-self.ice[t_it,:,:]



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
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=1000,grid_name='time_grid',from_file=(False,)) :
        """
        Parameters
        ----------
        """
        # initialize the ocean with no grid and no coefficient. The saved is also initialized. The volumetric mass of water is set to 1000 
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file)
        self.saved=np.array([])
        self.rho=rho
        
    def update_0(self):
        self.grd_0=self.grd.copy()
        
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
        return self
    
    def sea_level_solver(self,load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,conv_lim):
        if conv_it==0 :
            if t_it==0:
                self.height_time_coeff[t_it,:]=self.prev/self.prev[0]*(TO.coeff[0]-TO.prev[0])-TO.coeff-TO.prev
            else :
                self.height_time_coeff[t_it,:]=self.prev/self.prev[0]*(-ice_time_grid.rho/self.rho*ice_time_grid.height_time_coeff[:t_it+1,:].sum(0) + TO.coeff[0]-TO.prev[0])-TO.coeff-TO.prev
        chi=20
        while chi>conv_lim:
            chi=self.sea_level_equation(load,ice_time_grid,sed_time_grid,love_number,TO,t_it)
            conv_it+=1
        return conv_it
    
    def sea_level_equation(self,load,ice_time_grid,sed_time_grid,love_number,TO,t_it):
        if t_it==0 :
            delSLcurl_fl=love_number.E* love_number.T.coeff *(ice_time_grid.height_time_coeff[0,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[0,:]*sed_time_grid.rho+self.height_time_coeff[0,:]*self.rho)
            self.delSLcurl=sphericalobject(coeff=delSLcurl_fl - ice_time_grid.height_time_coeff[0,:]- sed_time_grid.height_time_coeff[0,:]).coefftogrd()
            RO=sphericalobject(grd=self.delSLcurl.grd*self.grd).grdtocoeff()
            self.delPhi_g=np.real(1/self.coeff[0] * (- ice_time_grid.rho/self.rho*ice_time_grid.coeff[0] - RO.coeff[0] + TO.coeff[0]))
            chi = np.abs( (np.sum(np.abs(RO.coeff + self.delPhi_g*self.coeff -  TO.coeff)) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
            self.height_time_coeff[t_it,:]=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff
        else :
            if t_it == 1 : 
                load.calc_viscuous_load(ice_time_grid.height_time_coeff[0,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[0,:]*sed_time_grid.rho+self.height_time_coeff[0,:]*self.rho,love_number.beta_l,0)
            else : 
                load.calc_viscuous_load(ice_time_grid.height_time_coeff[:t_it,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[:t_it,:]*sed_time_grid.rho+self.height_time_coeff[:t_it,:]*self.rho,love_number.beta_l,t_it)
            delSLcurl_fl=love_number.E* love_number.T.coeff *(ice_time_grid.height_time_coeff[:t_it+1,:].sum(0)*ice_time_grid.rho+sed_time_grid.height_time_coeff[:t_it+1,:].sum(0)*sed_time_grid.rho+self.height_time_coeff[:t_it+1,:].sum(0)*self.rho)+love_number.T.coeff*load.V_lm.coeff
            self.delSLcurl=sphericalobject(coeff=delSLcurl_fl - ice_time_grid.height_time_coeff[:t_it+1,:].sum(0)- sed_time_grid.height_time_coeff[:t_it+1,:].sum(0)).coefftogrd()
            RO=sphericalobject(grd=self.delSLcurl.grd*self.grd).grdtocoeff()
            self.delPhi_g=np.real(1/self.coeff[0] * (- ice_time_grid.rho/self.rho*ice_time_grid.height_time_coeff[:t_it+1,0].sum() - RO.coeff[0] + TO.coeff[0]))
            #print(t_it,':',self.coeff[0])
            if t_it==1:
                chi = np.abs((np.sum(np.abs(RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[0,:])) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
                self.height_time_coeff[t_it,:]=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[0,:]
            else :
                chi = np.abs((np.sum(np.abs(RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[:t_it,:].sum(0))) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
                self.height_time_coeff[t_it,:]=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[:t_it,:].sum(0)
        return chi

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
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=0,grid_name='time_grid',from_file=(False,)) : 
        """
        Parameters
        ----------
        """
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file)
        # initialize coefficient and grid to 0 because no grid or coefficient has been created
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])