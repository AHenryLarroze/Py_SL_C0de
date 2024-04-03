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
import math

from scipy import io

def sph2cart(az, el, r):
    rsin_theta = r * np.sin(el)
    x = rsin_theta * np.cos(az)
    y = rsin_theta * np.sin(az)
    z = r * np.cos(el)
    return x, y, z

class GRID(object):
        """
        The class _`GRID` is used to represent the Gaussian Grid.

        Attributes
        ----------

            .. note::

                This grid have no attribute * this might be changed in future version to be used as stand alone. The absence of attribute is the result of the use of this function in :ref:`TIME_GRID <TIME_GRID>` that use that define all parameters used in the init function of GRID.

        Methods
        -------
            `interp_on`_ : 
                method used to interpolate the grid over the model grid.
            `smooth_on`_ : 
                method used to smooth the grid to reduce noise effect.
            `disk`_ : 
                method used to create a disk of certain shape to test parameters.
                     

        """
        def __init__(self):
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
            The method _`interp_on` interpolate a grid of data on the grid of the model calculated in the init function of `GRID`_. This function is using stripy library to pursue the interpolation (`Stripy library <https://underworldcode.github.io/stripy/2.0.5b2/FrontPage.html>`_).
        
            Attribute
            --------- 
                grd : np.array([m,n])
                    The grid data over space.
                lon : np.array([m,])
                    The longitude of the data.
                lat : np.array([n,])
                    The latitude of the data.   
                smoothing : bool
                    If a smoothing is applied to the grid before the interpolation, see `smooth_on`_.
                grid_type : str
                    grid type define if it's a grid over the whole world or over a small area. This is used to avoid long computation for the interpolation in the case of interpolating small areas over the world. There is two possible 'global' and 'local'. The global interpolation is using `stripy.sTriangulaion <https://underworldcode.github.io/stripy/2.0.5b2/SphericalMeshing/SphericalTriangulations/Ex3-Interpolation.html>`_. The local interpolation is using `stripyt.Triangulation <https://underworldcode.github.io/stripy/2.0.5b2/SphericalMeshing/CartesianTriangulations/Ex3-Interpolation.html>`_.        
                error : bool
                    If true, the method return the error of the interpolation calculated by stripy.
                
            Return
            ------
                grd : np.array([maxdeg*2,maxdeg])
                    The grid interpolated on the model grid.
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

            '''
            The method _`smooth_on` is smoothing the grid over the area whith `smoothing <https://underworldcode.github.io/stripy/2.0.5b2/SphericalMeshing/CartesianTriangulations/Ex5-Smoothing.html>`_. This function can be used to correct values before the interpolation over time and space. This way, you can have better results on topographic convergence.

            Attribute
            ---------
                grd : np.array([m,n])
                    Array containig the grid values over the space.
                lat : nparray([n,])
                    latitude.
                lon : nparray([m,])
                    longitude.
            
            Return
            ------
                grd : np.array([m,n])
                    smoothed array with the same shape then the initial grid.
            '''

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
            _`disk` is a method used to create a thickness grid. This grid can be used to test different parameters.

            Attribute
            ---------
                lat : nparray([1,])
                    Array contaning the latitudinal coordinate in degree of the center of the disk.
                lon : nparray(1,[])
                    Array containing the longitudinal coordinate in degree of the center of the disk.
                radius : double
                    The radius in degree of the disk in degree.
                high : double
                    The thickness in meter of the disk over the considered area.
            
            Return
            ------
                grd : np.array([maxdeg*2,maxdeg])
                    The grid as defined by the `GRID`_ class with a disk of thikness high at lon,lat position with the size of radius.
        
            '''

            grd=np.zeros((tx,self.lats.size,self.elons.size))
            lon_g,lat_g=np.meshgrid(self.elons,self.lats)
            grd[:,((lon-lon_g)**2+(lat-lat_g)**2)<radius]=high
            grd[:,:2,:]=grd[:,:2,:]*0
            return grd

        def zeros(self,tx=1):
            '''
            _`zeros` is a method used to generate a zero array with the caracteristics of the grid. 
            
            Attribute : 
            -----------
                tx : int
                    times the thickness is repeated

            Return :
            --------
                np.zeros((tx,self.lats.size,self.elons.size))
                    An array containing only zeros 
            '''
            return np.zeros((tx,self.lats.size,self.elons.size))
        
        def along_transect(self,coord=('lat_start','lon_start','lat_stop','lon_stop'),point_density=None,point_distance=None):
            if not(point_density is None):
                theta=np.linspace(coord[0],coord[2],point_density)
                phi=np.linspace(coord[1],coord[3],point_density)
            if not(point_distance is None):
                theta=np.arange(coord[0],coord[2],point_distance)
                phi=np.arange(coord[1],coord[3],point_distance)
            i=0
            Y_lm=np.expand_dims(pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[0]+pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[1]*1j,axis=1)
            
            for i in range(1,len(phi)):
                Y_lm=np.concatenate((Y_lm,np.expand_dims((pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[0]+pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[1]*1j),axis=1)),axis=1)
            return (Y_lm*np.repeat(np.expand_dims(self.coeff,axis=1),point_density,axis=1)).sum(0)

class TIME_GRID(GRID,sphericalobject):
    """
    The _`TIME_GRID` class is used to manage the mass grids. These grids have a time dimenssion this way wa can manage the time variation of the mass. We can define the mass grid by it's mass directly or by coupling a height with a density. If needed you can load spherical harmonics coefficient.Tis class is inheriting the methods from :ref:`sphericalobject <sphericalobject>` and :ref:`GRID <GRID>`. 

    Attributes
    ----------
        time_step : np.array([time_step_number,])
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        maxdeg : int
            Maximum harmonic coefficient degree of the data. this define the chape of the grid and coefficient arrays
        height_time_grid : np.array([maxedg*2,maxdeg])
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array([maxedg*2,maxdeg])
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step
        rho : float
            The density of the considered layer. 
        
        .. note::

            In future development the density may vary threw space and time. We'll have to make a variable object more then a constant density. 

        grid_name : str
            The name of the grid. We recommand you to choose a specific name for each grid you create. This name is used to save the grid in an nc file with `save`_. 
        from_file : (bool,way)
            This parameter define if the data are new or loaded from a previously saved model in a nc file. If the first element is False, the code will create a blank object, based on provided datas. If the first element is True, the method will get the data from the file way specified in the second element of this attribute.
        superinit : bool
            This parameter is used to specify if the object is used as herited method in an initialisation of a child class object.
    
    Methods
    -------
        `interp_on_time`_ 
            Interpolate a grid over the time considered in the model  
        `interp_on_time_and_space`_ :
            Interpolate the grid over time and space as in the defined Grid during the initialisation of the class
        `grid_from_step`_ :
            Get the grid for a defined time iteration
        `coeff_from_step`_ :
            Get the spherical harmonics coefficient for a defined time iteration
        `timegrdtotimecoeff`_ :
            Convert the grid into spherical harmonics coefficient for all time steps
        `timecoefftotimegrd`_ :
            Convert the spherical harmonics coefficient into grid for all time steps
        `zeros_time`_ :
            Generate a zero grid for all time steps
        `disk_time`_ :
            Generate a disk of a specified thickness at a specified location over all time steps
        `update_0`_ :
            Update the 0 time step data of the grid
        `save`_ :
            Save the grid in a specified nc file with the name of the grid

    """

    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=0,grid_name='time_grid',from_file=(False,),superinit=False):
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

            if not(height_time_grid is None):
                self.height_time_grid=height_time_grid
                self.mass_time_grid=height_time_grid*rho
                self.grd_0=self.mass_time_grid[0,:,:]
            elif not(mass_time_grid is None):
                self.mass_time_grid=mass_time_grid
                self.grd_0=self.mass_time_grid[0,:,:]
            elif not(height_time_coeff is None):
                self.height_time_coeff=height_time_coeff
                self.mass_time_coeff=height_time_coeff*rho
                self.coeff_0=self.mass_time_coeff[0,:]
            elif not(mass_time_coeff is None):
                self.mass_time_coeff=mass_time_coeff
                self.coeff_0=self.mass_time_coeff[0,:]

        elif from_file[0] :
            self.ncgrid = Dataset(from_file[1]+'.nc',mode='r',format='NETCDF4_CLASSIC') 
            self.time_grid_name=self.ncgrid.title
            self.maxdeg=len(self.ncgrid['lat'][:].data)
            super().__init__()

            self.time_step=self.ncgrid['time'][:].data
            self.time_step_number=len(self.time_step)
            self.maxdeg=self.ncgrid.dimensions['maxdeg'].size

            self.rho=self.ncgrid['rho'][:].data

            self.height_time_grid=self.ncgrid['thickness'][:].data
            self.mass_time_grid=np.zeros((self.time_step_number,self.nlats,self.nlons))
            self.height_time_coeff=self.ncgrid['coeff_real'][:].data+self.ncgrid['coeff_imag'][:].data*1j
            self.mass_time_coeff=np.zeros((self.time_step_number,int(self.maxdeg*(self.maxdeg+1)/2)))+0j

            self.mass_time_grid=self.height_time_grid*rho

            if not(superinit):
                self.ncgrid.close()

    def interp_on_time(self,grid_to_interp,grid_time_step,model_time_step,interp_type='Thickness_divide',backend='False',grid_type='regular'):
        """
        The function _`interp_on_time` is used for interpolation upon time and space it call the interpolation function of the :ref:`GRID <GRID>` parameter. This function adapt the order of time and space interpolation to reduce computation time. The temporal interpolation try to preserve the thickness of the overall time. Tis is down by cutting and merging time steps of the original grid to match the model time_step. 

        Attributes
        ----------
            grid_to_interp : np.array([k,n,m])
                The grid to be interpreted.
            grid_time_step : np.array([k,])
                The time value of each time step of the grid model.
            model_time_step ; np.array([time_step_number,])
                The time values of the model.
            interp_type : str
                No use of this parameter anymore
            backend : bool 
                Define if the function retur, backends. True it will return the backends, False (default value) don't give any backend. 
        
        Return :
        --------
            grid_interpolated : np.array([time_step_number,n,m])
                The interpolated grid over time. Depending of the model parameters. 
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
                #print(np.concatenate((d_grid_to_interp.repeat(count,axis=0),np.zeros((out,d_grid_to_interp.shape[1],d_grid_to_interp.shape[2])))).shape)
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
                print(type(to_merge))
                if type(to_merge) is 'np.array' :
                    for i in range(len(to_merge)):
                        if backend :
                            print('time slicing : ' + str(i))
                        grd_interpolated[to_merge[i]-1,:]+=grd_interpolated[to_merge[i],:]
                        grd_interpolated=np.delete(grd_interpolated,to_merge[i],0)
                        to_merge-=1
                else :
                    if backend :
                        print('time slicing : ' + str(0))
                    grd_interpolated[to_merge-1,:]+=grd_interpolated[to_merge,:]
                    grd_interpolated=np.delete(grd_interpolated,to_merge,0)
                    to_merge-=1
                return grd_interpolated[::-1,:]
            else :
                print('no such grid type avaiable. try local or global')

    def interp_on_time_and_space(self,grid_to_interp,grid_time_step,grid_lon,grid_lat,interp_type='Thickness_divide',backend=False,grid_type='global'):
        """
        The _`interp_on_time_and_space` function is used for interpolation upon time and space it call the interpolation function of the :ref:`GRID <GRID>` parameter. This function perform the temporal and spatial interpolation in different order to ameliorate the computation time. If the temporal resolution of the input grid is higher than the model time resolution the temporal resolution will be perform first. The spatial resolution is performed first in the other case.

        Attributes
        ----------
            grid_to_interp : np.array([k,n,m])
                The grid to be interpreted.
            grid_time_step : np.array([k,])
                The time value of each time step of the grid_to_interp.
            grid_lon : np.array([n])
                The longitudinal coordinates of the grid_to_interp.
            grid_lat : np.array([m])
                The latitudinal coordinate of the grid_to_interp.
            interp_type : str
                No use of this parameter anymore
            backend : bool 
                Define if the function retur, backends. True it will return the backends, False (default value) don't give any backend. 
        
        Return :
        --------
            grid_interpolated : np.array([time_step_number,maxdeg*2,maxdeg])
                The interpolated grid over time. Depending of the model parameters. 
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
        The _`grid_from_step` method is used to get the value of the grid at the defined time step.
        
        Attributes :
        ------------
            t_it : int
                This is the value of the time step iteration on wich you are trying to retreave the grid. It must inside the time_step interpolation you have used during the initialisation of the time grid. 
        
        Return :
        --------
            None
        """
        # if not(self.mass_time_grid is None) :
        #     self.grd=self.mass_time_grid[t_it,:,:]
        # elif not(self.height_time_grid is None):
        self.grd=self.height_time_grid[t_it,:,:].copy()
        return self

    def coeff_from_step(self,t_it):
        """
        The _`coeff_from_step` method is used to get the value of the coefficient at the requested time iteration.
        
        Attributes :
        ------------
            t_it : double
                This is the value of the time step iteration on wich you are trying to retreave the coefficient. It must be inside the time_step interpolation you have used during the initialisation of the time grid.

        Return :
        --------
            None
        """
        # if not(self.mass_time_coeff is None) :
        #     self.coeff=self.mass_time_coeff[t_it,:]
        # elif not(self.height_time_coeff is None):
        self.coeff=self.height_time_coeff[t_it,:].copy()
        return self
    
    def timegrdtotimecoeff(self):
        """
        The _`timegrdtotimecoeff` method transform for each time step the grid into spherical harmonics coefficient. 

        Attribute :
        -----------
            None
        
        Result :
        --------
            None

        """
        for i in range(self.time_step_number-1):
            self.grd=self.height_time_grid[i,:,:]
            self.height_time_coeff[i,:]=self.grdtocoeff().coeff
        return self
    
    def timecoefftotimegrd(self):
        """
        The _`timecoefftotimegrd` method transform for each time step the spherical harmonic coefficient into a grid. 

        Attribute :
        -----------
            None
        
        Result :
        --------
            None
        """
        for i in range(self.time_step_number-1):
            self.coeff=self.height_time_coeff[i,:]
            self.height_time_grid[i,:,:]=self.coefftogrd().grd
        return self
    
    def zeros_time(self,time_step_number):
        """
        The _`zeros_time` method is used to define a grid over time with only 0 value. It is based on :ref:`GRID.zeros <zeros>`. 

        Attribute :
        -----------
            time_step_number : int
                The number of time step on wich we apply the zeros grid.
        
        Return :
        --------
            None
        """
        self.height_time_grid=self.zeros(time_step_number)
    
    def disk_time(self,time_step_number,lat,lon,radius,high):
        """
        The _`disk_time` method is used to define a grid over time with a disk defined with it's center coordinate and the height. This function is based on :ref:`GRID.disk <disk>`.

        Attribute :
        -----------
            time_step_number : int
                The number of time step on wich the disk load will be applyed.
            lat : double
                The latitude of the center of the disk (°).
            lon : double
                The longitude of the center of the disk (°).
            radius : double
                The radius of the disk (°).
            high : double
                The high of the disk (m).
        Return :
        --------
            None
        """
        self.height_time_grid=self.disk(lat,lon,radius,high,time_step_number)
    
    def update_0(self):
        """
        The _`update_0` function is used to save the first time iteration of the object before it's modification to be called at any moment in the code without alteration.

        Attribute :
        -----------
            None
        Return :
        --------
            None

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

    def point_time(self,coords):
        points=np.zeros((coords.shape[0],self.time_step_number))
        for t_it in range(self.time_step_number):
            coeff=(self.height_time_coeff[t_it-2,:]-self.height_time_coeff[t_it-3,:])/(self.time_step[t_it-2]-self.time_step[t_it-3])
            coeff=np.stack((coeff.real,coeff.imag))
            coeff=pysh.shio.SHCindexToCilm(coeff)
            points[:,t_it]=-pysh.expand.MakeGridPoint(coeff,-coords[:,0],coords[:,1])
        return points
    
    def along_transect(self,coord=('lat_start','lon_start','lat_stop','lon_stop'),point_density=None,point_distance=None,backend=False):
        if not(point_density is None):
            theta=np.linspace(-coord[0],-coord[2],point_density)
            phi=np.linspace(coord[1],coord[3],point_density)
        if not(point_distance is None):
            theta=np.arange(-coord[0],-coord[2],point_distance)
            phi=np.arange(coord[1],coord[3],point_distance)
        coeff=np.stack((self.coeff.real,self.coeff.imag))
        coeff=pysh.shio.SHCindexToCilm(coeff)
        transect=pysh.expand.MakeGridPoint(coeff,theta,phi)
        # Y_lm=np.expand_dims(pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[0]+pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[1]*1j,axis=1)
        # for i in range(1,len(phi)):
        #     Y_lm=np.concatenate((Y_lm,np.expand_dims((pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[0]+pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[1]*1j),axis=1)),axis=1)
        # self.coeff_from_step(0)
        # transect=np.expand_dims((Y_lm*np.repeat(np.expand_dims(self.coeff,axis=1),point_density,axis=1)).sum(0),axis=1)
        # for t_it in range(1,self.height_time_coeff.shape[0]):
        #     if backend :
        #         print(t_it)
        #     self.coeff_from_step(t_it)
        #     transect=np.concatenate((transect,np.expand_dims((Y_lm*np.repeat(np.expand_dims(self.coeff,axis=1),point_density,axis=1)).sum(0),axis=1)),axis=1)
        return transect

    def save(self,save_way='',supersave=False):

        """
        The _`save` function is used to save the grid and all it's parameters inside a nc file. Parameters saved are, longitude, latitude, maximum degree, the time steps of the grid, the thickness of the grid, the harmonic coefficient, grid density. The harmonic coefficient, due to complexe data type management of nc, are saved separately in there complexe and real part. The created file will have the name of the grid. 

        Attribute :
        -----------
            save_way : str
                The filepath where the data will be saved. Default value is the current file
            supersave : bool
                Define if the save is called as a super method from an object that inherit the function. This precise if this method has to close the nc file (False) of if the herited class will do it (True). Default value is False.
        Return :
        --------
            None

        """

        self.ncgrid=Dataset(save_way+'/'+self.time_grid_name+'.nc','w','NETCDF4_CLASSIC')
        self.ncgrid.title=self.time_grid_name

        self.ncgrid.createDimension('maxdeg',self.maxdeg)
        self.ncgrid.createDimension('maxdeg_order',(self.maxdeg)*(self.maxdeg+1)/2)
        self.ncgrid.createDimension('lon',self.nlons)
        self.ncgrid.createDimension('lat',self.nlats)
        self.ncgrid.createDimension('time_diff',self.time_step_number-1)
        self.ncgrid.createDimension('time_step',self.time_step_number)

        lat=self.ncgrid.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:]=self.lats

        lon=self.ncgrid.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        lon[:]=self.elons

        time=self.ncgrid.createVariable('time', np.float32, ('time_step',))
        time.units = 'kyr'
        time.long_name = 'time'
        time[:]=self.time_step

        thickness=self.ncgrid.createVariable('thickness',np.float32,('time_diff','lat','lon'))
        thickness.units='m'
        thickness.long_name='layer_thickness'
        thickness[:,:,:]=self.height_time_grid

        coeff_real=self.ncgrid.createVariable('coeff_real',np.float32,('time_diff','maxdeg_order'))
        coeff_real.units='m'
        coeff_real.long_name='layer thickness in spherical harmonics'
        coeff_real[:,:]=np.real(self.height_time_coeff)

        coeff_imag=self.ncgrid.createVariable('coeff_imag',np.float32,('time_diff','maxdeg_order'))
        coeff_imag.units='m'
        coeff_imag.long_name='layer thickness in spherical harmonics'
        coeff_imag[:,:]=np.imag(self.height_time_coeff)

        rho=self.ncgrid.createVariable('rho',np.float32)
        rho.units='kg/m3'
        rho.long_name='density of the considered grid'
        rho[:]=self.rho
        
        if not(supersave) :
            self.ncgrid.close()

class SEDIMENT_TIME_GRID(TIME_GRID):
    """
    The _`SEDIMENT_TIME_GRID` class is used to represent sediment deposited thickness and time over time. It inherit from :ref:`TIME_GRID <TIME_GRID>` and just add a default value of sediment density at 2600 kg/m3. 

    .. note::
        This class must include the developement of :cite:`ferrier_2017` on sediment compaction and it's effect on water redistriution. 

    Attributes
    ----------
        time_step : np.array([time_step_number,])
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        maxdeg : int
            Maximum harmonic coefficient degree of the data. this define the chape of the grid and coefficient arrays
        height_time_grid : np.array([maxedg*2,maxdeg])
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array([maxedg*2,maxdeg])
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step
        rho : float
            The density of the considered layer. Default value is 2600.
        
        .. note::

            In future development the density may vary threw space and time. We'll have to make a variable object more then a constant density. 

        grid_name : str
            The name of the grid. We recommand you to choose a specific name for each grid you create. This name is used to save the grid in an nc file with `save`_. 
        from_file : (bool,way)
            This parameter define if the data are new or loaded from a previously saved model in a nc file. If the first element is False, the code will create a blank object, based on provided datas. If the first element is True, the method will get the data from the file way specified in the second element of this attribute.

    Method
    ------
        :ref:`save <sed_save>` 
            used to save data

    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=2600,grid_name='time_grid',from_file=(False,)) : 
        """
        Parameters
        ----------
        """
        # preset the grid and coefficient to none. The sediment density is set to 2300 
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True)
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        if from_file[0] :
            self.ncgrid.close()

    def save(self,save_way=''):
        """
        .. _sed_save:

        The save method is used to save the grid data to a file. It call the super method :ref:`save <save>` method with no additional saved parameters.   

        .. note::
            Because we need to include :cite:`ferrier_2017` works, this save function will be modified to include data about the sediment compaction.  

        Attributes
        ----------
            rho : float
                Density value of the sediment as a constant, default value is 2600.

        Return
        ------
            None

        """
        super().save(save_way,supersave=True)
        self.ncgrid.close()

class ICE_TIME_GRID(TIME_GRID):
    """
    The _`ICE_TIME_GRID` class is used to represent the ice thickness evolution threw time. This class inherit of :ref:`TIME_GRID <TIME_GRID>`. 

    ...

    Attributes
    ----------
        time_step : np.array([time_step_number,])
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        maxdeg : int
            Maximum harmonic coefficient degree of the data. this define the chape of the grid and coefficient arrays
        height_time_grid : np.array([maxedg*2,maxdeg])
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array([maxedg*2,maxdeg])
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step
        rho : float
            The density of the considered layer. Default is 916.7 kg/m3.
        
        .. note::

            In future development the density may vary threw space and time. We'll have to make a variable object more then a constant density. 

        grid_name : str
            The name of the grid. We recommand you to choose a specific name for each grid you create. This name is used to save the grid in an nc file with `save`_. 
        from_file : (bool,way)
            This parameter define if the data are new or loaded from a previously saved model in a nc file. If the first element is False, the code will create a blank object, based on provided datas. If the first element is True, the method will get the data from the file way specified in the second element of this attribute.
        
        
    Methods
    -------
        `ice_correction`_ 
            correct the grounded ice thickness from the created floating ice by ground vertical mouvement
        :ref:`save <ice_save>` 
            used to save data
        
    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=916.7,grid_name='time_grid',from_file=(False,)) : 
        """
        Parameters
        ----------
        """
        # initialize to false the grid and coefficient (no grid or coefficient pre loaded). The ice volumetric masse is set to 916.7 kg/m3.
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True)
        # self.isgrd=False
        # self.iscoeff=False
        self.saved=np.array([])
        self.sdeli_00=0
        self.deli_00_prev=0

        if from_file[0] :
            self.ice=self.ncgrid['ice'][:].data
            self.ncgrid.close()

    def ice_correction(self,topo,oc):
        """
        The _`ice_correction` method is used to correct the grounded ice thickness from the floating ice generted by vertical ground motion. This correction is done for each time step to remove the floating ice.  The correction of grounded ice is based on :ref:`Grounded ice correction <Ice_corr>`.

        Attributes
        ----------
            topo : :ref:`TOPOGRAPHIC_TIME_GRID <TOPOGRAPHIC_TIME_GRID>` class object
                A topogrphic grid object, used to check if grounded ice become floating ice. The topography is then modified to by this function. 
            oc : :ref:`OCEAN_TIME_GRID <OCEAN_TIME_GRID>` class object
                An oceanic time grid object used only to get the density of ocean set for the model. 

            .. note::
                To avoid memory consumption the oc parameter should be replace simply by oc_rho the ocean density. 
                 

        Return
        ------
            None

        """
        for t_it in range(self.time_step_number-1):
            if t_it==0 :
                topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]-self.height_time_grid[t_it,:,:]+self.ice[t_it,:,:]
            else :
                topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]-self.height_time_grid[:t_it+1,:,:].sum(0)+self.ice[:t_it+1,:,:].sum(0)
        for t_it in range(self.time_step_number-1): 
            if t_it==0:
                check1 = OCEAN_TIME_GRID().evaluate_ocean(-topo.height_time_grid[t_it,:,:]+self.ice[t_it,:,:]) # generate the ocean function for ice-topo
                check2 = OCEAN_TIME_GRID().evaluate_ocean(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:]).grd*(OCEAN_TIME_GRID().evaluate_ocean(-self.ice[t_it,:,:]*self.rho-(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:])*oc.rho).grd)
                self.height_time_grid[t_it,:,:] = check1.grd*self.ice[t_it,:,:]+check2*self.ice[t_it,:,:] # add the two part of ice over the check1 nd check2 positive area.
            else :
                check1 = OCEAN_TIME_GRID().evaluate_ocean(-topo.height_time_grid[t_it,:,:]+self.ice[:t_it+1,:,:].sum(0)) # generate the ocean function for ice-topo
                check2 = OCEAN_TIME_GRID().evaluate_ocean(topo.height_time_grid[t_it,:,:]-self.ice[:t_it+1,:,:].sum(0)).grd*(OCEAN_TIME_GRID().evaluate_ocean(-self.ice[:t_it+1,:,:].sum(0)*self.rho-(topo.height_time_grid[t_it,:,:]-self.ice[:t_it+1,:,:].sum(0))*oc.rho).grd)
                self.height_time_grid[t_it,:,:] = check1.grd*self.ice[:t_it+1,:,:].sum(0)+check2*self.ice[:t_it+1,:,:].sum(0)-self.height_time_grid[:t_it,:,:].sum(0) # add the two part of ice over the check1 nd check2 positive area.
        for t_it in range(self.time_step_number-1):    
            if t_it==0 :
                topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]+self.height_time_grid[t_it,:,:]-self.ice[t_it,:,:]
            else :
                topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]+self.height_time_grid[:t_it+1,:,:].sum(0)-self.ice[:t_it+1,:,:].sum(0)

    def save(self,save_way=''):

        """
        .. _ice_save:

        The save method is used to save the data of the ice grid. Because of the `ice_correction`_ method that update the grid we choosed to preserve the original ice thickness data in a ice parameter that is saved inside the nc file. Otherwise this method use the super method :ref:`save` to save the rest of the data.  

        Attributes
        ----------
            save_way :str
                The way where the nc file is saved. Default value is the current file (an empty str).              

        Return
        ------
            None

        """

        super().save(save_way,supersave=True)
        ice=self.ncgrid.createVariable('ice',np.float32,('time_diff','lat','lon'))
        ice.units='m'
        ice.long_name='initial_ice_thickness'

        if self.__dict__.keys().__contains__('ice'):
            ice[:,:,:]=self.ice
        else :
            ice[:,:,:]=self.height_time_grid.copy()

        self.ncgrid.close()

class OCEAN_TIME_GRID(TIME_GRID):
    """
    The _`OCEAN_TIME_GRID` class is used to represent the ocean thickness variation and contains the method to resolve the sea level equation. This method inherit from :ref:`TIME_GRID <TIME_GRID>`.

    ...

    Attributes
    ----------
        time_step : np.array([time_step_number,])
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        maxdeg : int
            Maximum harmonic coefficient degree of the data. this define the chape of the grid and coefficient arrays
        height_time_grid : np.array([maxedg*2,maxdeg])
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array([maxedg*2,maxdeg])
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step
        rho : float
            The density of the considered layer. Default is 1000 kg/m3.
        
        .. note::

            In future development the density may vary threw space and time. We'll have to make a variable object more then a constant density. 

        grid_name : str
            The name of the grid. We recommand you to choose a specific name for each grid you create. This name is used to save the grid in an nc file with `save`_. 
        from_file : (bool,way)
            This parameter define if the data are new or loaded from a previously saved model in a nc file. If the first element is False, the code will create a blank object, based on provided datas. If the first element is True, the method will get the data from the file way specified in the second element of this attribute.
           
    Methods
    -------
        :ref:`update_0 <update_0_oce>`
            set the grd_0 from the actual grd loaded in the grid
        `evaluate_ocean`_
            evaluate the ocean function based on the Gaussian grid of the topography 

    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=1000,grid_name='time_grid',from_file=(False,)) :
        """
        Parameters
        ----------
        """
        # initialize the ocean with no grid and no coefficient. The saved is also initialized. The volumetric mass of water is set to 1000 

        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True)
        self.saved=np.array([])
        if from_file[0] :
            self.ncgrid.close()
        
    def update_0(self):
        """
        .. _update_0_oce:

        The update_0 method update the grd_0 parameter of the object to the currend loaded grd. 

        Attributes
        ----------
            None
            
        Return
        -------
            None

        """
        self.grd_0=self.grd.copy()
        
    def evaluate_ocean(self,topo) :
        '''
        The _`evaluate_ocean` method evaluate the ocean function using the topography. It create a 0-1 matrix wich is 1 where topo<0 and 0 where topo>0.  The ocean function is described in :ref:`ocean function <oc_func>`.
    
        Attribute
        ---------
            topo : np.array(maxdeg,maxdegx2)
                topographic gaussian grid.
        
        Returns :
            None

        '''
        # use sign function to optimize the conversion from positive negative value to boolean value.
        out = -0.5*np.sign(topo)+0.5
        out = 0.5*np.sign(out-0.6)+0.5
        # set the grd to the output of the previous computation.
        self.grd=out
        return self
    
    def sea_level_solver(self,load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,conv_lim):
        '''
        The _`sea_level_solver` method solve the sea level equation until. Beacause of the iterative type of the resolution of the SLE, this method define also a first guess of the Sea level at the first iteration and the first time step. This function is based on the convergence iteration for the estimation of the variability defined in :ref:`Convergence parameter <conv>`.

        Attribute
        ---------
            load : :ref:`LOAD_TIME_GRID <LOAD_TIME_GRID>` class object
                The load time grid as specified in the class object. This grid needs to be of the same shape (maxdeg) then the one of the current object. 
            ice_time_grid : :ref:`ICE_TIME_GRID <ICE_TIME_GRID>` class object
                The ice time grid as specified in the class object. This grid needs to be of the same shape (maxdeg) then the one of the current object. 
            sed_time_grid : :ref:`SEDIMENT_TIME_GRID <SEDIMENT_TIME_GRID>` class object
                The sediment time grid as specified in the class object. This grid needs to be of the same shape (maxdeg) then the one of the current object. 
            love_number : :ref:`LOVE <LOVE>` class object 
                The love numbers as specified in the class object. The love numbers must have been set up with the same maximm degree thne the currend object.
            TO : :ref:`sphericalobject <sphericalobject>` class object
                The ocean contours variability area computed as a sphericalobject class computed for the previous iteration. !Trouver où définir ce calcul!. 
            t_it : int
                The time iteration of the current computation on wich apply the resolution of the SLE.
            conv_it : int
                convergence iteration set to 0 if it's for a simple resolution of the SLE on one time step. This is used when you are working on a topographic convergence. In the code, the first guess for the SLE will be if it's not the first topographic convergence iteration, the guess of the previuous one.
            conv_lim : float
                To stop the convergence of the solution, the conv_lim is usually set to 10^-3. The number of required step is then between 13 and 7. 
        Return 
        ------
            None 

        '''
        if conv_it==0 :
            if t_it==0:
                self.height_time_coeff[t_it,:]=self.prev/self.prev[0]*(TO.coeff[0]-TO.prev[0])-TO.coeff-TO.prev
            else :
                self.height_time_coeff[t_it,:]=self.prev/self.prev[0]*(-ice_time_grid.rho/self.rho*ice_time_grid.height_time_coeff[:t_it+1,:].sum(0) + TO.coeff[0]-TO.prev[0])-TO.coeff-TO.prev
        chi=np.inf
        while chi>conv_lim:
            chi=self.sea_level_equation(load,ice_time_grid,sed_time_grid,love_number,TO,t_it)
            conv_it+=1
        #if t_it >1 :
            #print(load.rot_pot[t_it-1,:],sed_time_grid.height_time_coeff[t_it,:].max())
        return conv_it
    
    def sea_level_equation(self,load,ice_time_grid,sed_time_grid,love_number,TO,t_it):
        '''
        The _`sea_level_equation` method calculate the Sea level variation following the SLE. Tis function is resolving both the conservation of mass equation and the SL variation. This follows the method described in :ref:`Resolution of SLE including the deconvolution<spec_sol>`.

        Attribute
        ---------
            load : :ref:`LOAD_TIME_GRID <LOAD_TIME_GRID>` class object
                The load time grid as specified in the class object. This grid needs to be of the same shape (maxdeg) then the one of the current object. 
            ice_time_grid : :ref:`ICE_TIME_GRID <ICE_TIME_GRID>` class object
                The ice time grid as specified in the class object. This grid needs to be of the same shape (maxdeg) then the one of the current object. 
            sed_time_grid : :ref:`SEDIMENT_TIME_GRID <SEDIMENT_TIME_GRID>` class object
                The sediment time grid as specified in the class object. This grid needs to be of the same shape (maxdeg) then the one of the current object. 
            love_number : :ref:`LOVE <LOVE>` class object 
                The love numbers as specified in the class object. The love numbers must have been set up with the same maximm degree thne the currend object.
            TO : :ref:`sphericalobject <sphericalobject>` class object
                The ocean contours variability area computed as a sphericalobject class computed for the previous iteration. !Trouver où définir ce calcul!. 
            t_it : int
                The time iteration of the current computation on wich apply the resolution of the SLE.
        Return 
        ------
            None 
            
        '''
        if t_it==0 :
            delSLcurl_fl=love_number.E* love_number.T.coeff *(ice_time_grid.height_time_coeff[0,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[0,:]*sed_time_grid.rho+self.height_time_coeff[0,:]*self.rho)
            self.delSLcurl=sphericalobject(coeff=delSLcurl_fl - ice_time_grid.height_time_coeff[0,:]- sed_time_grid.height_time_coeff[0,:]).coefftogrd()
            load.calc_rot_visc(ice_time_grid.height_time_coeff[1:t_it,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[1:t_it,:]*sed_time_grid.rho+self.height_time_coeff[1:t_it,:]*self.rho,0,love_number)
            RO=sphericalobject(grd=self.delSLcurl.grd*self.grd).grdtocoeff()
            self.delPhi_g=np.real(1/self.coeff[0] * (- ice_time_grid.rho/self.rho*ice_time_grid.coeff[0] - RO.coeff[0] + TO.coeff[0]))
            chi = np.abs( (np.sum(np.abs(RO.coeff + self.delPhi_g*self.coeff -  TO.coeff)) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
            self.height_time_coeff[t_it,:]=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff
        else :
            
            if t_it == 1 : 
                #load.calc_viscuous_load(ice_time_grid.height_time_coeff[0,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[0,:]*sed_time_grid.rho+self.height_time_coeff[0,:]*self.rho,love_number.beta_l,0)
                load.calc_rot_visc(ice_time_grid.height_time_coeff[1,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[1,:]*sed_time_grid.rho+self.height_time_coeff[1,:]*self.rho,t_it,love_number)
                load.rot_pot[t_it,:]=load.rot_pot[t_it,:]
                load.calc_rotational_viscuous(np.zeros(6),love_number.beta_l_tide,0)
                load.calc_viscuous(np.zeros(ice_time_grid.height_time_coeff[0,:].shape),love_number.beta_l,0)

            else : 
                load.calc_rot_visc(ice_time_grid.height_time_coeff[1:t_it,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[1:t_it,:]*sed_time_grid.rho+self.height_time_coeff[1:t_it,:]*self.rho,t_it,love_number)
                load.rot_pot[t_it,:]=load.rot_pot[t_it,:]-load.rot_pot[:t_it,:].sum(0) # calc small variation
                load.calc_viscuous(ice_time_grid.height_time_coeff[1:t_it,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[1:t_it,:]*sed_time_grid.rho+self.height_time_coeff[1:t_it,:]*self.rho,love_number.beta_l,t_it)
                load.calc_rotational_viscuous(load.rot_pot[1:t_it,:],love_number.beta_l_tide,t_it)
            # print(t_it)
            # print('rot_pot : ',load.rot_pot[:t_it])
            # print('m : ',load.sdelm[:t_it])
            # print('I : ', load.sdelI[:t_it])
            # print('V_lm_T : ', load.V_lm_tide.coeff)
            delSLcurl_tide_fl=1/load.g*load.V_lm_tide.coeff+1/load.g*love_number.E_T[:6]*load.rot_pot[:t_it,:].sum(0)
            
            delSLcurl_fl=love_number.E* love_number.T.coeff *(ice_time_grid.height_time_coeff[1:t_it+1,:].sum(0)*ice_time_grid.rho+sed_time_grid.height_time_coeff[1:t_it+1,:].sum(0)*sed_time_grid.rho+self.height_time_coeff[1:t_it+1,:].sum(0)*self.rho)+love_number.T.coeff*load.V_lm.coeff

            #print(load.V_lm_tide.coeff.shape,load.V_lm.coeff.shape)

            delSLcurl_tide_fl=np.concatenate((delSLcurl_tide_fl,np.zeros(delSLcurl_fl.shape[0]-delSLcurl_tide_fl.shape[0])))

            self.delSLcurl=sphericalobject(coeff=delSLcurl_fl + delSLcurl_tide_fl - ice_time_grid.height_time_coeff[1:t_it+1,:].sum(0)- sed_time_grid.height_time_coeff[1:t_it+1,:].sum(0)).coefftogrd()
            RO=sphericalobject(grd=self.delSLcurl.grd*self.grd).grdtocoeff()
            self.delPhi_g=np.real(1/self.coeff[0] * (- ice_time_grid.rho/self.rho*ice_time_grid.height_time_coeff[1:t_it+1,0].sum() - RO.coeff[0] + TO.coeff[0]))
            #print(t_it,':',self.coeff[0])
            if t_it==1:
                chi = np.abs((np.sum(np.abs(RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[0,:])) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
                self.height_time_coeff[t_it,:]=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[0,:]
            else :
                chi = np.abs((np.sum(np.abs(RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[:t_it,:].sum(0))) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
                self.height_time_coeff[t_it,:]=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.height_time_coeff[:t_it,:].sum(0)
            # print('chi : ',chi)
        return chi

class TOPOGRAPHIC_TIME_GRID(TIME_GRID):
    """
    The _`TOPOGRAPHIC_TIME_GRID` class is used to save and include all the topographic variations. This class inherits of :ref:`TIME_GRID <TIME_GRID>`. This class main difference with TIME_GRID is the presence of a parameter called topo_pres wich is the present topography. It is created using :ref:`Precomputation <Precomputation>`.

    Attributes
    ----------
        time_step : np.array([time_step_number,])
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        maxdeg : int
            Maximum harmonic coefficient degree of the data. this define the chape of the grid and coefficient arrays
        height_time_grid : np.array([maxedg*2,maxdeg])
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array([maxedg*2,maxdeg])
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step
        rho : float
            The density of the considered layer.
        
        .. note::

            In future development the density may vary threw space and time. We'll have to make a variable object more then a constant density. 

        grid_name : str
            The name of the grid. We recommand you to choose a specific name for each grid you create. This name is used to save the grid in an nc file with `save`_. 
        from_file : (bool,way)
            This parameter define if the data are new or loaded from a previously saved model in a nc file. If the first element is False, the code will create a blank object, based on provided datas. If the first element is True, the method will get the data from the file way specified in the second element of this attribute.

    Methods
    -------
        :ref:`save <topo_save>` 
            Method to save the topographic datas   

    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=0,grid_name='time_grid',from_file=(False,)) : 
        """
        Parameters
        ----------
        """
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True)
        if from_file[0] :
            self.topo_pres=self.ncgrid['topo_pres'][:].data
            self.ncgrid.close()
        
        # initialize coefficient and grid to 0 because no grid or coefficient has been created
        # self.isgrd=False
        # self.iscoeff=False
        #self.saved=np.array([])

    def save(self,save_way=''):
        """
        .. _topo_save:

        The save method is used to save the data of the topographic grid. Particularity of the topography is the present day topography used in the code to converge toward it. The function save is. Otherwise this method use the super method :ref:`save` to save the rest of the data.  

        Attributes
        ----------
            save_way :str
                The way where the nc file is saved. Default value is the current file (an empty str).              

        Return
        ------
            None

        """
        super().save(save_way,supersave=True)
        topo_pres=self.ncgrid.createVariable('topo_pres',np.float32,('lat','lon'))
        topo_pres.units='m'
        topo_pres.long_name='present topography'

        if self.__dict__.keys().__contains__('topo_pres'):
            topo_pres[:,:]=self.topo_pres

        self.ncgrid.close()


from .love import get_tlm

class LOAD_TIME_GRID(TIME_GRID) :
    """
    The _`LOAD_TIME_GRID` class is used to save and include all the topographic variations. This class inherits of :ref:`TIME_GRID <TIME_GRID>` and :ref:`LOAD <LOAD>`.

    Attributes
    ---------- 
        sdelL : np.array([time_step_number,maxdeg,maxdegx2])
            The load variation grid used to compute earth vertical motion.
        betal : np.array([time_step_number,time_step_number,maxdeg])
            The beta love number as described in :ref:`Variation of geoïd and ground Equations <geoid_ground_variation_theory>` section. There calculated in the :ref:`LOVE <LOVE>` class.
        E : np.array([(maxdeg+1)(maxdeg+2)/2,])
            The elastic component of the earth as love numbers computed form :ref:`LOVE <LOVE>` class.
        a : float
            The earth radius in meter, set by default to 7371000 meters.
        Me : float
            The earth mass in set by default to 5000. 
        time_step : np.array([time_step_number,])
            This array contains the time step of the data you are importing. They will be use for temporal interpolation.
        maxdeg : int
            Maximum harmonic coefficient degree of the data. this define the chape of the grid and coefficient arrays
        height_time_grid : np.array([maxedg*2,maxdeg])
            This array is the height grid at each time steps defined in grid_time_step
        mass_time_grid : np.array([maxedg*2,maxdeg])
            This array is the mass grid at each time steps defined in grid_time_step
        height_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the height spherical harmonic coefficient at each time steps defined in grid_time_step
        mass_time_coeff : np.array([(maxdeg+1)(maxedg+2)/2,])
            This array is the mass spherical harmonic coefficient at each time steps defined in grid_time_step
        rho : float
            The density of the considered layer.
        grid_name : str
            The name of the grid. We recommand you to choose a specific name for each grid you create. This name is used to save the grid in an nc file with `save`_. 
        from_file : (bool,way)
            This parameter define if the data are new or loaded from a previously saved model in a nc file. If the first element is False, the code will create a blank object, based on provided datas. If the first element is True, the method will get the data from the file way specified in the second element of this attribute.

    Methods
    -------
        `calc_viscuous`_:
            Compute the viscous motion of the geoïd and ground for one time step.
        `calc_viscuous_time`_ :
            Compute the viscous ground motion of the earth on all time steps.    
        `calc_elastic_time`_ :
            Compute the elastic ground motion of the earth on all time steps.
        :ref:`save <load_save>` :
            Save the load data

    """
    def __init__(self,sdelL=np.array([]),beta_l=np.array([]),E=np.array([]),E_T=np.array([]),a=7371000,Me=5000,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=0,grid_name='time_grid',from_file=(False,),g=9.80665):
        TIME_GRID.__init__(self,time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True)
        #LOAD.__init__(self,maxdeg,time_step)
        self.beta_counter=np.repeat(np.arange(0,self.maxdeg),np.arange(1,self.maxdeg+1))
        if from_file[0]:
            self.a=self.ncgrid['a'][:].data
            self.Me=self.ncgrid['Me'][:].data
            self.viscuous_deformation=self.ncgrid['viscuous_deformation_real'][:].data+self.ncgrid['viscuous_deformation_imag'][:].data*1j
            self.elastic_deformation=self.ncgrid['elastic_deformation_real'][:].data+self.ncgrid['elastic_deformation_imag'][:].data*1j
            #self.beta_l=self.ncgrid['beta_l'][:].data
            self.elastic_love=self.ncgrid['elastic_love'][:].data
            self.load=self.ncgrid['load_real'][:].data+self.ncgrid['load_imag'][:].data*1j
        else :
            self.g=g
            self.viscuous_deformation=np.zeros((self.time_step_number-2,int(maxdeg*(maxdeg+1)/2)))+0j
            self.elastic_deformation=np.zeros((self.time_step_number-2,int(maxdeg*(maxdeg+1)/2)))+0j
            self.load=sdelL[1:,:]
            self.elastic_love=E[self.beta_counter.astype(int)]
            self.elastic_tidal_love=E_T
            self.a=a
            self.Me=Me
            self.beta_l=beta_l
            self.viscuous_love=self.beta_l[:,:,self.beta_counter.astype(int)]
            self.V_lm=sphericalobject(coeff=np.zeros((int(maxdeg*(maxdeg+1)/2)))+0j)
            self.V_lm_tide=sphericalobject(coeff=np.zeros((6,))+0j)
            self.rot_pot=np.zeros((self.time_step_number-1,6))+0j
            self.sdelI=np.zeros((self.time_step_number-2,3))+0j
            self.sdelm=np.zeros((self.time_step_number-2,3))+0j
        self.T = sphericalobject(coeff=get_tlm(self.maxdeg-1,self.a,self.Me))

    def calc_viscuous(self,sdelL,beta,t_it):
        '''
        The _`calc_viscuous` method is used to calculate the ground and geoïd deformation based on viscuous love numbers.

        Attribute
        ---------
            sdelL : np.array([t_it,(maxdeg+1)(maxdeg+20/2)])
                The load grid used to estimate the ground vertical mouvement. This include all previous loading history because of the viscous comportment of earth. 
            beta : np.array([time_step_number,time_step_number,(maxdeg+1)(maxdeg+2)/2])
                The beta love numbers used to compute the earth deformation to include the viscous part. These love numbers are particularly heavy in the memory due to the representation of the time. 
            t_it : int
                The time iteration at wich the computation is performed.
            
        '''
        if t_it==1 :
            self.V_lm.coeff=np.zeros(beta.shape[2])
        else :
            self.V_lm.coeff=(beta[t_it-1,:t_it-1]*sdelL).sum(0)

    def calc_rotational_viscuous(self,rot_pot,beta_tide,t_it):
        if t_it==1 :
            self.V_lm_tide.coeff=np.zeros(6)
        else :
            self.V_lm_tide.coeff=(beta_tide[t_it-1,:t_it-1]*rot_pot).sum(0)

    def calc_viscuous_time(self,backend=False) :
        '''
        The _`calc_viscuous_time` method compute the vicuous vertical ground motion. This method call the :ref:`LOAD <LOAD>` method for this. 

        Attribute
        ---------
            backend : bool
                Specifie if the method give backend (True) or not (False). Default is False.

        Return
        ------
            None

        '''
        for t_it in range(1,self.time_step_number-1):
            if t_it+1==1 :
                self.calc_viscuous(np.zeros(self.height_time_coeff[0,:].shape),self.viscuous_love,0)
            else :
                self.calc_viscuous(self.load[1:t_it,:],self.viscuous_love,t_it)
            self.viscuous_deformation[t_it-1,:]=self.T.coeff*np.squeeze(self.V_lm.coeff.T)
            if backend:
                print(f'viscuous calculation at {self.time_step[t_it]} kyr done')
    
    def calc_elastic_time(self):
        '''
        The _`calc_elastic_time` method compute the elactic vertical ground motion. This method call the :ref:`LOAD <LOAD>` method for this. 

        Attribute
        ---------
            None

        Return
        ------
            None
             
        '''
        self.elastic_deformation=np.repeat(np.expand_dims(self.T.coeff,axis=0),self.time_step_number-2,axis=0)*np.repeat(np.expand_dims(self.elastic_love,axis=0),self.time_step_number-2,axis=0)*self.load

    def save(self,save_way=''):
        '''
        .. _load_save:

        The save method is used to save the data from the class. It is based on the inherited :ref:`save <save>` method. Because of the particularity of this TIME_GRID, we had to save the new parameters, and calculated data. This function save for each data the real and complex part of the data due to the nc file particularity. The saved data are, The load (load), the viscuous groud motion (viscous_deformation), the elastic ground motion (elastic_deformation), the elastic love numbers (elastic_love), earth radius (a), earth mass (Me). 

        Attribute
        ---------
            save_way : str
                file path to where the grid will be saved. 

        Return
        ------
            None
             
        '''
        super().save(save_way,supersave=True)

        self.ncgrid.createDimension('time_no_init',self.time_step_number-2)

        load_real=self.ncgrid.createVariable('load_real',np.float32,('time_no_init','maxdeg_order'))
        load_real.units='kg'
        load_real.long_name='load grid used to compute the earth deformation.'
        load_real[:,:]=np.real(self.load)

        load_imag=self.ncgrid.createVariable('load_imag',np.float32,('time_no_init','maxdeg_order'))
        load_imag.units='kg'
        load_imag.long_name='load grid used to compute the earth deformation.'
        load_imag[:,:]=np.imag(self.load)

        viscuous_deformation_real=self.ncgrid.createVariable('viscuous_deformation_real',np.float32,('time_no_init','maxdeg_order'))
        viscuous_deformation_real.units='mm/yr'
        viscuous_deformation_real.long_name='viscuous component of the earth deformation due to load.'
        viscuous_deformation_real[:,:]=np.real(self.viscuous_deformation)

        viscuous_deformation_imag=self.ncgrid.createVariable('viscuous_deformation_imag',np.float32,('time_no_init','maxdeg_order'))
        viscuous_deformation_imag.units='mm/yr'
        viscuous_deformation_imag.long_name='viscuous component of the earth deformation due to load.'
        viscuous_deformation_imag[:,:]=np.imag(self.viscuous_deformation)

        elastic_deformation_real=self.ncgrid.createVariable('elastic_deformation_real',np.float32,('time_no_init','maxdeg_order'))
        elastic_deformation_real.units='mm/yr'
        elastic_deformation_real.long_name='elastic component of the earth deformation due to load.'
        elastic_deformation_real[:,:]=np.real(self.elastic_deformation)

        elastic_deformation_imag=self.ncgrid.createVariable('elastic_deformation_imag',np.float32,('time_no_init','maxdeg_order'))
        elastic_deformation_imag.units='mm/yr'
        elastic_deformation_imag.long_name='elastic component of the earth deformation due to load.'
        elastic_deformation_imag[:,:]=np.imag(self.elastic_deformation)

        elastic_love=self.ncgrid.createVariable('elastic_love',np.float32,('maxdeg_order'))
        elastic_love.units='none'
        elastic_love.long_name='elastic load love numbers used to compute the earth deformation'
        elastic_love[:]=self.elastic_love

        # beta_l=self.ncgrid.createVariable('beta_l',np.float32,('time_diff','time_diff','maxdeg'))
        # beta_l.units='none'
        # beta_l.long_name='viscuous load love numbers used to compute the earth deformation'
        # beta_l[:,:]=self.beta_l

        a=self.ncgrid.createVariable('a',np.float32)
        a.units='m'
        a.long_name='earth radius'
        a=self.a

        Me=self.ncgrid.createVariable('Me',np.float32)
        Me.units='kg'
        Me.long_name='earth mass'
        Me=self.Me

        self.ncgrid.close()

    def clean_memory(self):
        '''
        This method is used to clean the memory to avoïd over charging RAM. 
        '''
        self.viscuous_love=0

    def calc_rot_visc(self,sdelL,t_it,love_number,G=6.67408E-11,a=6371000,C=8.034e37,k_f=0.942,omega=7.292E-5):
        '''
        .. note::
            This function must be updated to work with the new versions of the code.

        '''
        # extract degree 2 coefficient from the load
        CminA = (k_f*(a**5)*(omega)**2)/(3*G)

        if t_it==1 :
            L00 = sdelL[0]
            L20 = sdelL[3]
            L21 = sdelL[4]
        else :
            L00 = sdelL[:,0].sum(0)
            L20 = sdelL[:,3].sum(0)
            L21 = sdelL[:,4].sum(0)

        # calculate the load effect constant 
        I1=math.sqrt(32/15)*math.pi*a**4*np.real(L21)
        I2=math.sqrt(32/15)*math.pi*a**4*np.imag(L21)   
        I3=8/3*math.pi*a**4*(L00-L20/math.sqrt(5))
        I=np.array([I1,I2,I3])
        if t_it==1 :
            V_lm=np.zeros(3)
            V_lm_T=np.zeros(3)
        if t_it==2 :
            V_lm = love_number.beta_konly_l[t_it-2,:t_it-1]*self.sdelI[:t_it-1,:].squeeze()
            V_lm_T = love_number.beta_konly_l_tide[t_it-2,:t_it-1]*self.sdelm[:t_it-1,:].squeeze()
        else : #apply the visco elastic properties of the earth on the rotation using the load.
            V_lm = np.dot(love_number.beta_konly_l[t_it-2,:t_it-1],self.sdelI[:t_it-1,:].squeeze())
            V_lm_T = np.dot(love_number.beta_konly_l_tide[t_it-2,:t_it-1],self.sdelm[:t_it-1,:].squeeze())
        
        temp = 1/(1-love_number.k_tide_e[1]/k_f)*(1/CminA * ((1+love_number.k_e[1])*I + V_lm.squeeze()) + V_lm_T.squeeze()/k_f)
        # calculate the perturbation to the rotational potential from Milne 1998
        m1=temp[0]
        m2=temp[1]
        temp = -1/(1-love_number.k_tide_e[1]/k_f)*(1/C * ((1+love_number.k_e[1])*I + V_lm.squeeze()))
        m3=temp[2]

        m=np.array([m1,m2,m3])
        # print(t_it)
        # print('L : ', sdelL)
        # print('L : ', L00,L20,L21)
        # print('temp : ',(1/CminA * ((1+love_number.k_e[1])*I + V_lm.squeeze()) + V_lm_T.squeeze()/k_f))
        # print('CminA : ',CminA)
        # print(t_it)
        # print('V_lm : ',V_lm)
        # print('V_lm_T : ', V_lm_T)
        # print('I : ',I)
        # print('m : ',m)
        #update the rotational load using.
        #print(self.sdelI[:t_it-1,:].sum(0),self.sdelI[:t_it-1,:].sum(0).shape,I.T.squeeze(),I.T.squeeze().shape,self.sdelI[t_it-1,:],self.sdelI[t_it-1,:].shape)
        self.sdelI[t_it-1,:] = I.T.squeeze() - self.sdelI[:t_it-2,:].sum(0)
        self.sdelm[t_it-1,:] = m.T.squeeze() - self.sdelm[:t_it-2,:].sum(0)

        #print('m : ',self.sdelm[t_it-1,:])
        # calculate the rotational perturbation of the earth potential, just for the 6 first coefficient (no use to calculate further)
        self.rot_pot[t_it-1,0] = a**2 * omega**2/3 * (np.sum(m**2) + 2*m3)+0j
        self.rot_pot[t_it-1,3] = a**2 * omega**2/(6*np.sqrt(5)) * (m1**2 + m2**2 - 2*m3**2 - 4*m3)+1j*0
        self.rot_pot[t_it-1,4] = a**2 * omega**2/np.sqrt(30) * (m1*(1+m3) - 1j*m2*(1+m3))+1j*0
        self.rot_pot[t_it-1,5] = a**2 * omega**2/(np.sqrt(5) * np.sqrt(24)) * ( (m2**2-m1**2) + 1j*2*m1*m2 )+1j*0
        
        #print('rot_pot : ',self.rot_pot[t_it-1])