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
import numpy as np
from scipy.interpolate import griddata,RegularGridInterpolator
from .love import calc_beta_counter

from scipy import io

def sph2cart(az, el, r):
    rsin_theta = r * np.sin(el)
    x = rsin_theta * np.cos(az)
    y = rsin_theta * np.sin(az)
    z = r * np.cos(el)
    return x, y, z

class GRID(object):
        """
        The GRID object is the initial element of our model. This object manage grided object and it's associated elements
        """
        def __init__(self,maxdeg=None,P_lm=None):
            """
            The class _`GRID` is used to represent the Gaussian Grid.

            Attributes
            ----------

                .. note::

                    This grid have no attribute * this might be changed in future version to be used as stand alone. The absence of attribute is the result of the use of this function in :ref:`TIME_GRID <TIME_GRID>` that define all parameters used in the init function of GRID.

            Methods
            -------
                `interp_on`_ : 
                    method used to interpolate a grid over the model grid.
                `smooth_on`_ : 
                    method used to smooth the grid to reduce noise effect.
                `disk`_ : 
                    method used to create a disk of certain shape.
                `zeros`_ : 
                    method used to create a 0 grid.
            """
            if maxdeg==None : 
                self.maxdeg=self.maxdeg
            else :
                self.maxdeg=maxdeg
            self.nlons = (self.maxdeg) * 2 
            self.nlats = (self.maxdeg)
            x,w=pysh.expand.SHGLQ(self.maxdeg-1)
            # x=x[::-1]
            x_GL = np.arccos(x)*180/math.pi - 90
            lon_GL = np.linspace(0,360,2*self.maxdeg+1)
            lon_GL = lon_GL[:-1]
            self.lats=x_GL.copy()
            self.elons=lon_GL.copy()
            self.colats = 90 - self.lats
            if P_lm is None :
                self.P_lm=np.zeros((self.maxdeg+1,self.maxdeg+1,len(x)))
                for i,x_i in enumerate(x) :
                    self.P_lm[:,:,i]=pysh.legendre.legendre(self.maxdeg,x_i,csphase=1,cnorm=1)#*correction
            else : 
                self.P_lm=P_lm.copy()
            self.N=(self.maxdeg+1)*(self.maxdeg+2)/2
            
        def interp_on(self,grd,lon,lat,smoothing=False,grid_type='global',error=False,interp_method='stripy'):
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
            if interp_method=='stripy':
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
                    grd,err=spherical_triangulation.interpolate_linear(vertices_elons,vertices_lats,data=grd.flatten())
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
                    # print(elons_in,lats_in,lon,lat)
                    grd_in,err=spherical_triangulation.interpolate_linear(elons_in,lats_in,grd)
                    # I need to complete the grid where the is no data
                    grd=np.zeros(elons.shape)
                    grd[point_in]=grd_in
                else:
                    print('No such grid type, please select one between : global or local')
                if error : 
                    return grd.reshape((self.nlats,self.nlons)), err.reshape((self.nlats,self.nlons))
                return grd.reshape((self.nlats,self.nlons))
            elif interp_method=='Old' :
                # # Extract current ice grid
                # grided_lon,grided_lat=np.meshgrid(self.elons,self.lats)
                # nointerp = grd.copy()  # Transpose to match MATLAB behavior
                # lat_extended = np.concatenate(([90], lat, [-90]))
                # long_extended = np.concatenate(([-0.5], lon, [360.5])) 
                # # Extend rows at the top, bottom, and sides
                # top_row = np.zeros((1, len(lon)))
                # bottom_row = nointerp[0, -1] * np.ones((1, len(lon)))
                # extended = np.vstack((top_row, nointerp, bottom_row))

                # left_column = extended[:, -1].reshape(-1, 1)
                # right_column = extended[:, 0].reshape(-1, 1)
                # ice_extended_2 = np.hstack((left_column, extended, right_column))

                # # Create meshgrid for interpolation
                # ice_long_grid, ice_lat_grid = np.meshgrid(long_extended, lat_extended)
                
                # # Flatten for griddata interpolation
                # points = np.array([ice_long_grid.ravel(), ice_lat_grid.ravel()]).T
                # values = ice_extended_2.ravel()

                # # Interpolation on common grid
                # lon_out_flat = grided_lon.flatten()
                # lat_out_flat = grided_lat.flatten()
                # grid_points = np.array([lon_out_flat, lat_out_flat]).T
                # interp_flat = griddata(points, values, grid_points, method='linear',rescale=True) # Use RegularGridInterpolator

                # # Reshape interpolated data back into grid
                # interp = interp_flat.reshape(grided_lon.shape)

                grided_lon, grided_lat = np.meshgrid(self.elons, self.lats)
                nointerp = grd.copy()[::-1,:]

                # Extend latitude and longitude arrays
                lat_extended = np.concatenate(([90], lat[::-1], [-90]))
                lon_extended = np.concatenate(([-0.5], lon, [360.5]))

                # Extend data array
                top_row = np.zeros((1, len(lon)))
                bottom_row = nointerp[0, -1] * np.ones((1, len(lon)))
                extended = np.vstack((top_row, nointerp, bottom_row))

                left_column = extended[:, -1].reshape(-1, 1)
                right_column = extended[:, 0].reshape(-1, 1)
                ice_extended_2 = np.hstack((left_column, extended, right_column))

                # Create RegularGridInterpolator
                interp_func = RegularGridInterpolator(
                    (lat_extended, lon_extended),
                    ice_extended_2,
                    method='linear',
                    bounds_error=False,
                    fill_value=None
                )

                # Prepare points for interpolation
                points = np.array([grided_lat.flatten(), grided_lon.flatten()]).T

                # Interpolate
                interp_flat = interp_func(points)

                # Reshape interpolated data back into grid
                interp = interp_flat.reshape(grided_lat.shape)


                # Store results
                return interp[::-1,:]

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
        
        # def along_transect(self,coord=('lat_start','lon_start','lat_stop','lon_stop'),point_density=None,point_distance=None):
        #     if not(point_density is None):
        #         theta=np.linspace(coord[0],coord[2],point_density)
        #         phi=np.linspace(coord[1],coord[3],point_density)
        #     if not(point_distance is None):
        #         theta=np.arange(coord[0],coord[2],point_distance)
        #         phi=np.arange(coord[1],coord[3],point_distance)
        #     i=0
        #     Y_lm=np.expand_dims(pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[0]+pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[1]*1j,axis=1)
            
        #     for i in range(1,len(phi)):
        #         Y_lm=np.concatenate((Y_lm,np.expand_dims((pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[0]+pysh.expand.spharm(self.maxdeg-1,theta[i],phi[i],packed=True)[1]*1j),axis=1)),axis=1)
        #     return (Y_lm*np.repeat(np.expand_dims(self.coeff,axis=1),point_density,axis=1)).sum(0)

class TIME_GRID(GRID,sphericalobject):
    """
    The TIME_GRID object is used to store the data of a grid over time. This object inherits the GRID object to use several method like the conversion from coeff to grd.
    """

    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=np.array([0]),grid_name='time_grid',from_file=(False,),superinit=False,P_lm=None):
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
            `point_time`_ :
                Compute over time the value of the height_time_coeff component at a group of points
            `along_transect`_ :
                Compute the value at one time step for the grid along a cross section
            `save`_ :
                Save the grid in a specified nc file with the name of the grid

        """
        if not(from_file[0]) :
            self.maxdeg=maxdeg
            super().__init__(P_lm=P_lm)
            self.isgrd=False
            self.iscoeff=False


            self.time_grid_name=grid_name

            self.saved=np.array([]) # initialize the save of the spherical harmonic object

            if len(rho)==1:
                self.rho=np.repeat(rho,len(time_step))
            else :
                self.rho=rho.copy()

            # print(self.rho)

            self.time_step=time_step
            self.time_step_number=len(time_step)

            self.height_time_grid=np.zeros((self.time_step_number,self.nlats,self.nlons))
            self.mass_time_grid=np.zeros((self.time_step_number,self.nlats,self.nlons))
            self.height_time_coeff=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j
            self.mass_time_coeff=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j

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
            super().__init__(P_lm=P_lm)

            self.time_step=self.ncgrid['time'][:].data
            self.time_step_number=len(self.time_step)
            self.maxdeg=self.ncgrid.dimensions['maxdeg'].size
            # print(self.maxdeg)

            self.rho=self.ncgrid['rho'][:].data

            self.height_time_grid=self.ncgrid['thickness'][:].data
            self.mass_time_grid=np.zeros((self.time_step_number,self.nlats,self.nlons))
            self.height_time_coeff=self.ncgrid['coeff_real'][:].data+self.ncgrid['coeff_imag'][:].data*1j
            self.mass_time_coeff=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j

            # self.mass_time_grid=self.height_time_grid*rho
            # print(self.rho)

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
        grid_0=grid_to_interp[0]
        grid_to_interp=np.diff(grid_to_interp,axis=0)

        t_step, n, m = grid_to_interp.shape
        new_t_step = len(model_time_step)

        # Create a new grid for the temporal variation redistribution
        delta_new_grid = np.zeros((new_t_step-1, n, m))
        rho_size_right= self.rho.shape[0]!=new_t_step
        if rho_size_right :
            rho_temp=np.zeros((new_t_step,))
        else :
            rho_temp=self.rho.copy()
        for k in range(new_t_step - 1):
            # Finf the corresponding interval in the model input. 
            time_start_new = model_time_step[k]
            time_end_new = model_time_step[k + 1]
            print('t_step :' + str(k)+'for time :' + str(model_time_step[k]))

            for i in range(t_step):

                time_start_original = grid_time_step[i]
                time_end_original = grid_time_step[i + 1]
                delta_time_original = time_start_original-time_end_original
                
                # Determine the overlap proportion between time interval 
                overlap_start = min(time_start_new, time_start_original)
                overlap_end = max(time_end_new, time_end_original)
                overlap_duration = overlap_start - overlap_end
                # print(overlap_duration)
                # print(time_start_original,time_end_original)
                if rho_size_right :
                    if model_time_step[k]>=time_end_original and model_time_step[k]<time_start_original:
                        rho_temp[k]=self.rho[i]
                        # print(rho_temp[k],i,self.rho[i])
                if overlap_duration > 0:
                    # print(overlap_duration/delta_time_original)
                    proportion = overlap_duration / delta_time_original
                    delta_new_grid[k,:,:] += grid_to_interp[i, :, :] * proportion
        rho_temp[-1]=self.rho[-1]
        self.rho=rho_temp.copy()
                    

        return np.concatenate((grid_0[np.newaxis,:,:],delta_new_grid),axis=0).cumsum(0) # return the cumulative grid.

    def interp_on_time_and_space(self,grid_to_interp,grid_time_step,grid_lon,grid_lat,interp_type='Thickness_divide',backend=False,grid_type='global',interp_method='stripy'):
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
        # time_grid_interp_0=self.interp_on(grid_to_interp[0],grid_lon,grid_lat,grid_type=grid_type)
        # grid_to_interp=np.diff(grid_to_interp,axis=0)
        time_grid_pre_interp=np.zeros((len(grid_time_step),self.nlats,self.nlons))
        for i in range(len(grid_time_step)):
            if backend :
                print('interpolation number : ' + str(i))
            time_grid_pre_interp[i,:,:]=self.interp_on(grid_to_interp[i],grid_lon,grid_lat,grid_type=grid_type,interp_method=interp_method)
        grid_type='global'
        # print(time_grid_pre_interp.shape)
        # time_grid_pre_interp=np.concatenate((time_grid_interp_0[np.newaxis,:,:],time_grid_pre_interp),axis=0).cumsum(0) # sum the time step to obtain the total thickness grid

        self.height_time_grid=time_grid_pre_interp.copy()

        self.height_time_grid=self.interp_on_time(time_grid_pre_interp,grid_time_step,self.time_step,interp_type,backend=backend,grid_type=grid_type)
    
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
        for i in range(self.time_step_number):
            self.grd=self.height_time_grid[i,:,:]
            print(self.height_time_coeff.shape,self.grdtocoeff().coeff.shape)
            self.height_time_coeff[i,]=self.grdtocoeff().coeff
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
        for i in range(self.time_step_number):
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

    # def plot_step_on_sphere(self,time_step,cmap=cm.inferno,vmin=None,vmax=None,clip=False,inverse=False):
    #     colormap=cmap
    #     if vmin==None and vmax==None :
    #         normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
    #     elif vmax==None and not(vmin==None):
    #         normaliser = mpl.colors.Normalize(vmin=vmin, vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
    #     elif vmin==None and not(vmax==None):
    #         normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=vmax,clip=clip)
    #     if inverse:
    #         normaliser.inverse()
        
    #     elons,lats=np.meshgrid(self.elons,self.lats)
        
    #     u=elons/360*2*math.pi
    #     v=(lats+90)/180*math.pi
    #     x,y,z=sph2cart(u.flatten(),v.flatten(),np.ones((u.flatten().shape)))
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection="3d")
    #     ax.plot_surface(np.reshape(x,(self.lats.shape[0],self.elons.shape[0])),np.reshape(y,(self.lats.shape[0],self.elons.shape[0])),np.reshape(z,(self.lats.shape[0],self.elons.shape[0])),facecolors=colormap(normaliser(self.height_time_grid[time_step,:,:])),cmap=colormap)
    #     ax.set_aspect('equal')

    # def scatter_step_on_sphere(self,time_step,cmap=cm.inferno,vmin=None,vmax=None,clip=False,inverse=False,marker='.',s=0.2):
    #     colormap=cmap
    #     if vmin==None and vmax==None :
    #         normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
    #     elif vmax==None and not(vmin==None):
    #         normaliser = mpl.colors.Normalize(vmin=vmin, vmax=np.max(self.height_time_grid[time_step,:,:]),clip=clip)
    #     elif vmin==None and not(vmax==None):
    #         normaliser = mpl.colors.Normalize(vmin=np.min(self.height_time_grid[time_step,:,:]), vmax=vmax,clip=clip)
    #     if inverse:
    #         normaliser.inverse()

    #     elons,lats=np.meshgrid(self.elons,self.lats)
        
    #     u=elons/360*2*math.pi
    #     v=(lats+90)/180*math.pi
    #     x,y,z=sph2cart(u.flatten(),v.flatten(),np.ones((u.flatten().shape)))
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection="3d")
    #     ax.scatter(x,y,z,marker=marker,c=colormap(normaliser(self.height_time_grid[time_step,:,:].flatten())),s=s)
    #     ax.set_aspect('equal')

    # def plot_step(self,time_step):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     elons,lats=np.meshgrid(self.elons,self.lats)
    #     ax.pcolor(elons,lats,self.height_time_grid[time_step,:,:])

    def point_time(self,coords):
        """
        The _`point_time` function is used compute a value at one point over each time step for a set of spherical harmonics coeffisient in self.coeff. Parameters saved are, longitude, latitude, maximum degree, the time steps of the grid, the thickness of the grid, the harmonic coefficient, grid density. The harmonic coefficient, due to complexe data type management of nc, are saved separately in there complexe and real part. The created file will have the name of the grid. 

        Attribute :
        -----------
            coords : np.array([n,2])
                Coordinates at wich the spherical harmonics coefficient are converted to real values. The shape of the coordinate must be [latitude, longitude]. 
        Return :
        --------
            points : np.array([n,self.time_step_number])
                The output values for each time steps, 
        """
        points=np.zeros((self.time_step_number))
        for t_it in range(self.time_step_number):
            self.coeff_from_step(t_it)
            points[t_it]=self.calc_at_point(coords[0],coords[1])
            # coeff=self.height_time_coeff[t_it,:].copy()
            # coeff=np.stack((coeff.real,coeff.imag))
            # coeff=pysh.shio.SHCindexToCilm(coeff).copy()
            # points[:,t_it]=-pysh.expand.MakeGridPoint(coeff,coords[:,0],coords[:,1]).copy()
        return points
    
    def along_transect(self,coord=('lat_start','lon_start','lat_stop','lon_stop'),point_density=None,point_distance=None,backend=False):
        """
        The _`along_transect` function is used compute the values of the grid along a transect.

        Attribute :
        -----------
            coord : tuple(4)
                Coordinates of the starting and final point of the transect. 
            point_density : int 
                The number of points along the transec
            point_distance : float
                The distance between points along the transect
            backend : bool 
                If you want a backend
        Return :
        --------
            points : np.array([n,self.time_step_number])
                The output values for each time steps, 
        """
        if not(point_density is None):
            theta=np.linspace(-coord[0],-coord[2],point_density)
            phi=np.linspace(coord[1],coord[3],point_density)
        if not(point_distance is None):
            theta=np.arange(-coord[0],-coord[2],point_distance)
            phi=np.arange(coord[1],coord[3],point_distance)
        coeff=np.stack((self.coeff.real,self.coeff.imag))
        coeff=pysh.shio.SHCindexToCilm(coeff)
        transect=pysh.expand.MakeGridPoint(coeff,theta,phi)
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
        self.ncgrid.createDimension('maxdeg_order',(self.maxdeg+1)*(self.maxdeg+2)/2)
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

        if self.height_time_grid.shape[0]==self.time_step_number :
            thickness=self.ncgrid.createVariable('thickness',np.float32,('time_step','lat','lon'))
        else : 
            thickness=self.ncgrid.createVariable('thickness',np.float32,('time_diff','lat','lon'))

        thickness.units='m'
        thickness.long_name='layer_thickness'
        thickness[:,:,:]=self.height_time_grid


        if self.height_time_coeff.shape[0]==self.time_step_number :
            coeff_real=self.ncgrid.createVariable('coeff_real',np.float32,('time_step','maxdeg_order'))
        else : 
            coeff_real=self.ncgrid.createVariable('coeff_real',np.float32,('time_diff','maxdeg_order'))
        coeff_real.units='m'
        coeff_real.long_name='layer thickness in spherical harmonics'
        coeff_real[:,:]=np.real(self.height_time_coeff)

        
        if self.height_time_coeff.shape[0]==self.time_step_number :
            coeff_imag=self.ncgrid.createVariable('coeff_imag',np.float32,('time_step','maxdeg_order'))
        else : 
            coeff_imag=self.ncgrid.createVariable('coeff_imag',np.float32,('time_diff','maxdeg_order'))
        coeff_imag.units='m'
        coeff_imag.long_name='layer thickness in spherical harmonics'
        coeff_imag[:,:]=np.imag(self.height_time_coeff)

        if len(self.rho)>1 :
            rho=self.ncgrid.createVariable('rho',np.float32,('time_step'))
        else :
            rho=self.ncgrid.createVariable('rho',np.float32)
        rho.units='kg/m3'
        rho.long_name='density of the considered grid'
        print(self.rho)
        rho[:]=self.rho
        
        if not(supersave) :
            self.ncgrid.close()

class SEDIMENT_TIME_GRID(TIME_GRID):
    """
    The SEDIMENT_TIME_GRID object is used to contain the data about the sediment grid.
    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=np.array([2300]),grid_name='time_grid',from_file=(False,),P_lm=None) : 
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
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True,P_lm=P_lm)
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        self.sed=np.zeros((self.time_step_number,self.maxdeg,self.maxdeg*2))
        if from_file[0] :
            self.sed=self.ncgrid['sed'][:].data
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
        sed=self.ncgrid.createVariable('sed',np.float32,('time_step','lat','lon'))
        sed.units='m'
        sed.long_name='initial_sed_thickness'
        sed[:,:,:]=self.sed

        self.ncgrid.close()

class ICE_TIME_GRID(TIME_GRID):
    """
    The ICE_TIME_GRID object is used to store the data about the ice and the potential modification of the model.
    """    
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=np.array([920]),grid_name='time_grid',from_file=(False,),P_lm=None) : 
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
        # initialize to false the grid and coefficient (no grid or coefficient pre loaded). The ice volumetric masse is set to 916.7 kg/m3.
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True,P_lm=P_lm)
        self.saved=np.array([])
        self.sdeli_00=0
        self.deli_00_prev=0
        self.ice=np.zeros((self.time_step_number,self.maxdeg,self.maxdeg*2))
        self.corrected_ice=np.zeros((self.time_step_number,self.maxdeg,self.maxdeg*2))

        if from_file[0] :
            self.ice=self.ncgrid['ice'][:].data
            self.corrected_ice=self.ncgrid['corrected_ice'][:].data
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
        for t_it in range(self.time_step_number):
            topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]-self.corrected_ice[t_it,:,:]+self.ice[t_it,:,:]
        for t_it in range(self.time_step_number): 
            check1 = OCEAN_TIME_GRID().evaluate_ocean(-topo.height_time_grid[t_it,:,:]+self.ice[t_it,:,:]) # generate the ocean function for ice-topo
            check2 = OCEAN_TIME_GRID().evaluate_ocean(topo.height_time_grid[t_it,:,:]-(self.ice[t_it,:,:])).grd*(OCEAN_TIME_GRID().evaluate_ocean(-self.ice[t_it,:,:]*self.rho[t_it]-(topo.height_time_grid[t_it,:,:]-self.ice[t_it,:,:])*oc.rho[t_it]).grd)# add the two part of ice over the check1 nd check2 positive area.
            self.corrected_ice[t_it,:,:] =  check1.grd*(self.ice[t_it,:,:])+check2*(self.ice[t_it,:,:])# add the two part of ice over the check1 nd check2 positive area.
        self.height_time_grid=np.concatenate((np.zeros(self.corrected_ice[0,:,:].shape)[np.newaxis,:,:],np.diff(self.corrected_ice,axis=0)),axis=0)
        for t_it in range(self.time_step_number): 
            topo.height_time_grid[t_it,:,:]=topo.height_time_grid[t_it,:,:]+self.corrected_ice[t_it,:,:]-self.ice[t_it,:,:]
        return topo

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
        ice=self.ncgrid.createVariable('ice',np.float32,('time_step','lat','lon'))
        ice.units='m'
        ice.long_name='initial_ice_thickness'
        ice[:,:,:]=self.ice

        corrected_ice=self.ncgrid.createVariable('corrected_ice',np.float32,('time_step','lat','lon'))
        corrected_ice.units='m'
        corrected_ice.long_name='initial_ice_thickness'
        corrected_ice[:,:,:]=self.corrected_ice

        self.ncgrid.close()

class OCEAN_TIME_GRID(TIME_GRID):
    """
    The OCEAN_TIME_GRID object is used to store and manage the ocean model. This object is also used to resolve the SLE.
    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=np.array([1000]),grid_name='time_grid',from_file=(False,),P_lm=None) :
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
            `sea_level_solver`_ 
                apply the convergence toward the sea level solution for a particular time step
            `sea_level_equation`_ 
                compute the sea level equation resolution for one iteration

        """
        # initialize the ocean with no grid and no coefficient. The saved is also initialized. The volumetric mass of water is set to 1000 
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True,P_lm=P_lm)
        self.S=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j
        self.RO=np.zeros((self.time_step_number,int(maxdeg),int(maxdeg*2)))
        self.TO=np.zeros((self.time_step_number,int(maxdeg),int(maxdeg*2)))
        self.ESL=np.zeros((self.time_step_number,))
        self.saved=np.array([])
        if from_file[0] :
            self.S=self.ncgrid['S_real'][:].data+self.ncgrid['S_imag'][:].data*1j
            self.RO=self.ncgrid['RO'][:].data
            self.TO=self.ncgrid['TO'][:].data
            self.ESL=self.ncgrid['ESL'][:].data
            self.ncgrid.close()
        
    def update_0(self,type='grd'):
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
        if type is 'grd' :
            self.grd_0=self.grd.copy()
        else :
            self.coeff_0=self.coeff.copy()
        
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
        self.grd=out.copy()
        return self
    
    def sea_level_solver(self,load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,conv_lim,topo_it):
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
        chi=np.inf
        while chi>=conv_lim and conv_it<=100:
            chi=self.sea_level_equation(load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,topo_it)
            conv_it+=1
        return conv_it
    
    def sea_level_equation(self,load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,topo_it):
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
        if conv_it==0 and topo_it==0 :
            print(int(self.maxdeg*(self.maxdeg+1)/2))
            self.height_time_coeff[t_it,:] =self.prev/self.prev[0]*(-ice_time_grid.rho[t_it]/self.rho[t_it]*ice_time_grid.sdeli_00+TO.coeff[0]-TO.prev[0])-TO.coeff-TO.prev
            # print(TO.coeff[0])
            # print(self.prev[0])
        self.S[t_it,:]=self.S[t_it-1,:]+self.height_time_coeff[t_it,:]
        if t_it == 1 : 
            load.calc_rot_visc(ice_time_grid.height_time_coeff[1,:]*ice_time_grid.rho[1]+sed_time_grid.height_time_coeff[1,:]*sed_time_grid.rho[1]+self.height_time_coeff[1,:]*self.rho[1],t_it,love_number)
            load.sdelLa[t_it-1]=load.rot_pot[t_it,:]-load.rot_pot[t_it-1,:]
            load.V_lm_tide.coeff=np.zeros(6)
            load.V_lm.coeff=np.zeros(ice_time_grid.height_time_coeff[0,:].shape)
        else : 
            load.calc_rot_visc(ice_time_grid.height_time_coeff[1:t_it+1,:]*np.repeat(ice_time_grid.rho[1:t_it+1,np.newaxis],self.N,axis=1)+sed_time_grid.height_time_coeff[1:t_it+1,:]*np.repeat(sed_time_grid.rho[1:t_it+1,np.newaxis],self.N,axis=1)+self.height_time_coeff[1:t_it+1,:]*np.repeat(self.rho[1:t_it+1,np.newaxis],self.N,axis=1),t_it,love_number)
            load.sdelLa[t_it-1]=load.rot_pot[t_it,:]-load.rot_pot[t_it-1,:]
            load.calc_viscuous(ice_time_grid.height_time_coeff[1:t_it,:]*np.repeat(ice_time_grid.rho[1:t_it,np.newaxis],self.N,axis=1)+sed_time_grid.height_time_coeff[1:t_it,:]*np.repeat(sed_time_grid.rho[1:t_it,np.newaxis],self.N,axis=1)+self.height_time_coeff[1:t_it,:]*np.repeat(self.rho[1:t_it,np.newaxis],self.N,axis=1),love_number.beta_l,t_it)
            load.calc_rotational_viscuous(load.sdelLa[:t_it-1,:],love_number.beta_l_tide,t_it)
        # print(load.V_lm_tide.coeff)
        # print(load.sdelLa[t_it-1])
        delSLcurl_tide_fl=1/load.g*load.V_lm_tide.coeff+1/load.g*love_number.E_T[:6]*load.rot_pot[t_it,:]
        # print(load.V_lm.coeff.shape)
        # print(love_number.E.shape)
       
        delSLcurl_fl=love_number.E* love_number.T.coeff *((ice_time_grid.height_time_coeff[:t_it+1,:]*np.repeat(ice_time_grid.rho[:t_it+1,np.newaxis],self.N,axis=1)).sum(0)+(sed_time_grid.height_time_coeff[:t_it+1,:]*np.repeat(sed_time_grid.rho[:t_it+1,np.newaxis],self.N,axis=1)).sum(0)+(self.height_time_coeff[:t_it+1,:]*np.repeat(self.rho[:t_it+1,np.newaxis],self.N,axis=1)).sum(0))+love_number.T.coeff*load.V_lm.coeff
        # print(conv_it)
        # # print(load.sdelLa[t_it-1])
        # print(load.rot_pot[t_it])
        # print(load.sdelI[t_it-1])
        # print(load.sdelm[t_it-1])
        # print(ice_time_grid.height_time_coeff[:t_it+1,0].sum(),self.S[t_it,0])
        # print(ice_time_grid.height_time_coeff[:t_it+1,0].sum()*ice_time_grid.rho,self.S[t_it,0]*self.rho[0])
        # print(ice_time_grid.height_time_coeff[:t_it+1,0].sum()*ice_time_grid.rho+self.S[t_it,0]*self.rho[0])
        # print((ice_time_grid.height_time_coeff[:t_it+1,:].sum(0)*ice_time_grid.rho+sed_time_grid.height_time_coeff[:t_it+1,:].sum(0)*sed_time_grid.rho[t_it]+self.S[t_it,:]*self.rho[0])[:6])
        # print((love_number.E* love_number.T.coeff *(ice_time_grid.height_time_coeff[:t_it+1,:].sum(0)*ice_time_grid.rho+sed_time_grid.height_time_coeff[:t_it+1,:].sum(0)*sed_time_grid.rho[t_it]+self.S[t_it,:]*self.rho))[:6])
        # print((love_number.T.coeff*load.V_lm.coeff)[:6])
        
        # print(1/load.g*love_number.E_T[:6]*load.rot_pot[t_it,:])
        # print(1/load.g*load.V_lm_tide.coeff)
        delSLcurl_tide_fl=np.concatenate((delSLcurl_tide_fl,np.zeros(delSLcurl_fl.shape[0]-delSLcurl_tide_fl.shape[0])))
        # print(sphericalobject(coeff=delSLcurl_fl + delSLcurl_tide_fl).coeff)
        # print((delSLcurl_fl + delSLcurl_tide_fl)[:6])
        self.delSLcurl=sphericalobject(coeff=delSLcurl_fl + delSLcurl_tide_fl,P_lm=ice_time_grid.P_lm).coefftogrd().grd - ice_time_grid.height_time_grid[1:t_it+1,:,:].sum(0)- sed_time_grid.height_time_grid[1:t_it+1,:,:].sum(0)

        self.RO[t_it,:,:]=self.delSLcurl*self.grd
        RO=sphericalobject(grd=self.RO[t_it,:,:].copy(),P_lm=ice_time_grid.P_lm).grdtocoeff()
        self.delPhi_g=np.real(1/self.coeff[0] * (- ice_time_grid.rho[t_it]/self.rho[t_it]*ice_time_grid.height_time_coeff[:t_it+1,0].sum(0) - RO.coeff[0] + TO.coeff[0]))

        sdelS_new=RO.coeff + self.delPhi_g*self.coeff -  TO.coeff - self.S[t_it-1,:]

        chi = np.abs((np.sum(np.abs(sdelS_new)) - np.sum(np.abs(self.height_time_coeff[t_it,:]))) / np.sum(np.abs(self.height_time_coeff[t_it,:])))
        # print(chi)
        self.height_time_coeff[t_it,:]=sdelS_new.copy()
        return chi
    
    def calculate_ESL(self,t_it,ice,topo,type='delPhi_g'):
        '''
        The _`calculate_dESL` function compute the ESL variation over time including a varaible ocean surface.
        
        Attribute :
        ----------- 
            Input_way : str
                way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
            type : str 
                Define the type of world mean sea level you want calculate. delPhi_g type include shore migration, ESL include only variations of ocean volume due to melting ice.

        Returns : 
        --------- 
            dESL : np.array([time_step_number,])
                The variation of the ESL at each time step.
            
        '''
        if t_it==0 :
            ESL=0
        else :
            if type=='delPhi_g' :
                self.evaluate_ocean(topo.height_time_grid[t_it,:,:]).grdtocoeff()
                TO=sphericalobject(grd=self.TO[t_it,:,:],P_lm=ice.P_lm).grdtocoeff().coeff.copy()
                RO=sphericalobject(grd=self.RO[t_it,:,:],P_lm=ice.P_lm).grdtocoeff().coeff.copy()
                self.ESL=np.real(1/self.coeff[0] * (- ice.rho[0]/self.rho[0]*ice.height_time_coeff[:t_it+1,0].sum(0)-RO[0]+TO[0]))
            elif type=='ESL' :
                self.evaluate_ocean(topo.height_time_grid[-1,:,:]).grdtocoeff()
                self.ESL=np.real(1/self.coeff[0] * (- ice.rho[0]/self.rho[0]*ice.height_time_coeff[:t_it+1,0].sum(0)))
        # return ESL
    
    def calculate_ESL_time(self,ice,topo,type='delPhi_g'):
        '''
        The _`calculate_dESL` function compute the ESL variation over time including a varaible ocean surface.
        
        Attribute :
        ----------- 
            Input_way : str
                way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
            type : str 
                Define the type of world mean sea level you want calculate. delPhi_g type include shore migration, ESL include only variations of ocean volume due to melting ice.

        Returns : 
        --------- 
            dESL : np.array([time_step_number,])
                The variation of the ESL at each time step.
            
        '''
        ESL=np.zeros(ice.time_step_number)
        for t_it in range(ice.time_step_number):
            ESL[t_it]=self.calculate_ESL(t_it,ice,topo,type=type)
        return ESL

    
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
        S_real=self.ncgrid.createVariable('S_real',np.float32,('time_step','maxdeg_order'))
        S_real.units='none'
        S_real.long_name='Sea level variation'
        S_real[:,:]=np.real(self.S)

        S_imag=self.ncgrid.createVariable('S_imag',np.float32,('time_step','maxdeg_order'))
        S_imag.units='none'
        S_imag.long_name='Sea level variation'
        S_imag[:,:]=np.imag(self.S)

        TO=self.ncgrid.createVariable('TO',np.float32,('time_step','lat','lon'))
        TO.units='none'
        TO.long_name='Sea level variation'
        TO[:,:,:]=self.TO

        RO=self.ncgrid.createVariable('RO',np.float32,('time_step','lat','lon'))
        RO.units='none'
        RO.long_name='Sea level variation'
        RO[:,:,:]=self.RO

        ESL=self.ncgrid.createVariable('ESL',np.float32,('time_step'))
        ESL.units='m'
        ESL.long_name='Eustatic sea level variations'
        ESL[:]=self.ESL


class TOPOGRAPHIC_TIME_GRID(TIME_GRID):
    """
    The TOPOGRAPHIC_TIME_GRID is used to manage the topographic data over time.
    """
    def __init__(self,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=np.array([0]),grid_name='time_grid',from_file=(False,),it_max=10,P_lm=None) : 
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
        super().__init__(time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True,P_lm=P_lm)
        self.it_max=it_max
        self.topo_pres=np.zeros((self.maxdeg,self.maxdeg*2))
        self.topo_initial=np.zeros((self.it_max+1,self.maxdeg,self.maxdeg*2))
        if from_file[0] :
            self.topo_initial=self.ncgrid['topo_initial'][:].data
            self.topo_pres=self.ncgrid['topo_pres'][:].data
            # print(self.topo_initial.shape)
            self.it_max=self.topo_initial.shape[0]-1
            self.ncgrid.close()
        else : 
            self.height_time_grid=np.zeros((self.time_step_number,self.maxdeg,self.maxdeg*2))
            self.height_time_coeff=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j

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
        topo_pres[:,:]=self.topo_pres

        self.ncgrid.createDimension('it_max',self.it_max+1)

        topo_initial=self.ncgrid.createVariable('topo_initial',np.float32,('it_max','lat','lon'))
        topo_initial.units='m'
        topo_initial.long_name='present topography'
        topo_initial[:,:,:]=self.topo_initial

        self.ncgrid.close()


from .love import get_tlm

class LOAD_TIME_GRID(TIME_GRID) :
    """
    The LOAD_ITME_GRID object is used to manage and compute the deformation of the solid earth and geoid. 
    """
    def __init__(self,sdelL=np.array([]),beta_l=np.array([]),beta_l_tide=np.array([]),E=np.array([]),E_T=np.array([]),a=7371000,Me=5.9742e24,time_step=np.array([1,2]),maxdeg=64,height_time_grid=None,mass_time_grid=None,mass_time_coeff=None,height_time_coeff=None,rho=np.array([0]),grid_name='time_grid',from_file=(False,),g=9.80665,P_lm=None):
        """
        The _`LOAD_TIME_GRID` class is used to save and include all the topographic variations. This class inherits of :ref:`TIME_GRID <TIME_GRID>` and :ref:`LOAD <LOAD>`.

        Attributes
        ---------- 
            sdelL : np.array([time_step_number,maxdeg,maxdegx2])
                The load variation grid used to compute earth vertical motion.
            betal : np.array([time_step_number,time_step_number,maxdeg])
                The beta love number as described in :ref:`Variation of geoïd and ground Equations <geoid_ground_variation_theory>` section. There calculated in the :ref:`LOVE <LOVE>` class.
            betal_tide : np.array([time_step_number,time_step_number,maxdeg])
                The tidal beta love number as described in :ref:`Variation of geoïd and ground Equations <geoid_ground_variation_theory>` section. There calculated in the :ref:`LOVE <LOVE>` class.
            E : np.array([(maxdeg+1)(maxdeg+2)/2,])
                The elastic component of the earth as love numbers computed form :ref:`LOVE <LOVE>` class.
            E_T : np.array([(maxdeg+1)(maxdeg+2)/2,])
                The tidal elastic component of the earth as love numbers computed form :ref:`LOVE <LOVE>` class.
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
            `calc_rotational_viscuous`_:
                Compute tje rotational viscous motion of the geoid and grounf for one time step.
            `calc_viscuous_time`_ :
                Compute the viscous ground motion of the earth on all time steps.    
            `calc_tidal_viscuous_time`_:
                Compute the tidal viscuous ground motion of the reath on all time steps
            `calc_elastic_time`_ :
                Compute the elastic ground motion of the earth on all time steps.
            :ref:`save <load_save>` :
                Save the load data

        """
        TIME_GRID.__init__(self,time_step,maxdeg,height_time_grid,mass_time_grid,height_time_coeff,mass_time_coeff,rho,grid_name,from_file,superinit=True,P_lm=P_lm)
        calc_beta_counter(self,maxdeg+1)
        if from_file[0]:
            # self.a=self.ncgrid['a'][:].data
            # self.Me=self.ncgrid['Me'][:].data
            self.viscuous_deformation=self.ncgrid['viscuous_deformation_real'][:].data+self.ncgrid['viscuous_deformation_imag'][:].data*1j
            self.elastic_deformation=self.ncgrid['elastic_deformation_real'][:].data+self.ncgrid['elastic_deformation_imag'][:].data*1j
            self.viscuous_tidal_deformation=self.ncgrid['viscuous_tidal_deformation_real'][:].data+self.ncgrid['viscuous_tidal_deformation_imag'][:].data*1j
            self.elastic_tidal_deformation=self.ncgrid['elastic_tidal_deformation_real'][:].data+self.ncgrid['elastic_tidal_deformation_imag'][:].data*1j
            self.elastic_love=self.ncgrid['elastic_love'][:].data
            self.load=self.ncgrid['load_real'][:].data+self.ncgrid['load_imag'][:].data*1j
            self.a=a
            self.Me=Me
        else :
            self.g=g
            self.viscuous_deformation=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j
            self.elastic_deformation=np.zeros((self.time_step_number,int((maxdeg+1)*(maxdeg+2)/2)))+0j
            self.viscuous_tidal_deformation=np.zeros((self.time_step_number,6))+0j
            self.elastic_tidal_deformation=np.zeros((self.time_step_number,6))+0j
            self.load=sdelL
            self.elastic_love=E
            self.elastic_tidal_love=E_T
            self.a=a
            self.Me=Me
            self.beta_l=beta_l
            self.beta_l_tide=beta_l_tide
            self.viscuous_love=self.beta_l[:,:,self.beta_counter.astype(int)]
            self.viscuous_tidal_love=self.beta_l_tide[:,:,self.beta_counter[:6].astype(int)]
            self.V_lm=sphericalobject(coeff=np.zeros((int((maxdeg+1)*(maxdeg+2)/2)))+0j)
            self.V_lm_tide=sphericalobject(coeff=np.zeros((6,))+0j)
            self.rot_pot=np.zeros((self.time_step_number,6))+0j
            self.sdelI=np.zeros((self.time_step_number,3))+0j
            self.sdelm=np.zeros((self.time_step_number,3))+0j
            self.sdelLa=np.zeros((self.time_step_number,6))+0j
        self.T = sphericalobject(coeff=get_tlm(self.maxdeg,self.a,self.Me))

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
        self.V_lm.coeff = np.einsum('ij,ij->j', beta[t_it-1,:t_it-1], sdelL)

    def calc_rotational_viscuous(self,rot_pot,beta_tide,t_it):
        '''
        The _`calc_rotational_viscuous` method is used to calculate the ground and geoïd deformation based on viscuous love numbers for the rotational component.

        Attribute
        ---------
            sdelL : np.array([t_it,(maxdeg+1)(maxdeg+20/2)])
                The load grid used to estimate the ground vertical mouvement. This include all previous loading history because of the viscous comportment of earth. 
            beta : np.array([time_step_number,time_step_number,(maxdeg+1)(maxdeg+2)/2])
                The beta love numbers used to compute the earth deformation to include the viscous part. These love numbers are particularly heavy in the memory due to the representation of the time. 
            t_it : int
                The time iteration at wich the computation is performed.
        
        '''
        self.V_lm_tide.coeff=np.einsum('ij,ij->j', beta_tide[t_it-1,:t_it-1,:6], rot_pot)

            

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
        for t_it in range(1,self.time_step_number):
            if t_it==1 :
                self.viscuous_deformation[t_it,:]=np.zeros(self.viscuous_deformation[t_it,:].shape)
            else :
                self.calc_viscuous(self.load[1:t_it,:],self.viscuous_love,t_it)
                self.viscuous_deformation[t_it,:]=self.T.coeff*np.squeeze(self.V_lm.coeff.T)
            if backend:
                print(f'viscuous calculation at {self.time_step[t_it]} kyr done')

    def calc_tidal_viscuous_time(self,backend=False):
        '''
        The _`calc_tidal_viscuous_time` method compute the vicuous vertical ground motion. This method call the :ref:`LOAD <LOAD>` method for this. 

        Attribute
        ---------
            backend : bool
                Specifie if the method give backend (True) or not (False). Default is False.

        Return
        ------
            None

        '''
        for t_it in range(1,self.time_step_number):
            if t_it==1 :
                self.viscuous_tidal_deformation[t_it,:]=np.zeros(self.viscuous_tidal_deformation[t_it,:].shape)
            else :
                self.calc_rotational_viscuous(self.sdelLa[:t_it-1,:],self.viscuous_tidal_love,t_it)
                self.viscuous_tidal_deformation[t_it,:]=1/self.g*np.squeeze(self.V_lm_tide.coeff.T)
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
        self.elastic_deformation=np.repeat(np.expand_dims(self.T.coeff,axis=0),self.time_step_number,axis=0)*np.repeat(np.expand_dims(self.elastic_love,axis=0),self.time_step_number,axis=0)*self.load.cumsum(0).squeeze()
 
    def calc_tidal_elastic_time(self):
        '''
        The _`calc_tidal_elastic_time` method compute the elactic vertical ground motion due to the earth rotation. This method call the :ref:`LOAD <LOAD>` method for this. 

        Attribute
        ---------
            None

        Return
        ------
            None
             
        '''
        # print(self.time_grid_name)
        # print(np.repeat(np.expand_dims(self.T.coeff[:6],axis=0),self.time_step_number,axis=0)[-1])
        # print(np.repeat(np.expand_dims(self.elastic_tidal_love[:6],axis=0),self.time_step_number,axis=0)[-1])
        # print(self.rot_pot.squeeze()[-1])
        self.elastic_tidal_deformation=1/self.g*np.repeat(np.expand_dims(self.elastic_tidal_love[:6],axis=0),self.time_step_number,axis=0)*self.rot_pot.squeeze()
        # print(self.elastic_tidal_deformation[-1,:])
        # self.elastic_tidal_deformation=1/self.g*np.concatenate(((self.T.coeff[:6]*self.elastic_tidal_love[:6]*self.rot_pot[0,:6])[np.newaxis,:],np.diff(np.repeat(np.expand_dims(self.T.coeff[:6],axis=0),self.time_step_number,axis=0)*np.repeat(np.expand_dims(self.elastic_tidal_love[:6],axis=0),self.time_step_number,axis=0)*self.rot_pot[:,:6].squeeze(),axis=0)),axis=0)

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

        self.ncgrid.createDimension('maxdeg_tidal',6)

        load_real=self.ncgrid.createVariable('load_real',np.float32,('time_step','maxdeg_order'))
        load_real.units='kg'
        load_real.long_name='load grid used to compute the earth deformation.'
        load_real[:,:]=np.real(self.load)

        load_imag=self.ncgrid.createVariable('load_imag',np.float32,('time_step','maxdeg_order'))
        load_imag.units='kg'
        load_imag.long_name='load grid used to compute the earth deformation.'
        load_imag[:,:]=np.imag(self.load)

        viscuous_deformation_real=self.ncgrid.createVariable('viscuous_deformation_real',np.float32,('time_step','maxdeg_order'))
        viscuous_deformation_real.units='mm/yr'
        viscuous_deformation_real.long_name='viscuous component of the earth deformation due to load.'
        viscuous_deformation_real[:,:]=np.real(self.viscuous_deformation)

        viscuous_deformation_imag=self.ncgrid.createVariable('viscuous_deformation_imag',np.float32,('time_step','maxdeg_order'))
        viscuous_deformation_imag.units='mm/yr'
        viscuous_deformation_imag.long_name='viscuous component of the earth deformation due to load.'
        viscuous_deformation_imag[:,:]=np.imag(self.viscuous_deformation)

        elastic_deformation_real=self.ncgrid.createVariable('elastic_deformation_real',np.float32,('time_step','maxdeg_order'))
        elastic_deformation_real.units='mm/yr'
        elastic_deformation_real.long_name='elastic component of the earth deformation due to load.'
        elastic_deformation_real[:,:]=np.real(self.elastic_deformation)

        elastic_deformation_imag=self.ncgrid.createVariable('elastic_deformation_imag',np.float32,('time_step','maxdeg_order'))
        elastic_deformation_imag.units='mm/yr'
        elastic_deformation_imag.long_name='elastic component of the earth deformation due to load.'
        elastic_deformation_imag[:,:]=np.imag(self.elastic_deformation)

        viscuous_tidal_deformation_real=self.ncgrid.createVariable('viscuous_tidal_deformation_real',np.float32,('time_step','maxdeg_tidal'))
        viscuous_tidal_deformation_real.units='mm/yr'
        viscuous_tidal_deformation_real.long_name='viscuous component of the earth deformation due to load.'
        viscuous_tidal_deformation_real[:,:]=np.real(self.viscuous_tidal_deformation)

        viscuous_tidal_deformation_imag=self.ncgrid.createVariable('viscuous_tidal_deformation_imag',np.float32,('time_step','maxdeg_tidal'))
        viscuous_tidal_deformation_imag.units='mm/yr'
        viscuous_tidal_deformation_imag.long_name='viscuous component of the earth deformation due to load.'
        viscuous_tidal_deformation_imag[:,:]=np.imag(self.viscuous_tidal_deformation)

        elastic_tidal_deformation_real=self.ncgrid.createVariable('elastic_tidal_deformation_real',np.float32,('time_step','maxdeg_tidal'))
        elastic_tidal_deformation_real.units='mm/yr'
        elastic_tidal_deformation_real.long_name='elastic component of the earth deformation due to load.'
        elastic_tidal_deformation_real[:,:]=np.real(self.elastic_tidal_deformation)

        elastic_tidal_deformation_imag=self.ncgrid.createVariable('elastic_tidal_deformation_imag',np.float32,('time_step','maxdeg_tidal'))
        elastic_tidal_deformation_imag.units='mm/yr'
        elastic_tidal_deformation_imag.long_name='elastic component of the earth deformation due to load.'
        elastic_tidal_deformation_imag[:,:]=np.imag(self.elastic_tidal_deformation)

        elastic_love=self.ncgrid.createVariable('elastic_love',np.float32,('maxdeg_order'))
        elastic_love.units='none'
        elastic_love.long_name='elastic load love numbers used to compute the earth deformation'
        elastic_love[:]=self.elastic_love

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

    def calc_rot_visc_time(self,love_number):
        '''
        The _`calc_rot_visc_time` method apply the computation of the rotational potential at each time steps. 

        Attribute
        ---------
            love_number : :ref:`LOVE <LOVE>` object
                The love numbers object

        Return
        ------
            None
             
        '''
        for t_it in range(1,self.time_step_number):
            if t_it == 1 : 
                self.calc_rot_visc(self.load[1,:],t_it,love_number)
                self.sdelLa[t_it-1]=self.rot_pot[t_it,:]-self.rot_pot[t_it-1,:]
            else : 
                self.calc_rot_visc(self.load[1:t_it+1,:],t_it,love_number)
                self.sdelLa[t_it-1,:]=self.rot_pot[t_it,:]-self.rot_pot[t_it-1,:] # calc small variation


    def calc_rot_visc(self,sdelL,t_it,love_number,G=6.67408E-11,a=6371000,C=8.034e37,k_f=0.942,omega=7.292E-5):
        '''
        The _`calc_rot_visc` method compute the rotational potential induced by load redistribution. 

        Attribute
        ---------
           sdelL : np.array([])
                The small variations in load at earth surface
            t_it : int 
                The time iteration indice
            love_number : :ref:`LOVE <LOVE>` object
                The love numbers object
            G : float
                Constant of gravitation 
            a : float
                The radius of earth
            C : float
        
            k_f : float

            omega : float

        Return
        ------
            None
             
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
        I = np.zeros(3, dtype=complex)
        I[0] = np.sqrt(32/15)*np.pi*(a**4)*np.real(L21)
        I[1] = -np.sqrt(32/15)*np.pi*(a**4)*np.imag(L21)
        I[2] = (8/3)*np.pi*(a**4)*(L00 - L20/np.sqrt(5))

        if t_it==1 :
            V_lm=np.zeros(3)
            V_lm_T=np.zeros(3)
        else :
            V_lm = np.dot(love_number.beta_konly_l[t_it-1,:t_it-1],self.sdelI[:t_it-1,:])
            V_lm_T = np.dot(love_number.beta_konly_l_tide[t_it-1,:t_it-1],self.sdelm[:t_it-1,:])
        # print(V_lm,V_lm_T)
        temp = 1/(1-love_number.k_tide_e[1]/k_f)*(1/CminA * ((1+love_number.k_e[1])*I + V_lm.squeeze()) + V_lm_T.squeeze()/k_f)
        # calculate the perturbation to the rotational potential from Milne 1998
        m1=temp[0].copy()
        m2=temp[1].copy()
        temp = -1/(1-love_number.k_tide_e[1]/k_f)*(1/C * ((1+love_number.k_e[1])*I + V_lm.squeeze()))
        m3=temp[2].copy()

        m=np.array([m1,m2,m3])

        self.sdelI[t_it-1,:] = I.T.squeeze() - self.sdelI[:t_it-1,:].sum(0)
        self.sdelm[t_it-1,:] = m.T.squeeze() - self.sdelm[:t_it-1,:].sum(0)

        self.rot_pot[t_it,0] = (a**2 * omega**2/3 * (np.dot(m, m) + 2*m3)+0j)
        self.rot_pot[t_it,3] = (a**2 * omega**2/(6*np.sqrt(5)) * (m1**2 + m2**2 - 2*m3**2 - 4*m3)+0j)
        self.rot_pot[t_it,4] = (a**2 * omega**2/np.sqrt(30) * (m1*(1+m3) - 1j*m2*(1+m3))+0j)
        self.rot_pot[t_it,5] = (a**2 * omega**2/(np.sqrt(5) * np.sqrt(24)) * ( (m2**2-m1**2) + 1j*2*m1*m2 )+0j)
        