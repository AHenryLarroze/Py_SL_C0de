from spharm import GaussQuad
import numpy as np
import math
from scipy import interpolate
from spharm import sphericalobject
from sediment import spherical_sed
from ice import spherical_ice
from topography import spherical_topo
from love import LOVE
from SeaLevel import spherical_sea_level

class GRID(object):
        """
    A class used to represent the Gaussian Grid

    ...

    Attributes
    ----------
        model_p : Object (from class World_Model_Parameter)
            model parameter used as entry to create the grid
        maxdeg : int
            maximum degree of the spherical harmonic used in the model !!!! peut être ne pas reproduire les constante au travers de objets !!!! 
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
        def __init__(self,model_p):
            """
            Parameters
            ----------
            model_p : object (object from class World_Model_Parameter)
            """
            self.model_p=model_p
            self.maxdeg=model_p.maxdeg # maybe no use of defining the maximum degree here ofr in model_p
            self.nlons = (self.maxdeg) * 2 
            self.nlats = (self.maxdeg)
            x,w=GaussQuad(self.maxdeg)
            x_GL = np.arccos(x)*180/math.pi - 90
            lon_GL = np.linspace(0,360,2*self.maxdeg+1)
            lon_GL = lon_GL[:-1]
            self.lats=x_GL
            #intgrid=pysh.SHGrid.from_zeros(self.maxdeg-1, grid='GLQ') # here we use maxdeg -1 because it is the entry for the shtools model. # we can use the shtools here because it is coincident with SL. 
            #self.lats=intgrid.lats()
            #self.elons  =  intgrid.lons()
            self.elons=lon_GL
            self.colats = 90 - self.lats
            # for the calclation we need elons and colats.
        
        def add_time(self,time_step):
            '''
            
            This function add time field to the object.
        
        Parameters :  
            time_step (np.array): a vector array of length number of time step.
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            time_step (np.array): an array of each time step of length number of time step.
        
            '''
            self.time_step=time_step # à ajouter dans ice grid definition
            
        def interp_on(self,grd,lon,lat):
            '''
            
            This function interpolate a grid of data on the grid of the model.
        
        Parameters :  
            grd (np.array) : The grid data
            lon (np.array) : The longitude of the data
            lat (np.array) : The latitude of the data          
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            grd (np.array) : The interpolated grid with the shape of the gausse legendre grid.
        
            '''
            if grd.shape!=lon.shape :# if the latitude and longitude are not a mesh grid for the grid input then applay meshgrid
                lon,lat=np.meshgrid(lon,lat)
            elons,lats=np.meshgrid(self.elons,self.lats) # mesh grid the elons, elats !!!! for calculation efficiency, it coud be usefull to precompute the mesh.
            grd = interpolate.griddata((lon.flatten(),lat.flatten()),grd.flatten(),(elons.flatten(),lats.flatten()),method='linear') # interpolate the data
            grd=np.reshape(grd,(self.lats.shape[0],self.elons.shape[0])) # reshape the data to the shape of a grid
            return grd
            
        def create_ICE(self):
            return spherical_ice()
        
        def create_SED(self):
            return spherical_sed()
        
        def create_TOPO(self):
            return spherical_topo()
        
        def create_LOVE(self,model_p,way):
            return LOVE(model_p,way)
        
        def create_SL(self):
            return spherical_sea_level(self)

        def disk(self,model_p,lat,lon,radius,high):
            disk=sphericalobject(np.array([0]))
            disk.disk=np.zeros((model_p.time_step_number+1,model_p.grid.lats.size,model_p.grid.elons.size))
            lon_g,lat_g=np.meshgrid(model_p.grid.elons,model_p.grid.lats)
            disk.disk[:,((lon-lon_g)**2+(lat-lat_g)**2)<radius]=high
            disk.disk[:2,:,:]=disk.disk[:2,:,:]*0
            return disk
