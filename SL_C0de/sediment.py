from spharm import sphericalobject
import numpy as np
from scipy import io
import par_sed as par
import itertools

class spherical_sed(sphericalobject):
    """
    A class used to represent the sediment grid

    ...

    Attributes
    ----------
        rho : float
           Density value of the sediment as a constant
        sed : np.array (maxdeg, maxdeg x2)
           Gaussian grid with the sediment thickness
    Methods
    -------
        load(grid)
            load the sediment grid, it actually just create a grid of 0.
        calc_del_sed(t_it)
            calculate the sediment thickness deposited between the time step t_it, t_it+1
    """
    def __init__(self) : 
        """
        Parameters
        ----------
        """
        # preset the grid and coefficient to none. The sediment density is set to 2300 
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        self.rho = 2300.
    
    def load(self,model_p,way):
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
        sed_lon=data['lon']
        sed_lat=data['lat']
        sed_mat_inp=data['Sed_map_fin']
        # for i in range(len(model_p.grid.time_step)):
        #     sed_interp=model_p.grid.interp_on(sed_mat_inp[i,:,:],sed_lon,sed_lat)
        #     self.sed[i,:,:]=sed_interp
        results=model_p.pool.starmap(par.f_sed_quick_load,zip(sed_mat_inp,itertools.repeat(sed_lon),itertools.repeat(sed_lat),itertools.repeat(model_p.grid.lats),itertools.repeat(model_p.grid.elons)))
        self.sed=np.array(results)
        return self
        
    def calc_del_sed(self,t_it) :
        '''
        Calculate the sediment thickness variation between the time step t_it and the bigining of the modelisation. 
    
        Parameters : 
            grid (object): The output of the class GRID. 
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            grd (np.array): the grd matrix of shape maxdeg, maxdeg x2 of sediment thickness.

        '''
        self.grd=self.sed[0,:,:]-self.sed[t_it,:,:] # calculate the difference between actual sediment thickness and the one at t_it
        self.isgrd=True # because a grid was created, update the isgrd parameter.
        return self
    
    def zeros(self,model_p):
        self.sed=np.zeros((model_p.time_step_number+1,model_p.grid.lats.size,model_p.grid.elons.size))