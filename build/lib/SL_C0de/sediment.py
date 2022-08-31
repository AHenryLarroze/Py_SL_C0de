from spharm import sphericalobject
import numpy as np

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
    
    def load(self,SL_model):
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
        self.sed=np.zeros((len(SL_model.grid.time_step),SL_model.grid.nlats,SL_model.grid.nlons)) # create a matrix of 0 for the sediment thicness
        self.isgrd=True # set the value of grid to true. 
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