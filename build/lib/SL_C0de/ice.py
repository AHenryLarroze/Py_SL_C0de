from spharm import sphericalobject
import numpy as np
from scipy import io
import time
from multiprocessing import Pool,cpu_count,freeze_support # used for the multi processing calculation 
import itertools
import par_ice as par

class spherical_ice(sphericalobject):
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
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        self.rho=916.7
        
    def load(self,grid):
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
        ice_time=ice_time[0]
        grid.time_step=np.zeros(ice_time.shape)# using the ice time format to set the variable. It could be interesting to set the time vairable from an other parameters.
        self.ice=np.zeros((len(grid.time_step),grid.nlats,grid.nlons)) # precreate the ice matrix. 
        # same for the ice_corrected matrix. 
        #create an enlarged ice_lat and ice_long matrix
        ice_lat=np.zeros((data['ice_lat'].shape[0]+2,1)) 
        ice_long=np.zeros((data['ice_long'].shape[0]+2,1))
        ice_lat[0],ice_lat[-1],ice_lat[1:-1] =90,-90,data['ice_lat']
        ice_long[0],ice_long[-1],ice_long[1:-1] =-0.5,360.5,data['ice_long']
        # do the same thing with ice_extended
        # for each time step interpolate the ice grid on the Gaussian grid defined in GRID class
        for t_it in range(len(grid.time_step)):
            ice_nointerp=ice_in[t_it]
            ice_extended=np.concatenate((np.zeros((1,len(ice_long)-2)),ice_nointerp,ice_nointerp[0,-1]*np.ones((1,len(ice_long)-2))),axis=0)
            ice_extended_2=np.concatenate((np.expand_dims(ice_extended[:,-1],axis=1),ice_extended,np.expand_dims(ice_extended[:,0],axis=1)),axis=1)

            ice_interp = grid.interp_on(ice_extended_2,ice_long,ice_lat)

            self.ice[len(grid.time_step)-t_it-1,:,:] = ice_interp
            grid.time_step[t_it] = ice_time[len(ice_time)-t_it-1]
        
    def quick_load(self,grid,nb_workers):
        '''
        
        Function to load the ice data from the ice6g map. 
        This version was parallelized to accelerate the loading of the model.
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
        ice_time=ice_time[0]
        grid.time_step=np.zeros(ice_time.shape)# using the ice time format to set the variable. It could be interesting to set the time vairable from an other parameters.
        self.ice=np.zeros((len(grid.time_step),grid.nlats,grid.nlons)) # precreate the ice matrix. 
         # same for the ice_corrected matrix. 
        #create an enlarged ice_lat and ice_long matrix
        ice_lat=np.zeros((data['ice_lat'].shape[0]+2,1)) 
        ice_long=np.zeros((data['ice_long'].shape[0]+2,1))
        ice_lat[0],ice_lat[-1],ice_lat[1:-1] =90,-90,data['ice_lat']
        ice_long[0],ice_long[-1],ice_long[1:-1] =-0.5,360.5,data['ice_long']
        # do the same thing with ice_extended
        # for each time step interpolate the ice grid on the Gaussian grid defined in GRID class
        # In this version i created a function called par.py wich isa parallelised version of the for loop of the upper version. 
        # It need to define a pool object to precharge the core calculator.
        #workers = int(cpu_count()/2)# half of the core are used by the code.
        workers = nb_workers
        pool = Pool(workers)# set the bar (it don't work really well).
        results=pool.starmap(par.f_ice_quick_load,zip(itertools.repeat(ice_time),ice_in,range(len(ice_time)),itertools.repeat(ice_long),itertools.repeat(ice_lat),itertools.repeat(grid.lats),itertools.repeat(grid.elons)))
        grid.time_step=np.array([results[i][1] for i in range(len(results))])
        self.ice=np.array([results[len(results)-i-1][0] for i in range(len(results))])
        self.ice_corrected=self.ice.copy()
        
    def calc_del_ice(self,t_it) :# avec ça deli seras près définit et prev doit être calculer
        '''
        
        Calculate the variation of ice thickness between the biggiinig of the modelisation until the time step t_it. 
        
        Parameters :  
             t_it (int): The indices of the time step.
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            grd : A gaussian grid of the ice thickness variation at t_it, with the whape maxdeg, maxdeg x 2. 
        
        ''' 
        self.grd=self.ice_corrected[t_it,:,:]-self.ice_corrected[0,:,:] # calculate the difference of ice thickness between actual and the one at t_it.
        self.isgrd=True
        return self