from spharm import sphericalobject
import numpy as np


class spherical_ocean_function(sphericalobject):
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
        self.isgrd=False
        self.iscoeff=False
        self.saved=np.array([])
        self.rho=1000.
        
    def update_oc_0(self):
        '''
        Update the oc_0_grd and oc_0_coeff fields with the coeff and grd fields of the object. For the initialization of ocean function in the code
    
        Parameters : 
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            oc_0_grd (np.array): array of a gaussian grid of size  maxdeg, maxdeg x2.
            oc_0_coeff (np.array): array of the spherical coefficient of size maxdeg,maxdeg !!!! A vérifier !!!!.

        '''
        # The two types of data are used from oc_0 so we create the two of them.
        self.oc_0_grd=self.grd.copy()
        self.oc_0_coeff=self.coeff.copy()
        
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