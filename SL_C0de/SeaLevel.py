from spharm import sphericalobject
import numpy as np



class spherical_sea_level(object):
    """
    A class used to represent the sea level spherical harmonic coefficient component.

    ...

    Attributes
    ----------
        del_L : np.array (maxdeg, maxdeg x2)
            Variation between the actual and past load on earth resulting from the redistribution of sediment ice thickness variation and geoid shpe modification. 
        
        del_L_prev : object (class sphericalobject)
            Previous variation between the actual and past load on earth resulting from the redistribution of sediment ice thickness variation and geoid shpe modification.
        
        TO : object (class sphericalobject)
            Variation of the topography between the ocean borders of two time step.
        
        sdelS : np.array (maxdeg, maxdeg x2)
            Variation of the sea level between two time step.
        
        sdelL : np.array (n_time_step, maxedg+1, maxdeg+1)
            Variation of the load at between the two time step
             
        delLa : np.array (n_time_step, maxedg+1, maxdeg+1)
            !!!! Je n'arrive pas à savoir ce que c'est !!!!
        
        delLa_prev : np.array (n_time_step, maxedg+1, maxdeg+1) !!!! A vérifier !!!!
            !!!! Pareil, ça dépend de la même chose, c'est le précedant mais à un time step de moins !!!! 
        
        sdelLa_prev : np.array (n_time_step,maxedg+1, maxdeg+1)
            !!!! Pareil , mais c'est la variation entre deux time step !!!!
        
        delSLcurl_fl : 
        
        delSLcurl : object (from class spherical object)
        
        RO : np.array (maxdeg, maxdeg)
            Ocean elevation on ocean (the geoïd cut by the ocean function)
        
        sdelS_new : np.array (maxdeg+1,maxdeg+1)
            The new geoïd created to calculated the convergence creterion. 
        
        delSL : np.array (maxdeg+1, maxdeg+1)
            The variation of the geoïd between actual and past. 
        
        saved : np.array (n, maxdeg, maxdeg) with n the number of time the save method is activated.
            The saved parameters as defined in the function. 
        
    Methods
    -------
        save(grid)
            Save a parameter each time the function is used and stack it with the others. 
                     

    """
    
    def __init__(self,grid):
        """
        Parameters
        ----------
        grid : object (object from class GRID)
        """
        self.TO=sphericalobject(np.zeros((int((grid.maxdeg+1)*(grid.maxdeg+2)/2),)),'coeff') #: coeff mais créer avec des grids
        self.sdelS=sphericalobject(np.zeros((int((grid.maxdeg+1)*(grid.maxdeg+2)/2),)),'coeff') #: coeff
        self.delS=sphericalobject(np.zeros((int((grid.maxdeg+1)*(grid.maxdeg+2)/2),)),'coeff')  #: coeff
        self.delSLcurl_fl=sphericalobject(np.zeros((int((grid.maxdeg+1)*(grid.maxdeg+2)/2),)),'coeff') #: coeff
        self.delSLcurl=sphericalobject(np.zeros((int((grid.maxdeg+1)*(grid.maxdeg+2)/2),)),'coeff')#: coeff
        self.RO=sphericalobject(np.zeros((int((grid.maxdeg+1)*(grid.maxdeg+2)/2),)),'coeff') 
        self.sdelS_new=0
        self.delSL=sphericalobject(np.zeros((grid.maxdeg,grid.maxdeg*2)))
        self.saved=np.array([])
        
    def save(self): # is it still usefull ? Check in code
        if self.saved.shape[0]==0:
            self.saved=np.array([self.delSL.grd.copy()])
        else :
            self.saved=np.concatenate((self.saved,np.array([self.delSL.grd.copy()])),axis=0)
        return self
    