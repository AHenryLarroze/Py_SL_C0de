from .spharm import sphericalobject
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

    def calc_sdelS(self,model_p,k,topo_it=0):
        if k == 0 and topo_it==0:
            #print('oc',oc.coeff.shape,'TO',SL.TO.coeff.shape,'TO_prev',SL.TO.prev.shape,'ICE',ice.coeff.shape,'SED',sed.coeff.shape)
            model_p.SL.sdelS.modify(model_p.oc.prev/model_p.oc.prev[0]*(-model_p.ice.rho/model_p.oc.rho*model_p.ice.sdeli_00 + model_p.SL.TO.coeff[0]-model_p.SL.TO.prev[0])-model_p.SL.TO.coeff-model_p.SL.TO.prev,'coeff')
        model_p.SL.delS.modify(model_p.SL.delS.prev + model_p.SL.sdelS.coeff,'coeff')

    def calc_sdelS_new(self,model_p):
        self.delSLcurl_fl.modify(model_p.love.E.coeff * model_p.love.T.coeff *model_p.load.delL.coeff+model_p.love.T.coeff*model_p.load.V_lm.coeff+1/model_p.g*model_p.love.E_T.coeff*model_p.load.delLa.coeff+1/model_p.g*model_p.load.V_lm_T.coeff,'coeff')
        self.delSLcurl.modify(self.delSLcurl_fl.coeff - model_p.ice.coeff+ -model_p.sed.coeff,'coeff').coefftogrd(model_p)
        self.RO.modify(self.delSLcurl.grd*model_p.oc.grd).grdtocoeff(model_p)
        self.delPhi_g = np.real(1/model_p.oc.coeff[0] * (- model_p.ice.rho/model_p.oc.rho*model_p.ice.coeff[0] - self.RO.coeff[0] + self.TO.coeff[0]))
        self.sdelS_new=self.RO.coeff + self.delPhi_g*model_p.oc.coeff -  self.TO.coeff - self.delS.prev
    
    def save(self): # is it still usefull ? Check in code
        if self.saved.shape[0]==0:
            self.saved=np.array([self.delSL.grd.copy()])
        else :
            self.saved=np.concatenate((self.saved,np.array([self.delSL.grd.copy()])),axis=0)
        return self
    