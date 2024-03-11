import numpy as np
import math
#import pyshtools as pysh #pyshtools is still necessary, even if i don't use the expend function (still need the legendre  polynome function)
import sys
import logging
import pyshtools as pysh
import pyshtools.expand as expand
    
def get_coeffs(a_lm,n):
    '''
    The _`get_coeffs` function get the spherical harmonics coefficients of the nth order from a linearly Spharm array.

    Attribute : 
    -----------
        a_lm : np.array([maxdeg*(maxdeg+1)/2,])
            An array containing the spherical coefficient in a linear form
        n : int
            The order of the spherical harmonics coefficient you are trying to retrive

    Returns :
    ---------  
        a_n : np.array([n,])
            The spherical harmonic coefficient of the order n.
    '''
    if n == 0 :
        a_n = a_lm[0]
    else :
        gauss_sum = int(n*(n+1)/2)
        # position of current order in a_lm vector
        a_n = a_lm[gauss_sum:gauss_sum+n+1]
    return a_n

class sphericalobject(object):
    """
    The _`sphericalobject` class include any spherical object. This class is working with pyshtools (`pyshtools <https://shtools.github.io/SHTOOLS/>`_).

    Attributes
    ----------
        grd : np.array[(maxdeg,maxdeg x2)]
           Value of the spherical object on a Gaussian Grid. 
        isgrd : Bool
           A boolean to define if a Gaussian grid have been defined for this object. 
        coeff : np.array[(maxdeg,maxdeg)]
           Spherical harmonic coefficient array. 
        iscoeff : Bool
           A boolean to define if a spherical harmonic coefficient have been defined for this object.
        saved : np.array[(n, maxdeg, maxdeg)]
           An array wich contain the spherical harmonic coefficient each time th save method is applied.
        prev : np.array[(maxedg, maxdeg)] 
           save the spherical coefficient using save_prev. 

    Methods
    -------
        `grdtocoeff`_ : 
           Convert the Gaussian grid to spherical harmonic coefficient
        `coefftogrd`_ : 
           Convert spherical harmonic coefficient to Gaussian grid
        `coefftogrdhd`_ :
           Convert spherical harmonic coefficient to a Gaussian grid with a higher resolution then maxdeg
        `save_prev`_ :
           Save the spherical harmonic coefficient to the attribute prev
        
           
    """
    
    def __init__(self,grd=None,coeff=None,maxdeg=None):
        
        if coeff is None : # initialize the grid if the entry is a grid
            self.grd=grd.copy()
            self.maxdeg=grd.shape[0]
        elif grd is None: # initialize the coefficient if the entry is a coefficient
            self.coeff=coeff.copy()
            self.maxdeg=int(abs((-1+math.sqrt(1+8*len(coeff)))/2))
        else :
            self.maxdeg=maxdeg
        
    def grdtocoeff(self):
        '''
        The _`grdtocoeff` method convert a Gaussian grid into spherical harmonic coefficient array using a numerical method to create spherical harmonic coefficient
        self.coeff is updated usig these output.
        self.iscoeff defining if a coefficient have been created for this object. 
        If there is no grid created (self.isgrid == 0) then it returns an error.
    
        Attribute :
        ----------- 
            None 
        
        Returns :
        ---------
            None        
        '''
        zero , w = expand.SHGLQ(self.maxdeg)
        self.coeff=expand.SHExpandGLQ(self.grd,w,zero)
        self.coeff=pysh.shio.SHCilmToCindex(self.coeff)
        self.coeff=self.coeff[0]+self.coeff[1]*1j
        self.iscoeff=True
        return self
                                                                                  
    def coefftogrd(self):
        '''
        The _`coefftogrd` method convert spherical harmonic coefficient into a gird array using shtools.
        The output of pysh.SHCoeff are converted to real.
        self.grd is updated usig these output.
        self.isgrd defining if a grid have been created for this object. 
        If there is no grid created (self.iscoeff == 0) then it returns an error.
        A modifier pour pouvoir modifier les entrées de la fonction.
    
        Parameters : 
            
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
        
        Added fields : 

        '''       
        zero , w = expand.SHGLQ(self.maxdeg)
        coeff=np.stack((self.coeff.real,self.coeff.imag))
        coeff=pysh.shio.SHCindexToCilm(coeff)
        self.grd=expand.MakeGridGLQ(coeff,zero,extend=1)
        self.isgrd=True

        return self
    
    def coefftogrdhd(self,max_calc_deg):
        '''
        The _`coefftogrdhd` convert spherical harmonic coefficient into a gird array using shtools.
        The output of pysh.SHCoeff are converted to real.
        self.grd is updated usig these output.
        self.isgrd defining if a grid have been created for this object. 
        If there is no grid created (self.iscoeff == 0) then it returns an error.
        A modifier pour pouvoir modifier les entrées de la fonction.
    
        Attribute :
        -----------
            max_calc_deg : int
                The maximum spherical harmonic degree to calculate the grid. This function can be used to have a better rendering in output.
        
        Returns :
        ---------
            None 
        '''
        
        zero , w = expand.SHGLQ(max_calc_deg-1)
        x_GL = np.arccos(zero[::-1])*180/math.pi - 90
        lon_GL = np.linspace(0,360,2*(max_calc_deg)+1)
        lon_GL = lon_GL[:-1]
        lat_hd=x_GL
        lon_hd=lon_GL
        self.colats = 90 - self.lats
        coeff=np.stack((self.coeff.real,self.coeff.imag))
        coeff=pysh.shio.SHCindexToCilm(coeff)
        return expand.MakeGridGLQ(coeff,zero,lmax=max_calc_deg-1,extend=1),lon_hd,lat_hd[::-1]
    
    # def multiply(self,coeff2):
    #     '''
    #     This function is no longer used !!!!    
        
    #     Multiply the spherical coefficient of a gird with the others. It could be done using the pysh.SHcoeffs.Multiply.
    #     I choose to do it manually but i could test the efficiency of the pysh tool function. 
    
    #     Parameters : 
    #         flm2 (np.array): the harmonic coefficient matrix of size maxdeg, maxdeg !!!! A vérifier !!!!!
             
    #     See the documentation of the cited class object for more information on different parameters used in the function.
        
    #     Returns : 
    #         flm[0,:,:]+1j*flm[1,:,:] (np.array): the harmonic coefficient matrix resulting of the multiplication of the two grids. 
    #         size maxdeg, maxdeg !!!! A vérifier !!!!
            
    #     Added fields : 

    #     '''
    #     # convert two spherical coefficient into gaussian grid.
    #     coeff1=np.stack((self.coeff.real,self.coeff.imag))
    #     coeff1=pysh.shio.SHCindexToCilm(coeff1)
    #     coeff2=np.stack((coeff2.real,coeff2.imag))
    #     coeff2=pysh.shio.SHCindexToCilm(coeff2)
    #     coeff=pysh.expand.SHMultiply(coeff1,coeff2)
    #     coeff=pysh.shio.SHCilmToCindex(coeff)
    #     coeff=coeff[0]+coeff[1]*1j
    #     return coeff[:int((self.maxdeg*(self.maxdeg+1))/2)]
    
    # def calculate_value_at_point(self,phi,theta):
    #     Y_lm=(pysh.expand.spharm(self.maxdeg,theta,phi,packed=True)[0]+pysh.expand.spharm(self.maxdeg,theta,phi,packed=True)[1]*1j)
    #     return np.dot(self.coeff,Y_lm)
    
    def save_prev(self):
        '''
        The _`save_prev` create a new field for the object to save the spherical coefficient at the moment of the applied function. 
        This function make a clean copy of the array to avoid modification. 
    
        Attribute :
        -----------
            None 
        
        Returns :
        --------- 
            None
        '''
        if not(self.coeff is None) : # if there is spherical coefficient for this object copy these coefficient to self.prev
            self.prev=self.coeff.copy()
        else : # else retrun an error
            logging.error('error: ', "No coeff created for this spherical object. Check if you have created the object with coeff or haven't run the grdtocoeff() method")
            sys.exit(1)
        return self
