import numpy as np
import math
from scipy import special as scsp
from scipy.special import lpmn as legendre_a
import scipy
import numpy.matlib as npmlib #used to add repmat 
#import pyshtools as pysh #pyshtools is still necessary, even if i don't use the expend function (still need the legendre  polynome function)
import sys
import logging
import pyshtools as pysh
import par_spharm as par
import itertools


def GaussQuad(N):
    '''
    Calculate the Gaussian grid parameters for the computation of the gaussian quadrature.
    
        Parameters : 
            N : int 
                The number of latitude in the grid. It's equal to the maximum degree
            
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
        
        Added fields : 
            w (np.array): The weight of the node in the gaussian grid
            x (np.array):The nodes coefficient in the gaussian grid
    '''
    beta=0.5/np.sqrt(1-(2*np.linspace(1,N-1,N-1))**(-2))
    T=np.diag(beta,1)+np.diag(beta,-1)
    D,V=np.linalg.eig(T)
    x=D
    i=np.argsort(x)
    x=x[i]
    w=2*V[0,i]**2
    return x,w

"""
def legendre(N):
    x,w=GaussQuad(N) # compute the parameter of the gaussian grid
    P_lm=[np.expand_dims(np.array([]),axis=0) for i in range(N+1)]# initialize the legendre polynôme !!!! calculer ces polynômes plutot que le faire chaque fois ira plus vite !!!!
    for l in range (N+1):
        # precalculate the normalization, like the Matlab code, we use the fully normalize version of the Legendre Polynôme
        m=np.linspace(0,l,l+1,dtype=int)
        #c=np.expand_dims(((-1)**(l))*np.sqrt((l+1/2)*(scsp.factorial(l-m))/scsp.factorial(l+m)),axis=1)
        for i in range(len(x)): 
            if i==0 : # if it's the first order of the degree we initialize the array
                #print('c',c.shape)
                #print('leg',np.expand_dims(pysh.legendre.legendre_lm(l,np.linspace(0,l,l+1,dtype=int),x[i],normalization='unnorm',csphase=-1),axis=1).shape)
                #P_lm[l]=c*np.expand_dims(pysh.legendre.legendre_lm(l,np.linspace(0,l,l+1,dtype=int),x[len(x)-i-1],normalization='4pi',csphase=-1),axis=1)
                P_lm[l]=np.expand_dims(pysh.legendre.legendre_lm(l,np.linspace(0,l,l+1,dtype=int),x[len(x)-i-1],normalization='4pi',csphase=-1),axis=1)
            else :
                #P_lm[l]=np.concatenate((P_lm[l],c*np.expand_dims(pysh.legendre.legendre_lm(l,np.linspace(0,l,l+1,dtype=int),x[len(x)-i-1],normalization='4pi',csphase=-1),axis=1)),axis=1) # use the pyshtools function to calclate the legendre polynôme
                P_lm[l]=np.concatenate((P_lm[l],np.expand_dims(pysh.legendre.legendre_lm(l,np.linspace(0,l,l+1,dtype=int),x[len(x)-i-1],normalization='4pi',csphase=-1),axis=1)),axis=1)
        #P_lm[l]=P_lm[l].transpose() # The output of concatenate is of the wrong shape so we transpose the P_lm
    return P_lm
"""

# def legendre(N):
#     x,w=GaussQuad(N) # compute the parameter of the gaussian grid
#     P_lm=[0 for i in range (N+1)]
#     for n in range (N+1):
#         if n==0:
#             P_lm[n]=np.expand_dims(scipy.special.lpmv(0,n,x)*math.sqrt((n+1/2)*(math.factorial(n-0)/math.factorial(n+0))),axis=0)
#         else :
#             for m in range (n+1):
#                 if m==89 and n==89 :
#                     print(m,n,math.sqrt((n+1/2)*(math.factorial(n-m)/math.factorial(n+m)))*(-1)**m,scipy.special.lpmv(m,n,x))
#                 if m==0:
#                     P_lm[n]=np.expand_dims(scipy.special.lpmv(m,n,x)*math.sqrt((n+1/2)*(math.factorial(n-m)/math.factorial(n+m)))*(-1)**m,axis=0)
#                 else :
#                     P_lm[n]=np.concatenate((P_lm[n] ,np.expand_dims(scipy.special.lpmv(m,n,x)*math.sqrt((n+1/2)*(math.factorial(n-m)/math.factorial(n+m)))*(-1)**m,axis=0)),axis=0)
#     return P_lm

def legendre(N,pool):
    x,w=GaussQuad(N) # compute the parameter of the gaussian grid
    P_lm=[0 for i in range (N+1)]
    P_lm_shtools=np.zeros((N,N+1,N+1))
    results=pool.starmap(par.f_legendre,zip(itertools.repeat(N),x))
    P_lm_shtools=np.array([results[i] for i in range(N)])
    for i in range(N+1):
        P_lm[i]=P_lm_shtools[:,i,:i+1].transpose()
    return P_lm
    


def get_coeffs(a_lm,n):
    if n == 0 :
        a_n = a_lm[0]
    else :
        gauss_sum = int(n*(n+1)/2)
        # position of current order in a_lm vector
        a_n = a_lm[gauss_sum:gauss_sum+n+1]
    return a_n

def calc_Y_lm(N,theta,phi):
    return pysh.expand.spharm(N,theta,phi,packed=True)

# def grd_to_coeff(grd,N,P_lm,w): # extract the shape of the grid
#     F_ym=np.fft.fft(grd,axis=1,norm='backward') # use the fourrier transform to sum over the exponential
#     ind_a=0
#     coeff=np.zeros((int((N+1)*(N+2)/2),))*1j # initialize the spherical harmonic coefficient vector
#     for l in range(0,N+1): # calculate for all degree the spherical harmonic coefficient
#         coeff[ind_a:ind_a+l+1]=np.sum(P_lm[l][:l+1,:].transpose()*npmlib.repmat(np.expand_dims(w,axis=1),1,l+1)*F_ym[:,:l+1],axis=0)
#         ind_a+=l+1
#     # apply a coerrection to the coefficient
#     coeff = coeff/(2*N)/math.sqrt(2)
#     return coeff

# def coeff_to_grd(coeff,N,P_lm):
#     F_ym=np.zeros((N,N*2))*1j
#     for l in range(N+1):
#         F_ym[:,:l+1] = F_ym[:,:l+1] + P_lm[l][:l+1,:].transpose()*npmlib.repmat(get_coeffs(coeff,l),N,1)
#     F_ym[:,N+1:] = np.conj(F_ym[:,N-1:0:-1]) 
#     F_yx = np.fft.ifft(F_ym,axis=1,norm='backward') 
#     grd = F_yx*(2*N)*math.sqrt(2)
#     grd=np.real(grd)
#     return grd


class sphericalobject(object):
    """
    A class difining any spherical object. This class is working with pyshtools.

    ...

    Attributes
    ----------
        grd : np.array (maxdeg,maxdeg x2)
           Value of the spherical object on a Gaussian Grid. 
        isgrd : Bool
           A boolean to define if a Gaussian grid have been defined for this object. 
        coeff : np.array (maxdeg,maxdeg) !!!! A vérifier !!!!
           Spherical harmonic coefficient array. 
        iscoeff : Bool
           A boolean to define if a spherical harmonic coefficient have been defined for this object.
        saved : np.array (n, maxdeg, maxdeg) !!!! A vérifier !!!! n is the number of time the save method is activated. 
           An array wich contain the spherical harmonic coefficient each time th save method is applied.
        prev : np.array (maxedg, maxdeg) !!!! A vérifier !!!!
           save the spherical coefficient using save_prev. 

    Methods
    -------
        grdtocoeff()
           Convert the Gaussian grid to spherical harmonic coefficient
        coefftogrd()
           Convert spherical harmonic coefficient to Gaussian grid
        multiply(flm2)
           Multiply the spherical coefficient of the object with the coefficient of flm2
        save_prev()
           Save the spherical harmonic coefficient to the attribute prev
        save()
           Save the spherical harmonic coefficient to the attribute saved and stack with previous saved coefficient
        modify(gc,t='grd')
           change the grid or coeff attribute with gc depending on t value    
    """
    
    def __init__(self,grd,t='grd'):
        """
        Parameters
        ----------
        grd : np.array
            array input a gaussian grid or spherical harmonic coefficient
        t : str
            value defining the type of the grd input, can be 'grd' or 'coeff' (default 'grd')
        """
        
        if t=='grd' : # initialize the grid if the entry is a grid
            self.grd=grd.copy()
            self.isgrd=True
            self.iscoeff=False
        elif t=='coeff': # initialize the coefficient if the entry is a coefficient
            self.coeff=grd.copy()
            self.isgrd=False
            self.iscoeff=True
        self.saved=np.array([]) # initialize the save of the spherical harmonic object
        
    def grdtocoeff(self,model_p):
        '''
    Convert a Gaussian grid into spherical harmonic coefficient array using a numerical method to create spherical harmonic coefficient
    self.coeff is updated usig these output.
    self.iscoeff defining if a coefficient have been created for this object. 
    If there is no grid created (self.isgrid == 0) then it returns an error.
    A modifier pour pouvoir modifier les entrées de la fonction.
    
        Parameters : 
            
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
        
        Added fields : 

        '''
        if self.isgrd: # if there is a grid convert the gaussian grid to spherical coefficient.
            nlat,nlon=self.grd.shape # extract the shape of the grid
            F_ym=np.fft.fft(self.grd,axis=1) # use the fourrier transform to sum over the exponential
            ind_a=0
            self.coeff=np.zeros((int((nlat+1)*(nlat+2)/2),))*1j # initialize the spherical harmonic coefficient vector
            # for l in range(0,nlat+1): # calculate for all degree the spherical harmonic coefficient
            #     self.coeff[ind_a:ind_a+l+1]=np.sum(model_p.P_lm[l][:l+1,:].transpose()*npmlib.repmat(np.expand_dims(model_p.w,axis=0).transpose(),1,l+1)*F_ym[:,:l+1],axis=0)
            #     ind_a+=l+1
            results=model_p.pool.starmap(par.f_grdtocoeff,zip([l for l in range(nlat+1)],model_p.P_lm,itertools.repeat(model_p.w),itertools.repeat(F_ym)))
            self.coeff=np.concatenate(results)
            # apply a coerrection to the coefficient
            self.coeff = self.coeff/(2*nlat)/math.sqrt(2)
            self.iscoeff=True # Set the spherical harmonic coefficient test to true.
        else : # else there is no grid and then return an error.
            logging.error('error: ', "No map created for this spherical object. Check if you have created the object with map or haven't run the coefftogrd() method")
            sys.exit(1)
        return self
                                                                                  
    def coefftogrd(self,model_p):
        '''
    Convert spherical harmonic coefficient into a gird array using shtools.
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
        if self.iscoeff :# if there is spherical coefficient convert the spherical coefficient to a gaussian grid.
            m=self.coeff.shape[0]
            N=int((math.sqrt(1+8*m)-1)/2-1) 
            #P_lm,w=legendre(N)
            F_ym=np.zeros((N,N*2))*1j
            # for l in range(N+1):
            #     F_ym[:,:l+1] = F_ym[:,:l+1] + model_p.P_lm[l][:l+1,:].transpose()*npmlib.repmat(get_coeffs(self.coeff,l),N,1)
            results=model_p.pool.starmap(par.f_coefftogrd,zip([l for l in range(N+1)],model_p.P_lm,itertools.repeat(self.coeff),itertools.repeat(N)))
            for l in range(N+1):
                F_ym[:,:l+1]=F_ym[:,:l+1]+results[l]
            F_ym[:,N:] = np.conj(F_ym[:,N:0:-1])
            F_yx = np.fft.ifft(F_ym,axis=1) 
            self.grd = F_yx*(2*N)*math.sqrt(2)
            self.grd=np.real(self.grd)
            self.isgrd=True 
        else : # else there is no coefficient and then return an error.
            logging.error('error: ', "No coeff created for this spherical object. Check if you have created the object with coeff or haven't run the grdtocoeff() method")
            sys.exit(1)
        return self

    def coefftogrd_hd(self,model_p):
        '''
        Convert spherical harmonic coefficient into a gird array using shtools.
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
        if self.iscoeff :# if there is spherical coefficient convert the spherical coefficient to a gaussian grid.
            #P_lm,w=legendre(N)
            self.coeff=np.concatenate((self.coeff,np.zeros((int((model_p.maxdeg_hd+1)*(model_p.maxdeg_hd+2)/2)-len(self.coeff),))))
            F_ym=np.zeros((model_p.maxdeg_hd,model_p.maxdeg_hd*2))*1j
            for l in range(model_p.maxdeg_hd+1):
                F_ym[:,:l+1] = F_ym[:,:l+1] + model_p.P_lm_res[l][:l+1,:].transpose()*npmlib.repmat(get_coeffs(self.coeff,l),model_p.maxdeg_hd,1)
            F_ym[:,model_p.maxdeg_hd:] = np.conj(F_ym[:,model_p.maxdeg_hd:0:-1]) 
            F_yx = np.fft.ifft(F_ym,axis=1) 
            grd = F_yx*(2*model_p.maxdeg_hd)*math.sqrt(2)
            grd=np.real(grd)
        else : # else there is no coefficient and then return an error.
            logging.error('error: ', "No coeff created for this spherical object. Check if you have created the object with coeff or haven't run the grdtocoeff() method")
            sys.exit(1)
        return grd
    
    def multiply(self,flm2):
        '''
    This function is no longer used !!!!    
    
    Multiply the spherical coefficient of a gird with the others. It could be done using the pysh.SHcoeffs.Multiply.
    I choose to do it manually but i could test the efficiency of the pysh tool function. 
    
        Parameters : 
            flm2 (np.array): the harmonic coefficient matrix of size maxdeg, maxdeg !!!! A vérifier !!!!!
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
            flm[0,:,:]+1j*flm[1,:,:] (np.array): the harmonic coefficient matrix resulting of the multiplication of the two grids. 
            size maxdeg, maxdeg !!!! A vérifier !!!!
            
        Added fields : 

        '''
        # convert two spherical coefficient into gaussian grid.
        flm1=np.array([np.real(self.coeff),np.imag(self.coeff)])# for the object
        Clm1=pysh.shio.SHctor(flm1)
        flm2=np.array([np.real(flm2),np.imag(flm2)])# for the matrix in entry
        Clm2=pysh.shio.SHctor(flm2)
        grd1=pysh.SHCoeffs.from_array(Clm1).expand(grid='GLQ', lat=None, colat=None, lon=None, degrees=True, zeros=None, lmax=None, lmax_calc=None, extend=True, backend=None, nthreads=0).to_array()
        grd2=pysh.SHCoeffs.from_array(Clm2).expand(grid='GLQ', lat=None, colat=None, lon=None, degrees=True, zeros=None, lmax=None, lmax_calc=None, extend=True, backend=None, nthreads=0).to_array()                                         
        grd=grd1*grd2 # multiply the two grids.
        # reconvert the resulting gaussian grid into spherical harmonic coefficient using pyshtools.
        Clm=pysh.SHGrid.from_array(grd,'GLQ').expand(normalization='4pi', csphase=1, lmax_calc=None, backend=None, nthreads=0).to_array()
        flm= pysh.shio.SHrtoc(Clm)
        return flm[0,:,:]+1j*flm[1,:,:]
    
    def save_prev(self):
        '''
    Create a new field for the object to save the spherical coefficient at the moment of the applied function. 
    This function make a clean copy of the array to avoid modification. 
    
        Parameters : 
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
            
        Added fields : 
            prev (np.array): a copy of the self.coeff spherical coeffcient matrix of size maxdeg, maxdeg !!!! A vérifier !!!!

        '''
        if self.iscoeff : # if there is spherical coefficient for this object copy these coefficient to self.prev
            self.prev=self.coeff.copy()
        else : # else retrun an error
            logging.error('error: ', "No coeff created for this spherical object. Check if you have created the object with coeff or haven't run the grdtocoeff() method")
            sys.exit(1)
        return self
    
    def save(self):  
        '''
    Create a 3D matrix containing the spherical coefficient of each time this function is activated or add a copy of the
    spherical coefficient to the already existing matrix.
    
        Parameters : 
            flm2 (np.array): the harmonic coefficient matrix of size maxdeg, maxdeg !!!! A vérifier !!!!!
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
            flm[0,:,:]+1j*flm[1,:,:] (np.array): the harmonic coefficient matrix resulting of the multiplication of the two grids. 
            size maxdeg, maxdeg !!!! A vérifier !!!!
            
        Added fields : 

        '''
        if self.iscoeff : # if there is spherical coeffcient created
            if self.saved.shape[0]==0: # if the matrix don't exist create one
                self.saved=np.array([self.coeff])
            else : # else add to the matrix the spherical coefficient of the object
                self.saved=np.concatenate((self.saved,np.array([self.coeff.copy()])),axis=0)
        else : # else retrun an error
            logging.error('error: ', "No coeff created for this spherical object. Check if you have created the object with coeff or haven't run the grdtocoeff() method")
            sys.exit(1)
        return self
    
    def modify(self,gc,t='grd'):
        '''
    Modify the coeff or grid field of the object depending on the t variable. 
    Default value for t is grd. 
    
        Parameters : 
            gc (np.array): spherical coefficient matrix of size maxdeg, maxdeg !!!! A vérifier !!!! or
                           Gaussian grid of the data of size maxdeg, maxdegx2. 
            t (str): string value to set the type of the gc entry. 'grd' for grid, 'coeff' for spherical coefficient.
             
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns :
            
        Added fields : 
            grd (np.array): array of size maxdeg, maxdeg. !!!! A vérifier !!!!
            coeff (np.array): array of spherical coefficient of size maxdeg,maxdeg. !!!! A vérifier !!!!
            iscoeff (bool): boolean value, true if there is a coefficient false otherwise
            isgrd (bool): boolean value, true if there is a grid false otherwise

        '''
        if t=='grd':# if t is 'grd' modify self.grd with gc and update the iscoeff and isgrd field to false and true respectively.
            self.grd=gc.copy()
            self.iscoeff=False
            self.isgrd=True
        elif t=='coeff': # if t is 'coeff' modify self.coeff with gc and update the iscoeff and isgrd field to true and false respectively.
            self.coeff=gc.copy()
            self.isgrd=False
            self.iscoeff=True
        else : # else there is no recgnize value of t and then return an error.
            logging.error('error: ', "Wrong entry for t parameter. Verify the entry is 'grd' or 'coeff'.")
            sys.exit(1)
        return self