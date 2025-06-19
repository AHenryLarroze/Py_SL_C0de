import numpy as np
import math
import sys
import logging
import pyshtools as pysh
import pyshtools.expand as expand
import random
from scipy.special import lpmv, factorial
from scipy.fft import ifft
    
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
    if n == 0:
        a_n = a_lm[0]
    else:
        gauss_sum = n * (n + 1) // 2
        # Indices for the current order in the a_lm vector
        vec = slice(gauss_sum, gauss_sum + n + 1)
        a_n = a_lm[vec]

    return a_n

class sphericalobject(object):
    """
    s"""
    
    def __init__(self,grd=None,coeff=None,maxdeg=None,P_lm=None):
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
        # print(coeff,grd,maxdeg)
        if not(maxdeg is None) :
            self.maxdeg=maxdeg
        if not(coeff is None) : # initialize the grid if the entry is a grid
            self.coeff=coeff.copy()
            self.maxdeg=int(abs((-1+math.sqrt(1+8*len(coeff)))/2)-1)
            self.grd=np.array([0])
        if not(grd is None): # initialize the coefficient if the entry is a coefficient
            self.grd=grd.copy()
            self.coeff=np.array([0+0j])
            self.maxdeg=grd.shape[0]

        x,w=pysh.expand.SHGLQ(self.maxdeg-1)
        # x=x[::-1]
        x_GL = np.arccos(x)*180/math.pi - 90
        lon_GL = np.linspace(0,360,2*self.maxdeg+1)
        lon_GL = lon_GL[:-1]
        self.lats=x_GL.copy()
        self.elons=lon_GL.copy()
        self.colats = 90 - self.lats
        #correction=np.repeat((-1)**(np.arange(0,self.maxdeg+1)[:,np.newaxis]),self.maxdeg+1,axis=1)
        if P_lm is None :
            self.P_lm=self.calc_Plm(self.maxdeg)
        else :
            self.P_lm = P_lm.copy()
        
    def rng(self, num_gaussians):
        """
        Generate a 2D array of values from a randomly defined polynomial function.

        Parameters:
        - x_range: tuple, (x_min, x_max), range of x values.
        - y_range: tuple, (y_min, y_max), range of y values.
        - degree: int, maximum degree of the polynomial terms.
        - array_shape: tuple, (rows, cols), shape of the output array.

        Returns:
        - poly_array: 2D numpy array with values of the polynomial.
        - coefficients: list of random coefficients for each term in the polynomial.
        """
        # Generate x and y grids
        x = self.elons.copy()
        y = self.lats.copy()
        X, Y = np.meshgrid(x, y)

        # Generate random coefficients for each term in the polynomial
        grid=np.zeros(X.shape)
        # Parameters
        num_gaussians = 20      # Number of Gaussian curves
        amplitude_range = (1, 10)  # Range of Gaussian amplitudes
        width_range = (5, 20)      # Range of Gaussian widths (standard deviation)
        seed = 42                 # Seed for reproducibility

        # Initialize grid
        np.random.seed(seed)

        for _ in range(num_gaussians):
            # Random parameters for the Gaussian
            amplitude = np.random.uniform(*amplitude_range)
            width = np.random.uniform(*width_range)
            center_x = np.random.randint(0, len(self.elons))
            center_y = np.random.randint(0, len(self.lats))

            # Gaussian equation
            gaussian = amplitude * np.exp(-(((X - center_x)**2 + (Y - center_y)**2) / (2 * width**2)))

            # Add the Gaussian to the grid
            grid += gaussian

        self.grd=grid.copy()
        
    def grdtocoeff(self,type='shtools'):
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
        
        # zero=zero[::-1]
        if type == 'shtools':
            zero , w = expand.SHGLQ(self.maxdeg)
            # self.grd=self.grd[:,::-1]
            # self.coeff=expand.SHExpandGLQC(self.grd,w,zero)
            # print(np.concatenate((self.grd,(self.grd[:,-1])[:,np.newaxis]),axis=1).shape)
            # self.coeff=expand.SHExpandGLQC(np.concatenate((self.grd,(self.grd[:,-1])[:,np.newaxis]),axis=1),w,zero,csphase=1,lmax_calc=self.maxdeg)
            # self.coeff=expand.SHExpandGLQC(np.concatenate((np.concatenate((self.grd,(self.grd[:,0])[:,np.newaxis]),axis=1),(np.concatenate((self.grd,(self.grd[:,0])[:,np.newaxis]),axis=1)[0,:])[np.newaxis,:]),axis=0),w,zero,csphase=1)
            self.coeff=expand.SHExpandGLQC(np.concatenate((np.concatenate(((self.grd[:,-1])[:,np.newaxis],self.grd),axis=1),(np.concatenate(((self.grd[:,-1])[:,np.newaxis],self.grd),axis=1)[-1,:])[np.newaxis,:]),axis=0),w,zero,csphase=1)
            # self.grd=self.grd[:,::-1]
            # self.coeff[1,:,0]=self.coeff[0,:,0]

            coeff_real=pysh.shio.SHCilmToCindex(np.real(self.coeff))
            coeff_imag=pysh.shio.SHCilmToCindex(np.imag(self.coeff))
            
            # self.coeff=coeff_real[0]+coeff_real[1]-(coeff_imag[0]+coeff_imag[1])*1j
            self.coeff=coeff_real[0]-coeff_imag[0]*1j
            # self.coeff=np.concatenate((self.coeff,np.zeros((self.maxdeg+1,))+0j),axis=0)

            # self.coeff=pysh.shio.SHCilmToCindex(self.coeff)
        else :
            zero , w = expand.SHGLQ(self.maxdeg-1)
            # Use FFT to sum over exponential
            F_ym = np.fft.fft(self.grd, axis=1)
            
            # Initialize output array
            num_coefficients = (self.maxdeg + 1) * (self.maxdeg + 2) // 2
            a_lm = np.zeros(num_coefficients, dtype=np.complex128)
            
            # Perform Gauss-Legendre quadrature
            ind_a = 0
            for l in range(self.maxdeg + 1):
                legendre_part = self.P_lm[l,:l+1,:].T  # Select relevant part of Legendre polynomials
                weight_part = np.tile(w[:, np.newaxis], (1, l+1))  # Broadcast weights
                fft_part = F_ym[:, :l+1]  # Corresponding FFT coefficients
                
                # Sum over quadrature points
                a_lm[ind_a:ind_a + l + 1] = np.sum(legendre_part * weight_part * fft_part, axis=0)
                ind_a += l + 1
            
            # Apply additional factors
            self.coeff = a_lm / (2 * self.maxdeg) / 2#np.sqrt(2)
        

        
        self.iscoeff=True
        return self
                                                                                  
    def coefftogrd(self,type='shtools'):
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
        
        # zero=zero[::-1]      
        if type == 'shtools':
            zero , w = expand.SHGLQ(self.maxdeg-1)
            coeff_imag=-pysh.shio.SHCindexToCilm(np.stack((np.imag(self.coeff[:int(self.maxdeg*(self.maxdeg+1)/2)]),np.imag(self.coeff[:int(self.maxdeg*(self.maxdeg+1)/2)]))))
            coeff_real=pysh.shio.SHCindexToCilm(np.stack((np.real(self.coeff[:int(self.maxdeg*(self.maxdeg+1)/2)]),np.real(self.coeff[:int(self.maxdeg*(self.maxdeg+1)/2)]))))

            coeff=coeff_real.copy()+coeff_imag.copy()*1j
            correction=np.repeat((-1)**(np.arange(0,self.maxdeg)[:,np.newaxis]),self.maxdeg,axis=1)
            correction_2=-correction.copy()
            # coeff_imag=pysh.shio.SHCindexToCilm(np.stack((np.imag(self.coeff),np.imag(self.coeff))))
            # coeff_real=pysh.shio.SHCindexToCilm(np.stack((np.real(self.coeff),np.real(self.coeff))))
            # coeff_imag[0]=coeff_imag[0]*correction_2.T
            coeff_real[0]=coeff_real[0]*correction.T
            coeff_imag[0]=coeff_imag[0]*correction_2.T
            # coeff_real[0]=coeff_real[0]*correction.T
            coeff_neg=coeff_real.copy()+coeff_imag.copy()*1j

            # coeff[0,:,:]=coeff_neg[0,:,:]
            coeff[1,:,:]=coeff_neg[0,:,:]
            coeff[1,:,0]=0+0j
            # print(coeff[:,2,:])
            # print(coeff[:,:3,:3])
            lon,lat=np.meshgrid(self.elons,self.lats[::-1])
            grd_points=pysh.expand.MakeGridPointC(coeff,lat.flatten(),lon.flatten())
            self.grd=np.reshape(np.real(grd_points),(len(self.lats),len(self.elons)))
            # self.grd=np.real(expand.MakeGridGLQC(coeff,zero,extend=True,csphase=1))
            # self.grd=self.grd[:-1,1:]
            # self.grd=np.concatenate((self.grd,(self.grd[:,-1])[:,np.newaxis]),axis=1)
            self.isgrd=True
            # self.grd=self.grd[:,::-1]
        else : 
            zero , w = expand.SHGLQ(self.maxdeg-1)
            F_ym = np.zeros((self.maxdeg, 2 * self.maxdeg), dtype=complex)
    
            for l in range(0, self.maxdeg + 1):
                # Placeholder for P_lm[l], which should be of shape (l+1, N)
                # P_lm_l = self.P_lm[l]  # You need to define P_lm appropriately
                # Transpose and select first l+1 rows
                P_lm_l_T = self.P_lm[l,:l+1, :].T  # Shape: (N, l+1)
                # Get coefficients for degree l
                coeffs = get_coeffs(self.coeff, l)  # Define get_coeffs based on your data structure
                # Expand coefficients to shape (N, l+1)
                coeffs_expanded = np.tile(coeffs, (self.maxdeg, 1))
                # Multiply element-wise and accumulate
                F_ym[:, :l+1] += P_lm_l_T * coeffs_expanded

            # Use inverse FFT to sum over exponentials
            # F_ym[:, self.maxdeg:] = np.conj(F_ym[:, self.maxdeg-2::-1])
            F_ym[:, self.maxdeg+1:] = np.conj(F_ym[:, self.maxdeg-1:0:-1])
            F_yx = np.fft.ifft(F_ym, axis=1)

            # Apply additional factors
            F_yx *= (2 * self.maxdeg)#np.sqrt(2)

            self.grd=np.real(F_yx.copy())
        
        return self
    
    def coefftogrdhd(self,max_calc_deg,P_lm_hd=None,type='Jacky'):
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
        x_GL = np.arccos(zero)*180/math.pi - 90
        lon_GL = np.linspace(0,360,2*(max_calc_deg))
        lat_hd=x_GL.copy()
        lon_hd=lon_GL.copy()
        if type == 'shtools':
        # zero=zero[::-1]
            self.colats = 90 - self.lats
            coeff_imag=pysh.shio.SHCindexToCilm(np.stack((np.imag(self.coeff),np.imag(self.coeff))))
            coeff_real=pysh.shio.SHCindexToCilm(np.stack((np.real(self.coeff),np.real(self.coeff))))
            coeff=coeff_real.copy()+coeff_imag.copy()*1j
            correction=np.repeat((-1)**(np.arange(0,self.maxdeg)[:,np.newaxis]),self.maxdeg,axis=1)
            correction_2=-correction.copy()
            coeff_imag=pysh.shio.SHCindexToCilm(np.stack((np.imag(self.coeff),np.imag(self.coeff))))
            coeff_real=pysh.shio.SHCindexToCilm(np.stack((np.real(self.coeff),np.real(self.coeff))))
            coeff_imag[0]=coeff_imag[0]*correction_2.T
            coeff_real[0]=coeff_real[0]*correction.T
            coeff_neg=coeff_real.copy()+coeff_imag.copy()*1j

            coeff[1,:,:]=coeff_neg[0,:,:]
            coeff[1,:,0]=0
            grd_hd=np.real(expand.MakeGridGLQC(coeff,zero,lmax=max_calc_deg-1,extend=1,csphase=-1))
        else : 

            if P_lm_hd is None :
                P_lm_hd=self.calc_Plm(max_calc_deg)
            else :
                P_lm_hd = P_lm_hd.copy()

            if max_calc_deg>self.maxdeg :
                coeff_hd=np.concatenate((self.coeff,np.zeros(int((max_calc_deg+1)*(max_calc_deg+2)/2-self.coeff.shape[0]),dtype=np.complex128)),axis=0)
            elif max_calc_deg<self.maxdeg : 
                coeff_hd=self.coeff[:int((max_calc_deg+1)*(max_calc_deg+2)/2)]
            F_ym = np.zeros((max_calc_deg, 2 * max_calc_deg), dtype=complex)
    
            for l in range(0, max_calc_deg + 1):
                # Placeholder for P_lm[l], which should be of shape (l+1, N)
                # P_lm_l = self.P_lm[l]  # You need to define P_lm appropriately
                # Transpose and select first l+1 rows
                P_lm_l_T = P_lm_hd[l,:l+1, :].T  # Shape: (N, l+1)
                # Get coefficients for degree l
                coeffs = get_coeffs(coeff_hd, l)  # Define get_coeffs based on your data structure
                # Expand coefficients to shape (N, l+1)
                coeffs_expanded = np.tile(coeffs, (max_calc_deg, 1))
                # Multiply element-wise and accumulate
                F_ym[:, :l+1] += P_lm_l_T * coeffs_expanded

            coeff_matrix = np.zeros((max_calc_deg + 1, max_calc_deg + 1), dtype=complex)
            # idx = 0
            # for l in range(max_calc_deg + 1):
            #     coeff_matrix[l, :l + 1] = coeff_hd[idx:idx + l + 1]
            #     idx += l + 1

            # # Multiply P_lm_hd with coefficients across all degrees
            # # P_lm_hd shape: (max_calc_deg + 1, l+1, N)
            # # coeff_matrix shape: (max_calc_deg + 1, l+1)
            # # Result shape: (max_calc_deg + 1, N)
            # F_ym = np.sum(P_lm_hd[:, :, :] * coeff_matrix[:, :, None], axis=1)

            # # Enforce conjugate symmetry for FFT
            # F_ym = np.hstack([F_ym, np.conj(F_ym[:, -2:0:-1])])

            # # Perform inverse FFT
            # F_yx = ifft(F_ym, axis=1)

            # Use inverse FFT to sum over exponentials
            # F_ym[:, self.maxdeg:] = np.conj(F_ym[:, self.maxdeg-2::-1])
            F_ym[:, max_calc_deg+1:] = np.conj(F_ym[:, max_calc_deg-1:0:-1])
            F_yx = np.fft.ifft(F_ym, axis=1)

            # Apply additional factors
            F_yx *= (2 * max_calc_deg)#np.sqrt(2)
            grd_hd=np.real(F_yx.copy())

            
        return grd_hd,lon_hd,lat_hd
    
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
    
    def calc_Plm(self,deg):
        x,w=pysh.expand.SHGLQ(deg-1)
        P_lm=np.zeros((deg+1,deg+1,len(x)))
        for i,x_i in enumerate(x) :
                P_lm[:,:,i]=pysh.legendre.legendre(deg,x_i,csphase=1,cnorm=1)#*correction
        return P_lm
    
    # def calc_sphfunc_atpoint(self,lon,lat):
    #     theta=(lat+90)*np.pi/180
    #     phi=lon*np.pi/180  
    #     P_lm=pysh.legendre.legendre(self.maxdeg,np.cos(theta),csphase=1,cnorm=1,packed=False)
    #     # Y_lm=P_lm*np.exp(phi*np.repeat(np.arange(0,self.maxdeg+1),np.arange(1,self.maxdeg+2))*1j)
    #     Y_lm=np.zeros(P_lm.shape,dtype=np.complex128)
    #     for m in range(self.maxdeg+1):
    #         Y_lm[m,:]=P_lm[m,:]*np.exp(phi*np.arange(0,self.maxdeg+1)*1j)
    #     return Y_lm
    
    # def calc_at_point(self,lon,lat,Y_lm=None):
    #     if Y_lm is None :
    #         Y_lm = self.calc_sphfunc_atpoint(lon,lat)
    #     # result=np.sum(self.coeff*Y_lm)
    #     result=self.coeff.copy()
    #     ind_a=0
    #     for l in range(self.maxdeg+1):
    #         result[ind_a:ind_a + l + 1] = self.coeff[ind_a:ind_a + l + 1]*Y_lm[l,:l+1]
    #         ind_a += l + 1
    #     return np.real(result.sum())*2

    def calc_sphfunc_atpoint(self, lon, lat):
            """Calcule les spherical harmonics 4π-normalized à un point donné."""
            # Conversion de la latitude en colatitude
            theta = (lon) * np.pi / 180
            phi = (lat+90) * np.pi / 180

            # Calcul des fonctions associées de Legendre
            P_lm = pysh.legendre.legendre(self.maxdeg, np.cos(phi), csphase=1, cnorm=0, packed=False)

            # Calcul des spherical harmonics complexes
            Y_lm = np.zeros(P_lm.shape, dtype=np.complex128)
            for l in range(self.maxdeg + 1):
                for m in range(l + 1):
                    Y_lm[l, m] = P_lm[l, m] * np.exp(1j * m * theta)
            return Y_lm

    def calc_at_point(self, lon, lat, Y_lm=None,deg=None):
        # Convert geographic coordinates to spherical (radians)
        # Note: colatitude = 90° - lat
        theta = np.radians(lat+90)  
        phi = np.radians(lon)

        # Compute associated Legendre functions at cos(theta)
        P_lm = pysh.legendre.legendre(self.maxdeg, np.cos(theta), csphase=1, cnorm=1)

        f_val = 0.0
        ind = 0
        if deg is None :
            deg=self.maxdeg
        for l in range(deg + 1):
            for m in range(l + 1):
                # Compute the complex spherical harmonic value:
                
                Y_lm = P_lm[l, m] * np.exp(1j * m * phi)
                # If m == 0, it contributes directly.
                if m == 0:
                    f_val += (self.coeff[ind]  * Y_lm).real
                else:
                    # For real functions, the coefficient for negative m is 
                    # the conjugate of that for positive m. Thus the contribution is 2 * Re(coeff * Y_lm)
                    f_val += 2 * (self.coeff[ind] * Y_lm).real
                ind += 1
        return f_val
    
    import numpy as np

    def calc_multi_point(self, points):
        """
        Compute the value of the function at multiple points.

        Parameters:
            points (list of tuples): List of (lon, lat) points where the function will be evaluated.

        Returns:
            np.ndarray: Array of computed values at the given points.
        """
        results = np.zeros((len(points),))

        for i, (lon, lat) in enumerate(points):
            results[i] = self.calc_at_point(lon, lat)

        return results


            