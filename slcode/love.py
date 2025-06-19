import numpy as np
from numpy.matlib import repmat #used to add repmat 
from .spharm import sphericalobject
import math
import pandas as pd
from scipy import io
import matplotlib.pyplot as plt
import re
import os



def love_lm(num,maxdeg):
    '''
    the _`love_lm` funtion get from love numbers the h_lm spherical coefficient. 

    Attribute 
    ---------
        num : np.array([n,]) 
            LLN coefficient of the size of the entry file
        maxdeg : int
            The maximum harmonic coefficient degree. 

    Returns
    -------
        h_lm : np.array([(maxdeg+1)(maxdeg+2)/2,])
            Array of the love number repeated on harmonic degree orders. 

    '''
    h = np.concatenate(([0], num)) # add zeros as the first order of the LLN. 

    # Génération de h_lm
    h_lm = np.concatenate([np.repeat(h[n], n+1) for n in range(maxdeg)]) #repeat the l LLN upon the m (from 0 to l-1). 

    return h_lm 

def get_tlm(maxdeg,a,Me):
    '''
    The _`get_lm` function generate the T spherical harmonic coefficient as defined in :ref:`Theory <T_definition>`.

    Attribute
    --------- 
        maxdeg : int
            maximum degree of spherical harmonic to compute T_lm. 
        a : float
            The earth radius in meters. 
        Me : float
            The earth mass.

    Returns
    ------- 
        T_lm : np.array([(maxdeg+1)(maxdeg+2)/2])
            The T harmonic coefficient 

    '''
    T_lm=np.array([]) # preset the output
    T = np.zeros(maxdeg+1) # preset an array of the size maxdeg who will obtaine the coefficient and then be added to T_lm
    const = 4*np.pi*a**3/Me # setting the tide constant for earth
    for n in range(maxdeg+1) :
            T[n]=const/(2*n+1) # for each nth add to T the earth constant modulated by the nth step
            T_add=repmat(T[n],1,n+1) # prepare T_add as the repetition of T
            if n==0 :
                T_lm=T_add # if there is nothing in T_lm preset T_lm as T_add
            else :
                T_lm=np.concatenate((T_lm,T_add),1) # else add T_add to T_lm
    return np.squeeze(T_lm) # we have to squeeze the array so that the requested indices are directly on the good axis

def calc_beta_counter(self,maxdeg):
    '''
    The _`calc_beta_counter` define the indices of spherical harmonics order towards the l degree must be replecated.
    
    Attribute
    ---------
        self : :ref:`LOVE <LOVE>` class object
            The LOVE class object on wich the betacounter is calculated
        maxdeg : int
            maximum spherical harmonic coefficient
    
    Returns
    -------
        None

    '''
    self.beta_counter=np.repeat(np.arange(0,maxdeg),np.arange(1,maxdeg+1))#create an array where it's values are the indices towards wich an array of LLN upon the l degree are reproduce to the m order. 
    return

class LOVE(object):
    """
    The _`LOVE` class is used to keep the love numbers values and prepare them for the computation of geoïd and ground vertical motion. The love number are calculted and loaded from a file as described in :ref:`Implementation of Love numbers <love>`. This class also include the possibility to compute the love numbers from normal modes love numbers parameters. 

    Methods
        -------
            `dev_beta`_ :
                This method is used to calculate the beta love numbers upon spherical harmonic order m.
            `dev_beta_tide`_ :
                This method is used to calculate the beta tidal love numbers upon spherical harmonic order m.
            `clean_memory`_ :
                This method is used to clean the memory. 
    """
    
    def __init__(self,maxdeg=512,time_step=np.arange(122,-1,-1),way=None,a=6371000,Me=5.9742e24,type='time',T=None):
        """
        The _`__init__` function 

        Attributes
        ----------
            maxdeg : int
                The spherical harmonic maximum degree on wich the LLN are loaded/computed.
            way : str
                The file path to the LLN. If the love_type is 'normal' it must be a .mat file that contain 10 elements : 'mode_found', 'k_amp', 'h_amp', 'k_amp_tide', 'h_amp_tide', 'spoles', 'h_el', 'k_el', 'h_el_tide', 'k_el_tide'. To know more on this type of love numbers see ... . If the love_type is 'time', the folder must constain the following love numbers files : 'h_e', 'k_e', 'k_ve', 'h_ve', 'h_e_T', 'k_e_T', 'h_ve_T', 'k_ve_T', 'time'. 
            time_step : np.array([time_step_number,])
                The time step of the model, used to compute the LLN.
            a : float
                The earth radius in meter.
            Me : float
                The earth mass.
            type : str
                The type of love number in input. could be 'time' or 'normal', where 'time' is for love numbers from ALMA3 code, 'normal' is for love number in normal mode from MIT server can be found at :url:`https://github.com/jaustermann/SLcode/tree/master/SavedLN`. Default is 'time'.
        """
        self.maxdeg=maxdeg
        if not(way is None):
            self.type=type
            self.time_step=time_step
            self.maxdeg=maxdeg
            self.time_step_number=len(time_step)
            # Load the Load Love Numbers
            if type is 'time': # case where the LLN are computed from ALMA3
                # Loading the LLN :
                #Load the elastic love numbers
                self.h_e=np.loadtxt(way+'/h_e.dat',unpack=True)[1,:maxdeg+1]
                self.k_e=np.loadtxt(way+'/k_e.dat',unpack=True)[1,:maxdeg+1]
                
                #Load the viscous love numbers
                self.k_ve=np.loadtxt(way+'/k_ve.dat',unpack=True)[:,:]
                self.h_ve=np.loadtxt(way+'/h_ve.dat',unpack=True)[:,:]

                #Load the elastic tidal love numbers
                self.h_tide_e=np.loadtxt(way+'/h_e_T.dat',unpack=True)[1,:6]
                self.k_tide_e=np.loadtxt(way+'/k_e_T.dat',unpack=True)[1,:6]

                #Load the viscous tidal love numbers
                self.k_tide_ve=np.loadtxt(way+'/k_ve_T.dat',unpack=True)[:,:6]
                self.h_tide_ve=np.loadtxt(way+'/h_ve_T.dat',unpack=True)[:,:6]

                #Load the time of the LLN
                self.love_time=np.loadtxt(way+'/time.dat',unpack=True)

                # locate the time indices in love_time and the corresponding indices of the derivation of time step. 
                time_step_diff=-(time_step.reshape(-1,1)-time_step)
                time_step_diff[time_step_diff<=0]=0 #create the derivative of time_step
                diff_matrix = np.abs(time_step_diff[:, :, np.newaxis] - self.love_time)
                data = np.argmin(diff_matrix, axis=2) #Determine the indices corresponding to the diff_matrix. 

                #Develop the elastic LLN.
                self.k_e[0]=0
                self.h=love_lm(self.h_e,maxdeg+1)
                self.k=love_lm(self.k_e,maxdeg+1)

                #Transform the viscous geoid LLN following the data matice containing the indices of the corresponding time_diff.
                self.beta_G_l=self.k_ve[data,:maxdeg+1].squeeze()-np.repeat(self.k_e[:,np.newaxis],self.time_step_number,1).T#Retreave to the modified viscous LLN the elastic part.
                self.beta_G_l[data<=0]=self.beta_G_l[data<=0]*0# Set to 0 all Geoid viscous LLN that are associated witha negative time_diff. 
                self.beta_G_l=self.beta_G_l[:-1,:-2,:]#Correct the shape of the love numbers

                #Transform the viscous solid earth LLN following the data matice containing the indices of the corresponding time_diff.
                self.beta_R_l=self.h_ve[data,:maxdeg+1].squeeze()-np.repeat(self.h_e[:maxdeg+1,np.newaxis],self.time_step_number,1).T#Retreave to the modified viscous LLN the elastic part.
                self.beta_R_l[data<=0]=self.beta_R_l[data<=0]*0# Set to 0 all Geoid viscous LLN that are associated witha negative time_diff.
                self.beta_R_l=self.beta_R_l[:-1,:-2,:]#Correct the shape of the love numbers

                self.beta_l=self.beta_G_l-self.beta_R_l# Calculate the total viscuous LLN for sea level, i.e. G-R
                self.beta_konly_l=-(self.k_ve[data,1]-self.k_e[1]-self.h_ve[data,1]+self.h_e[1])# Compute the konly LLN : the 1 order beta component. 
                
                # Note that, as the polar wonderer computation is only on the 3 first harmonics degree we compute only the tidal LLN upon 6 degree of spherical harmonics.
                #Develop the tidal elastic LLN. 
                self.k_tide_e[0]=0
                self.h_tide_e[0]=0
                self.h_tide=love_lm(self.h_tide_e,6)
                self.k_tide=love_lm(self.k_tide_e,6)

                #Transform the tidal viscous geoid LLN following the data matice containing the indices of the corresponding time_diff.
                self.beta_G_l_tide=self.k_tide_ve[data,:6].squeeze()-np.repeat(self.k_tide_e[:6,np.newaxis],self.time_step_number,1).T#Retreave to the modified viscous LLN the elastic part.
                self.beta_G_l_tide[data<=0]=self.beta_G_l_tide[data<=0]*0# Set to 0 all Geoid viscous LLN that are associated witha negative time_diff. 
                self.beta_G_l_tide=self.beta_G_l_tide[:,:-2,:]#Correct the shape of the love numbers
                self.beta_G_l_tide=np.concatenate(((self.beta_G_l_tide[:,:,1]*0)[:,:,np.newaxis],self.beta_G_l_tide),axis=2)[:,:,:6]

                #Transform the tidal viscous solid earth LLN following the data matice containing the indices of the corresponding time_diff.
                self.beta_R_l_tide=self.h_tide_ve[:,:6]-np.repeat(self.h_tide_e[:6,np.newaxis],self.h_tide_ve.shape[0],1).T#Retreave to the modified viscous LLN the elastic part.
                self.beta_R_l_tide=self.beta_R_l_tide[data,:6].squeeze()
                self.beta_R_l_tide[data<=0]=self.beta_R_l_tide[data<=0]*0# Set to 0 all Geoid viscous LLN that are associated witha negative time_diff.
                self.beta_R_l_tide=self.beta_R_l_tide[:,:-2,:]#Correct the shape of the love numbers
                self.beta_R_l_tide=np.concatenate(((self.beta_R_l_tide[:,:,1]*0)[:,:,np.newaxis],self.beta_R_l_tide),axis=2)[:,:,:6]

                self.beta_l_tide=(self.beta_G_l_tide-self.beta_R_l_tide)# Calculate the tidal total viscuous LLN for sea level, i.e. G-R
                self.beta_konly_l_tide=-(self.k_tide_ve[data,1]-self.k_tide_e[1]-self.h_tide_ve[data,1]+self.h_tide_e[1])# Compute the konly tidal LLN : the 1 order beta component. 
                
                # Compute the elastic Load and Tidal LN for solid earth and geoid. 
                self.E = 1+self.k - self.h
                self.E_T = 1 + self.k_tide - self.h_tide
                self.T = sphericalobject(coeff=get_tlm(maxdeg,a,Me))#Compute T the earth constant correction.

                calc_beta_counter(self,maxdeg+1)# Compute the corresponding indices for each degree to the order associated to degree. 
                self.beta_konly_l=np.array(self.beta_konly_l)
                self.beta_konly_l_tide=np.array(self.beta_konly_l_tide)

            elif type is 'normal': #Case where the LLN are in normal mode. 
                tide_size=6 # define the length of tidal love numbers. As the rotational potential deduced from the polar wanderer is affecting on ly the 3 first degree of the spherical harmonics we just need to compute a small part of the tidal LLN. 
                # print(way)
                love = io.loadmat(way)# loading love numbers file from a .mat file.
                self.mode_found=love['mode_found']
                self.k_amp=love['k_amp'][:maxdeg+1,:]
                self.h_amp=love['h_amp'][:maxdeg+1,:]
                self.k_amp_tide=love['k_amp_tide'][:maxdeg,:]
                self.h_amp_tide=love['h_amp_tide'][:maxdeg,:]
                self.spoles=love['spoles'][:maxdeg+1,:]
                self.h_e=love['h_el'][:maxdeg+1].squeeze()
                self.k_e=love['k_el'][:maxdeg+1].squeeze()
                self.h_tide_e=love['h_el_tide'][:tide_size].squeeze()
                self.k_tide_e=love['k_el_tide'][:tide_size].squeeze()
                
                # Preparing a the elastic love numbers
                self.h = love_lm(self.h_e,maxdeg+1)
                self.k = love_lm(self.k_e,maxdeg+1)
                self.h_e=self.h_e.squeeze()
                self.k_e=self.k_e.squeeze()

                # Preparing a the tidal elastic love numbers
                self.h_tide = love_lm(self.h_tide_e,tide_size)
                self.k_tide = love_lm(self.k_tide_e,tide_size)

                # Presetting the viscous LLN as zeros matrices.
                self.beta_l = np.zeros((self.time_step_number-1,self.time_step_number-2,maxdeg+1))
                self.beta_G_l = np.zeros((self.time_step_number-1,self.time_step_number-2,maxdeg+1))
                self.beta_R_l = np.zeros((self.time_step_number-1,self.time_step_number-2,maxdeg+1))
                self.beta_konly_l = np.zeros((self.time_step_number-1,self.time_step_number-2))
                self.beta_G_konly_l = np.zeros((self.time_step_number-1,self.time_step_number-2))
                self.beta_R_konly_l = np.zeros((self.time_step_number-1,self.time_step_number-2))

                # computing the viscous LLN from theyre normal mode to theyre time shape
                for t_it in range(1,self.time_step_number): # for each time spet we compute the LLN
                    for n in range(1,t_it):# for each time step we compute the LLN for the differential time step self.time_step[n]-self.time_step[t_it]
                        beta = np.zeros((maxdeg,))
                        beta_G = np.zeros((maxdeg,))
                        beta_R = np.zeros((maxdeg,))
                        delta_time = -self.time_step[t_it] + self.time_step[n]
                        for lm in range(maxdeg):# The computation is applied for each degree of spherical harmonics up to maximum degree
                            num_mod = self.mode_found.squeeze()[lm]
                            numerator = self.k_amp[lm, :num_mod] - self.h_amp[lm, :num_mod]
                            denominator = self.spoles[lm, :num_mod]
                            exponentials = 1 - np.exp(-denominator * delta_time)
                            terms = numerator / denominator * exponentials
                            beta[lm] = np.sum(terms)
                            beta_G[lm]=np.sum((self.k_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                            beta_R[lm]=np.sum((self.h_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                        self.beta_l[t_it-1,n-1,:]=np.concatenate((np.zeros((1,)),beta))
                        self.beta_G_l[t_it-1,n-1,:]=np.concatenate((np.zeros((1,)),beta_G))
                        self.beta_R_l[t_it-1,n-1,:]=np.concatenate((np.zeros((1,)),beta_R))
                        lm=1
                        num_mod=self.mode_found[lm][0]
                        self.beta_konly_l[t_it-1,n-1]=np.sum((self.k_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n])))) # Compute the konly LLN : the 1 order beta component. 
                
                # Presetting the tidal viscous LLN as zeros matrices.
                self.beta_l_tide = np.zeros((self.time_step_number-1,self.time_step_number-2,tide_size))
                self.beta_G_l_tide = np.zeros((self.time_step_number-1,self.time_step_number-2,tide_size))
                self.beta_R_l_tide = np.zeros((self.time_step_number-1,self.time_step_number-2,tide_size))
                self.beta_konly_l_tide = np.zeros((self.time_step_number-1,self.time_step_number-2))
                self.beta_G_konly_l_tide = np.zeros((self.time_step_number-1,self.time_step_number-2))
                self.beta_R_konly_l_tide = np.zeros((self.time_step_number-1,self.time_step_number-2))

                # computing the tidal viscous LLN from theyre normal mode to theyre time shape
                for t_it in range(1,self.time_step_number): # for each time spet we compute the tidal LLN
                    for n in range(1,t_it):# for each time step we compute the tidal LLN for the differential time step self.time_step[n]-self.time_step[t_it]
                        beta_tide = np.zeros((tide_size-1,))
                        beta_G_tide = np.zeros((tide_size-1,))
                        beta_R_tide = np.zeros((tide_size-1,))
                        for lm in range(tide_size-1): # The computation is applied for each degree of spherical harmonics up to maximum degree
                            num_mod = self.mode_found[lm][0] 
                            beta_tide[lm] = np.sum((self.k_amp_tide[lm,:num_mod] - self.h_amp_tide[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                            beta_G_tide[lm]=np.sum((self.k_amp_tide[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                            beta_R_tide[lm]=np.sum((self.h_amp_tide[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                        self.beta_l_tide[t_it-1,n-1,:]=np.concatenate((np.zeros((1,)),beta_tide))
                        self.beta_G_l_tide[t_it-1,n-1,:]=np.concatenate((np.zeros((1,)),beta_G_tide))
                        self.beta_R_l_tide[t_it-1,n-1,:]=np.concatenate((np.zeros((1,)),beta_R_tide))
                        lm=1
                        num_mod=self.mode_found[lm][0]
                        self.beta_konly_l_tide[t_it-1,n-1]=np.sum((self.k_amp_tide[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))# Compute the konly tidal LLN : the 1 order beta component. 

                # Compute the elastic Load and Tidal LN for solid earth and geoid. 
                self.E = 1+self.k - self.h
                self.E_T = 1 + self.k_tide - self.h_tide
                self.T = sphericalobject(coeff=get_tlm(maxdeg,a,Me)) #Compute T the earth constant correction.

                calc_beta_counter(self,maxdeg+1)# Compute the corresponding indices for each degree to the order associated to degree.

            

    def dev_beta(self,applied='beta'):
        '''
        The _`dev_beta` method can be used to rise the love numbers to their full shape to fit the computation method. This is required before using the beta in the methods from :ref:`LOAD <LOAD>`.

        Attribute
        ---------
            applied : str
                The kinf of beta you are calculating. Three values are possible : beta, beta_R et beta_G. 'beta' define the love numbers used to caculate the variation of ocean thickness. 'beta_R' define the love numbers used to calculate the variation of groud. 'beta_G' define the love numbers used to calculate the variation of geoïd. 
        Return
        ------
            None

        '''
        if self.type is 'time':
            if applied is 'beta':
                self.beta_G=0
                self.beta_R=0
                # self.beta_l=self.beta_G_l-self.beta_R_l
                self.beta_l=np.array(self.beta_l)[:,:,self.beta_counter.astype(int)]
            elif applied is 'beta_R':
                self.beta_l=0
                self.beta_G=0
                self.beta_R=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)]
            elif applied is 'beta_G':
                self.beta_l=0
                self.beta_R=0
                self.beta_G=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)]
        elif self.type is "normal" :
            if applied is 'beta':
                self.beta_G=0
                self.beta_R=0
                self.beta_l=np.array(self.beta_l)[:,:,self.beta_counter.astype(int)]
            elif applied is 'beta_R':
                self.beta_l=0
                self.beta_G=0
                self.beta_R=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)]
            elif applied is 'beta_G':
                self.beta_l=0
                self.beta_R=0
                self.beta_G=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)]

    def dev_beta_tide(self,applied='beta'):
        '''
        The _`dev_beta_tide` method can be used to rise the love numbers to their full shape to fit the computation method. This is required before using the beta in the methods from :ref:`LOAD <LOAD>`.

        Attribute
        ---------
            applied : str
                The kinf of beta you are calculating. Three values are possible : beta, beta_R et beta_G. 'beta' define the love numbers used to caculate the variation of ocean thickness. 'beta_R' define the love numbers used to calculate the variation of groud. 'beta_G' define the love numbers used to calculate the variation of geoïd. 
        Return
        ------
            None

        '''
        if self.type is 'time':
            if applied is 'beta':
                self.beta_G_tide=0
                self.beta_R_tide=0
                # self.beta_l_tide=self.beta_G_l_tide-self.beta_R_l_tide
                self.beta_l_tide=np.array(self.beta_l_tide)[:,:,self.beta_counter.astype(int)[:6]]
            elif applied is 'beta_R':
                self.beta_l_tide=0
                self.beta_G_tide=0
                self.beta_R_tide=np.array(self.beta_R_l_tide)[:,:,self.beta_counter.astype(int)[:6]]
            elif applied is 'beta_G':
                self.beta_l_tide=0
                self.beta_R_tide=0
                self.beta_G_tide=np.array(self.beta_G_l_tide)[:,:,self.beta_counter.astype(int)[:6]]
        elif self.type is "normal" :
            if applied is 'beta':
                # print(np.array(self.beta_l_tide)[:,:,self.beta_counter.astype(int)[:int(6*(6+1)/2)]])
                self.beta_G_tide=0
                self.beta_R_tide=0
                self.beta_l_tide=np.array(self.beta_l_tide)[:,:,self.beta_counter.astype(int)[:6]]
            elif applied is 'beta_R':
                self.beta_l_tide=0
                self.beta_G_tide=0
                self.beta_R_tide=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)[:6]]
            elif applied is 'beta_G':
                self.beta_l_tide=0
                self.beta_R_tide=0
                self.beta_G_tide=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)[:6]]


    def write_layer_file(self, filename, thickness, rigidity, density, viscosity, material_type,
                     earth_radius=6371e3, force_elastic=False, thickness_is_altitude=False):
        """
        Writes a file describing Earth layers and returns a DataFrame of the layer parameters.

        Parameters:
        - filename: path to output file
        - thickness: array of layer thicknesses (m) or altitudes/radii if thickness_is_altitude=True
        - rigidity: array of elastic moduli (Pa)
        - density: array of densities (kg/m³)
        - viscosity: array of viscosities (Pa·s)
        - material_type: array of strings (e.g., 'elastic', 'maxwell', 'fluid')
        - earth_radius: reference radius from surface (default is 6371e3 m)
        - force_elastic: if True, override all material_type with 'elastic' and write again
        - thickness_is_altitude: if True, interpret `thickness` as absolute radius, not layer thickness

        Returns:
        - pd.DataFrame containing all parameters, including computed radius
        """
        assert len(thickness) == len(rigidity) == len(density) == len(viscosity) == len(material_type), \
            "Arrays must have the same length"
        if force_elastic:
            self.write_layer_file(
                filename+'_elastic',
                thickness,
                rigidity,
                density,
                viscosity,
                material_type=['elastic'] * (len(material_type)-1)+['fluid'],
                earth_radius=earth_radius,
                force_elastic=False,
                thickness_is_altitude=thickness_is_altitude
            )

        if thickness_is_altitude:
            radius = np.array(thickness)
        else:
            radius = earth_radius - np.cumsum(thickness)

        # Build DataFrame to return
        self.earth = pd.DataFrame({
            'radius (m)': radius,
            'density (kg/m³)': density,
            'rigidity (Pa)': rigidity,
            'viscosity (Pa.s)': viscosity,
            'material_type': material_type
        })

        with open(filename+'.dat', 'w') as f:
            f.write('!------------------------------------------------------------\n')
            f.write('! radius,    density,      rigidity     viscosity\n')
            f.write('!  (m)       (kg/m^3)        (Pa)         (Pa.s) \n')
            f.write('!------------------------------------------------------------\n')
            for _, row in self.earth.iterrows():
                f.write(f'{row["radius (m)"]:10.1e}    {row["density (kg/m³)"]:10.5e}    '
                        f'{row["rigidity (Pa)"]:10.4e}    {row["viscosity (Pa.s)"]:10.2e}      {row["material_type"]}\n')
                

    

    def write_config_file(self, template_dir,output_config, output_dir, model_way, final_time, time_step,
                      earth_file=None, extra_params=False):
        """
        Generate multiple config files from templates, filling in placeholders.

        Parameters:
        - template_dir: folder containing the .dat templates
        - output_dir: folder where to save the modified config files
        - model_way: path to the model executable or identifier
        - final_time: float, simulation final time (e.g., in kyr)
        - time_step: float, time step (same units as final_time)
        - earth_file: path to the Earth model file (.dat)
        - extra_params: optional dict of extra key-value pairs to substitute

        Returns:
        - List of paths to written config files
        """
        templates = ['config.EarthLLNs.ELLN', 'config.EarthLLNs.ETLN',
                    'config.EarthLLNs.VLLN', 'config.EarthLLNs.VTLN']
        written_files = []

        os.makedirs(output_config, exist_ok=True)  # Ensure output directory exists
        os.makedirs(f'{output_dir}/{earth_file}', exist_ok=True) 

        for config_name in templates:
            input_path = os.path.join(template_dir, config_name + '.dat')
            output_path = os.path.join(output_config, config_name + '.dat')

            with open(input_path, 'r') as f:
                template = f.read()

            placeholders = set(re.findall(r"\$(.*?)\$", template))

            values = {
                'maximum_degree': f"{self.maxdeg}",
                'time_end': f"{np.log10(final_time)}",
                'time_start': f"{np.log10(time_step)}",
                'time_step_number': f"{int(final_time / time_step)}",
                'model_way': model_way
            }

            if earth_file:
                values['model_name'] = earth_file

            values['output_way'] = output_dir

            if self.earth is not None:
                try:
                    values['layers_number'] = str(len(self.earth))
                except KeyError:
                    pass

            if extra_params:
                for k, v in extra_params.items():
                    values[k] = str(v)

            for key in placeholders:
                if key in values:
                    template = template.replace(f"${key}$", values[key])

            with open(output_path, 'w') as f:
                f.write(template)

            written_files.append(output_path)

        return written_files
    
    def append_to_bash_script(self,script_path, config_path, model_name):
        if not(os.path.isfile(script_path)):
            with open(script_path, 'w') as bash_file:
                bash_file.write("#!/bin/bash\n\n")
        else :
            with open(script_path, 'a') as bash_file:
                bash_file.write(f"./alma.exe {config_path}/{model_name}/config.EarthLLNs.ELLN.dat\n")
                bash_file.write(f"./alma.exe {config_path}/{model_name}/config.EarthLLNs.ETLN.dat\n")
                bash_file.write(f"./alma.exe {config_path}/{model_name}/config.EarthLLNs.VLLN.dat\n")
                bash_file.write(f"./alma.exe {config_path}/{model_name}/config.EarthLLNs.VTLN.dat\n")




    def plot_layers(self, parameters=None, ax=None, figsize=(6, 8), colors=None,earth_radius=6371e3):
        """
        Plot selected Earth layer parameters as horizontal-layer steps.

        Parameters:
        - parameters: list of column names to plot (e.g., ['density (kg/m³)', 'viscosity (Pa.s)'])
        - ax: existing matplotlib axis to plot on (optional)
        - figsize: figure size if ax is None
        - colors: dict of colors for parameters, or list of colors

        Returns:
        - fig, ax: matplotlib figure and axes
        """
        if self.earth is None:
            raise ValueError("No layer data to plot. Please run write_layer_file first.")

        import numpy as np

        df = self.earth.copy()
        radius = df['radius (m)'].values
        depth = earth_radius - radius  # Compute depth from radius

        if parameters is None:
            parameters = ['density (kg/m³)', 'rigidity (Pa)', 'viscosity (Pa.s)']

        # Prepare stepped data
        def stepped(x):
            return np.repeat(x, 2)

        depth_edges = np.concatenate(([0], depth))
        depth_steps = stepped(depth_edges)
        depth_steps=depth_steps[1:]
        # depth_steps=np.concatenate((depth_steps, depth_steps[-1].flatten()))

        fig, ax = plt.subplots(figsize=figsize) if ax is None else (plt.gcf(), ax)

        for i, param in enumerate(parameters):
            values = df[param].values
            # values_edges = np.concatenate(([values[0]], values))  # duplicate first value for step start
            value_steps = stepped(values)
            # value_steps=value_steps[:-1]
            # value_steps = value_steps
            value_steps=np.concatenate((value_steps, value_steps[-1].flatten()))

            color = None
            if colors:
                if isinstance(colors, dict):
                    color = colors.get(param, None)
                elif isinstance(colors, list):
                    color = colors[i % len(colors)]

            ax.plot(value_steps, depth_steps, label=param, color=color)

        ax.invert_yaxis()
        ax.set_xlabel("Value")
        ax.set_ylabel("Depth (m)")
        ax.legend()
        ax.grid(True)

        return fig, ax


    
    def clean_memory(self):
        '''
        The _`clean_memory` method can be used to araise the beta_l, beta_R and beta_G to avoid memory issues. 

        Attribute
        ---------
            None

        Return
        ------
            None

        '''
        self.beta_l=0
        self.beta_R=0
        self.beta_G=0
