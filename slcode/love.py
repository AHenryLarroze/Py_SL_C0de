import numpy as np
from numpy.matlib import repmat #used to add repmat 
from .spharm import sphericalobject
import math
from scipy import io



def love_lm(num,maxdeg):
    '''
    the _`love_lm` funtion get from love numbers the h_lm spherical coefficient. 

    Attribute 
    ---------
        num : np.array([n,]) 
            Love number coefficient of the size of the entry file
        maxdeg : int
            The maximum harmonic coefficient degree. 

    Returns
    -------
        h_lm : np.array([(maxdeg+1)(maxdeg+2)/2,])
            Array of the love number repeated on harmonic degree orders. 

    '''
    num=np.concatenate((np.zeros((1,1)),num)) 
    h_lm=np.repeat(num,np.arange(1,maxdeg+2,1))
    return h_lm 

def get_tlm(maxdeg,a,Me):
    '''
    The _`get_lm` function generate the T spherical harmonic coefficient as defined in :ref:`Theory <T_definition>`.

    Attribute
    --------- 
        maxdeg : int
            maximum degree of spherical harmonic defined in the model parameters. 
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
    const = 4*np.pi*a**3/Me # setting the tide constant for earth !!!! CHECK !!!!
    for n in range(maxdeg+1) :
            T[n]=const/(2*n+1) # for each nth add to T the earth constant moduladed by the nth step
            T_add=repmat(T[n],1,n+1) # prepare T_add as the repetition of T
            if n==0 :
                T_lm=T_add # if there is nothing in T_lm preset T_lm as T_add
            else :
                T_lm=np.concatenate((T_lm,T_add),1) # else add T_add to T_lm
    return np.squeeze(T_lm) # we have to squeeze the array so that the requested indices are directly on the good axis

def calc_beta_counter(self,maxdeg):
    '''
    The _`calc_beta_counter` define the degree of the spherical harmonic for each beta coefficient 
    
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
    # self.beta_counter = np.ones((self.h.shape[0]-1,))
    self.beta_counter=np.repeat(np.arange(0,maxdeg),np.arange(1,maxdeg+1))
    # l_it = 1
    # for lm_it in range(1,len(self.h)-1):
    #     # for each time step if the coefficient indices is one degree then , all the next values in beta_counter will have the previous of this one +1. litteraly the degree associated to there indices.
    #     if lm_it == l_it*(l_it+1)/2:
    #         self.beta_counter[lm_it] = self.beta_counter[lm_it-1]+1
    #         l_it = l_it+1
    #     else:
    #         self.beta_counter[lm_it] = self.beta_counter[lm_it-1]
    return

class LOVE(object):
    """
    The _`LOVE` class is used to keep the love numbers values and prepare them for the computation of geoïd and ground vertical motion. The love number are calculted and loaded from a file as described in :ref:`Implementation of Love numbers <love>`. This function also include the possibility to compute the love numbers from normal modes love numbers parameters. 

    Attributes
    ----------
        maxedeg : int
            The spherical harmonic maximum degree.
        way : str
            The file path to the ALMA output file. This file must follow the described pattern in ...
        time_step : np.array([time_step_number,])
            The time step of the model, used to prepare the viscuous love numbers.
        a : float
            The earth radius in meter.
        Me : float
            The earth mass.
        type : str
            The type of love number in input. could be 'time', 'normal' or 'spectral', where 'time' is for love numbers from ALMA3 code, 'normal' is for love number in normale mode from MIT server and 'spectral' is the frequential decomposition of the Love numbers available in ALMA3. Default is 'time'.
    
    Methods
    -------
        `dev_beta`_ :
            This method is used to calculate the beta love numbers.
        `dev_beta_tide`_ :
            This method is used to calculate the beta tidal love numbers.
        `clean_memory`_ :
            This method is used to clean the memory of your computer. 

    """
    def __init__(self,maxdeg,way,time_step,a,Me,type='time',T=None):
        """
    Parameters
    ----------
    grid : object (from class GRID)
    way : !!!! A corriger !!!! 
        """
        self.type=type
        self.time_step=time_step
        self.maxdeg=maxdeg
        self.time_step_number=len(time_step)
        # Load the Load Love Numbers
        if type is 'time': 
            self.h_e=np.loadtxt(way+'/h_e.dat',unpack=True)[1,:maxdeg]#/(4*math.pi)
            self.k_e=np.loadtxt(way+'/k_e.dat',unpack=True)[1,:maxdeg]#/(4*math.pi)
            self.k_e[0]=0
            self.h=np.repeat(self.h_e,np.arange(1,maxdeg+1,1))
            self.k=np.repeat(self.k_e,np.arange(1,maxdeg+1,1))
            self.k_ve=np.loadtxt(way+'/k_ve.dat',unpack=True)[1:,:]#/(4*math.pi)
            self.h_ve=np.loadtxt(way+'/h_ve.dat',unpack=True)[1:,:]#/(4*math.pi)

            self.love_time=np.loadtxt(way+'/time.dat',unpack=True)


            time_step_diff=-(time_step.reshape(-1,1)-time_step)
            time_step_diff[time_step_diff<=0]=0
            data=np.arange(0,len(self.love_time),1)
            data=np.repeat(repmat(data,len(time_step_diff),1)[np.newaxis,:,:],len(time_step_diff),axis=0)
            time_interval=(np.repeat(time_step_diff[:,:,np.newaxis],len(self.love_time),axis=2)-np.repeat(repmat(self.love_time,len(time_step_diff),1)[np.newaxis,:,:],len(time_step_diff),axis=0))==0
            data[time_interval==False]=0
            data=data.sum(2)

            # time_step_diff=time_step_diff.flatten()
            # time_step_diff=time_step_diff[time_step_diff>0]


            # time_interval=np.where((repmat(time_step_diff,len(self.love_time),1).transpose()-repmat(self.love_time,len(time_step_diff),1))==0)[1]
            # time_interval_tri=np.zeros((self.time_step_number,self.time_step_number))
            # time_interval_tri[np.triu_indices(self.time_step_number,1)]=time_interval
            # time_interval_tri=time_interval_tri.transpose().astype(int)

            k_ve=np.concatenate((np.zeros((self.k_ve.shape[0],1)),self.k_ve),axis=1)
            h_ve=np.concatenate((np.zeros((self.h_ve.shape[0],1)),self.h_ve),axis=1)
            k_e=np.concatenate((np.zeros((1,)),self.k_e))
            h_e=np.concatenate((np.zeros((1,)),self.h_e))
            self.beta_G_l=k_ve[data,:maxdeg].squeeze()+np.repeat(self.k_e[:maxdeg,np.newaxis],self.time_step_number,1).T
            self.beta_G_l[0,:,0]=0
            self.beta_G_l[0,0,:]=0
            #self.beta_G_l=k_ve[time_interval_tri+1,:maxdeg].squeeze()-k_ve[time_interval_tri,:maxdeg].squeeze()
            self.beta_G_l[data<=0]=self.beta_G_l[data<=0]*0
            self.beta_R_l=h_ve[:,:maxdeg]+np.repeat(self.h_e[:maxdeg,np.newaxis],h_ve.shape[0],1).T
            self.beta_R_l=self.beta_R_l[data,:maxdeg].squeeze()#-h_e[:maxdeg]
            #self.beta_R_l=h_ve[time_interval_tri+1,:maxdeg].squeeze()-h_ve[time_interval_tri,:maxdeg].squeeze()
            self.beta_R_l[data<=0]=self.beta_R_l[data<=0]*0
            self.beta_l=self.beta_G_l-self.beta_R_l
            self.beta_konly_l=self.k_ve[data,1]-self.k_e[0]#-(self.h_ve[data,1]-self.h_e[0])


            # Load the Tide Love Numbers
            self.h_tide_e=np.loadtxt(way+'/h_e_T.dat',unpack=True)[1,1:7]
            self.k_tide_e=np.loadtxt(way+'/k_e_T.dat',unpack=True)[1,1:7]
            self.h_tide=np.repeat(self.h_tide_e,np.arange(1,7,1))
            self.k_tide=np.repeat(self.k_tide_e,np.arange(1,7,1))
            self.k_tide_ve=np.loadtxt(way+'/k_ve_T.dat',unpack=True)[1:,1:7]
            self.h_tide_ve=np.loadtxt(way+'/h_ve_T.dat',unpack=True)[1:,1:7]
            #print(self.h_tide_e.shape,self.k_tide_e.shape,self.h_tide.shape,self.k_tide.shape,self.k_tide_ve.shape,self.h_tide_ve.shape)

            k_tide_ve=np.concatenate((np.zeros((self.k_tide_ve.shape[0],1)),self.k_tide_ve),axis=1)
            h_tide_ve=np.concatenate((np.zeros((self.h_tide_ve.shape[0],1)),self.h_tide_ve),axis=1)
            k_tide_e=np.concatenate((np.zeros((1,)),self.k_tide_e))
            h_tide_e=np.concatenate((np.zeros((1,)),self.h_tide_e))

            self.beta_G_l_tide=k_tide_ve[data,:6].squeeze()+np.repeat(self.k_tide_e[:6,np.newaxis],self.time_step_number,1).T
            self.beta_G_l_tide[0,:,0]=0
            self.beta_G_l_tide[0,0,:]=0
            #self.beta_G_l=k_ve[time_interval_tri+1,:maxdeg].squeeze()-k_ve[time_interval_tri,:maxdeg].squeeze()
            self.beta_G_l_tide[data<=0]=self.beta_G_l_tide[data<=0]*0
            self.beta_R_l_tide=h_tide_ve[:,:6]+np.repeat(self.h_tide_e[:6,np.newaxis],h_tide_ve.shape[0],1).T
            self.beta_R_l_tide=self.beta_R_l_tide[data,:7].squeeze()#-h_e[:maxdeg]
            #self.beta_R_l=h_ve[time_interval_tri+1,:maxdeg].squeeze()-h_ve[time_interval_tri,:maxdeg].squeeze()
            self.beta_R_l_tide[data<=0]=self.beta_R_l_tide[data<=0]*0
            self.beta_l_tide=self.beta_G_l_tide-self.beta_R_l_tide
            self.beta_konly_l_tide=self.k_tide_ve[data,1]-self.k_tide_e[0]#-(self.h_tide_ve[data,1]-self.h_tide_e[0])

            
            # self.beta_tide=k_tide_ve[time_interval_tri,:maxdeg+1].squeeze()-k_tide_e[:maxdeg+1]-(h_tide_ve[time_interval_tri,:maxdeg+1].squeeze()-h_tide_e[:maxdeg+1])
            # self.beta_konly_tide=self.k_tide_ve[time_interval_tri,1]-self.k_tide_e[0]-(self.h_tide_ve[time_interval_tri,1]-self.h_tide_e[0])
            
            self.E = 1+self.k - self.h
            self.E_T = 1 + self.k_tide - self.h_tide
            self.T = sphericalobject(coeff=get_tlm(maxdeg-1,a,Me))

            calc_beta_counter(self,maxdeg)
            #self.beta_G_l=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)]
            #self.beta_R_l=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)]
            self.beta_konly_l=np.array(self.beta_konly_l)
            #self.beta_tide=np.array(self.beta_tide)[:,:,self.beta_counter.astype(int)]
            self.beta_konly_l_tide=np.array(self.beta_konly_l_tide)
        elif type is 'normal':
            a=0
            #retrouver comment on charge les love numbers. 
            love = io.loadmat(way)
            self.mode_found=love['mode_found']
            self.k_amp=love['k_amp']
            self.h_amp=love['h_amp']
            self.k_amp_tide=love['k_amp_tide']
            self.h_amp_tide=love['h_amp_tide']
            self.spoles=love['spoles']     
            self.k_el=love['k_el']
            self.k_el_tide=love['k_el_tide']
            self.h = love_lm(love['h_el'],maxdeg)
            self.k = love_lm(love['k_el'],maxdeg)
            self.h_tide = love_lm(love['h_el_tide'],maxdeg)
            self.k_tide = love_lm(love['k_el_tide'],maxdeg)
            self.E = 1 + self.k - self.h
            self.T = sphericalobject(coeff=get_tlm(maxdeg-1,a,Me))
            self.beta_l = np.zeros((self.time_step_number-1,self.time_step_number-1,maxdeg+1))
            self.beta_G_l = np.zeros((self.time_step_number-1,self.time_step_number-1,maxdeg+1))
            self.beta_R_l = np.zeros((self.time_step_number-1,self.time_step_number-1,maxdeg+1))
            self.beta_konly_l = np.zeros((self.time_step_number-1,self.time_step_number-1))
            self.beta_G_konly_l = np.zeros((self.time_step_number-1,self.time_step_number-1))
            self.beta_R_konly_l = np.zeros((self.time_step_number-1,self.time_step_number-1))
            for t_it in range(1,self.time_step_number):  # loop on the time step
                #initialise the temporary variable for calculation use
                beta_l_int=self.beta_l[t_it-1]
                beta_konly_l_int=self.beta_konly_l[t_it-1]
                beta_G_l_int=self.beta_l[t_it-1]
                beta_G_konly_l_int=self.beta_konly_l[t_it-1]
                beta_R_l_int=self.beta_l[t_it-1]
                beta_R_konly_l_int=self.beta_konly_l[t_it-1]
                for n in range(1,t_it):
                    # initialize the beta for each time step
                    beta = np.zeros((maxdeg,)) 
                    beta_G = np.zeros((maxdeg,)) 
                    beta_R = np.zeros((maxdeg,)) 
                    for lm in range(maxdeg):
                        # num_mod is the number of love number needed for each degree-order
                        num_mod = self.mode_found[lm][0] 
                        # calculate the beta coefficient for each time step and degree-order
                        beta[lm] = np.sum((self.k_amp[lm,:num_mod] - self.h_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                        beta_G[lm]=np.sum((self.k_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                        beta_R[lm]=np.sum((self.h_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-time_step[t_it] + time_step[n]))))
                    # add to the beta a coefficient of value 0 !!!!
                    beta_l_int[n-1,:]=np.concatenate((np.array([0]),beta))
                    beta_G_l_int[n-1,:]=np.concatenate((np.array([0]),beta_G))
                    beta_R_l_int[n-1,:]=np.concatenate((np.array([0]),beta_R))
                    # for rotation only needed for degree 2
                    lm=1
                    num_mod=self.mode_found[lm][0]
                    # calculate the beta konly 
                    beta_konly_l_int[n-1] = np.sum((self.k_amp[lm,:num_mod]-self.h_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
                    beta_G_konly_l_int[n-1] = np.sum((self.k_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
                    beta_R_konly_l_int[n-1] = np.sum((self.h_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-time_step[t_it] + time_step[n]))))
                #for each time step update the beta coeffcient
                self.beta_l[t_it-1]=beta_l_int
                self.beta_R_l[t_it-1]=beta_R_l_int
                self.beta_G_l[t_it-1]=beta_G_l_int
            self.beta_konly_l[t_it-1]=beta_konly_l_int
            self.beta_G_konly_l[t_it-1]=beta_G_konly_l_int
            self.beta_R_konly_l[t_it-1]=beta_R_konly_l_int
        elif type is 'spectral':
            self.h=np.loadtxt(way+'/h_ve.dat',unpack=True)[:,:maxdeg]
            self.h=self.h[1,:]+self.h[2,:]*1j
            self.k=np.loadtxt(way+'/k_ve.dat',unpack=True)[:,:maxdeg]
            self.k=self.k[1,:]+self.k[2,:]*1j
            self.h_e=np.loadtxt(way+'/h_e.dat',unpack=True)[1,:maxdeg]
            self.k_e=np.loadtxt(way+'/k_e.dat',unpack=True)[1,:maxdeg]
            self.h_e[0]=0
            self.k_e[0]=0
            # self.h_e=np.abs(self.h)*np.exp(np.imag(self.h)/np.real(self.h))#-np.loadtxt(way+'/k_e.dat',unpack=True)[1,:maxdeg]
            # self.k_e=np.abs(self.k)*np.exp(np.imag(self.k)/np.real(self.k))#-np.loadtxt(way+'/k_e.dat',unpack=True)[1,:maxdeg]
            time_diff=np.repeat(self.time_step[::-1,np.newaxis],self.time_step_number,axis=1)-np.repeat(self.time_step[::-1,np.newaxis],self.time_step_number,axis=1).T
            time_diff[time_diff<0]=0
            # self.beta_R_l=np.repeat(np.repeat(self.h[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)*np.exp(1j*2*np.pi/T*np.repeat(time_diff[:,:,np.newaxis],self.h.shape[0],axis=2))#-np.repeat(np.repeat(self.h_e[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)
            # self.beta_G_l=np.repeat(np.repeat(self.k[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)*np.exp(1j*2*np.pi/T*np.repeat(time_diff[:,:,np.newaxis],self.k.shape[0],axis=2))#-np.repeat(np.repeat(self.k_e[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)
            # self.beta_R_l=np.real(np.repeat(np.repeat(self.h[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)*np.exp(1j*2*np.pi/T*np.repeat(time_diff[:,:,np.newaxis],self.h.shape[0],axis=2)))#-np.repeat(np.repeat(self.h_e[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)
            # self.beta_G_l=np.real(np.repeat(np.repeat(self.k[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)*np.exp(1j*2*np.pi/T*np.repeat(time_diff[:,:,np.newaxis],self.k.shape[0],axis=2)))#-np.repeat(np.repeat(self.k_e[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)

            Phi=np.imag(np.repeat(np.repeat(self.h[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0))/np.real(np.repeat(np.repeat(self.h[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0))


            self.beta_R_l=np.real(np.abs(np.repeat(np.repeat(self.h[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0))*np.exp(1j*2*np.pi/T*np.repeat(time_diff[:,:,np.newaxis],self.h.shape[0],axis=2)+Phi))-np.repeat(np.repeat(self.h_e[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)-np.abs(self.h)*np.exp(np.imag(self.h)/np.real(self.h))

            self.beta_G_l=np.real(np.repeat(np.repeat(self.k[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)*np.exp(1j*2*np.pi/T*np.repeat(time_diff[:,:,np.newaxis],self.k.shape[0],axis=2)+Phi))-np.repeat(np.repeat(self.k_e[np.newaxis,:],self.time_step_number,axis=0)[np.newaxis,:,:],self.time_step_number,axis=0)-np.abs(self.k)*np.exp(np.imag(self.k)/np.real(self.k))

            self.beta_G_l[0,:,0]=0
            self.beta_G_l[0,0,:]=0
            self.beta_l=self.beta_G_l-self.beta_R_l
            calc_beta_counter(self,maxdeg)


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
        if applied is 'beta':
            self.beta_G=0
            self.beta_R=0
            self.beta_l=self.beta_G_l-self.beta_R_l
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
        if applied is 'beta':
            self.beta_G_tide=0
            self.beta_R_tide=0
            self.beta_l_tide=self.beta_G_l_tide-self.beta_R_l_tide
            self._tide=np.array(self.beta_l)[:,:,self.beta_counter.astype(int)[:6]]
        elif applied is 'beta_R':
            self.beta_l_tide=0
            self.beta_G_tide=0
            self.beta_R_tide=np.array(self.beta_R_l_tide)[:,:,self.beta_counter.astype(int)[:6]]
        elif applied is 'beta_G':
            self.beta_l_tide=0
            self.beta_R_tide=0
            self.beta_G_tide=np.array(self.beta_G_l_tide)[:,:,self.beta_counter.astype(int)[:6]]
    
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
