import numpy as np
from numpy.matlib import repmat #used to add repmat 
from .spharm import sphericalobject

def love_lm(num,maxdeg):
    '''
    Exrtact from love numbers the h_lm spherical coefficient. 

    Parameters : 
        num (np.array): love number coefficient of the size of the entry file
        grid (object): output of the class GRID
        
    See the documentation of the cited class object for more information on different parameters used in the function.

    Returns : 
        h_lm (np.array): array of the love number transformed of size 1, maxdeg. 

    Added fields : 
    '''
    num=np.concatenate((np.zeros((1,1)),num)) 
    h_lm=np.repeat(num,np.arange(1,maxdeg+2,1))
    return h_lm # because the output is of the wrong shape we correct it by transposing the matrix

def get_tlm(maxdeg,a,Me, group='l'):
    '''
    Generate the T spherical harmonic coefficient. Retrouver d'où viennet ces coeff. 

    Parameters : 
        maxdeg (int): maximum degree of spherical harmonic defined in the model parameters. 
        earth (object): output of the class World_Model_Parameter. 
        
    See the documentation of the cited class object for more information on different parameters used in the function.

    Returns : 
        T_lm (np.array): array of size maxdeg + 1 

    Added fields : 
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
    define the degree of the spherical harmonic for each beta coefficient 
    
    Parameters : 
        self (object): output of the class LOVE
        
    See the documentation of the cited class object for more information on different parameters used in the function.
    
    Returns : 
    
    Added fields : 
        beta_counter : np.array ((maxdeg+1)*(maxdeg+2)/2)
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
    A class used to represent the Love numbers parameters

    ...

    Attributes
    ----------
        mode_found : np.array (PREM_maxdeg)
        k_amp : np.array(PREM_maxdeg, maximum(mode_found))
        h_amp : np.array(PREM_maxdeg, maximum(mode_found))
        k_amp_tide : np.array(PREM_maxdeg, maximum(mode_found))
        h_amp_tide : np.array(PREM_maxdeg, maximum(mode_found))
        spoles : np.array(PREM_maxdeg, maximum(mode_found))
        k_el : np.array(PREM_maxdeg)
        k_el_tide : np.array(PREM_maxdeg)
        h : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        k : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        h_tide : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        k_tide : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        E : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        T : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        E_T : np.array (maxdeg+1, maxdeg+1) !!!! Corriger ce truc avec le maxdeg peut être le caller mieux !!!!
        beta_l : 
        beta_k_only :
        beta_tide :
        beta_konly_tide :
        beta_counter :

    """
    def __init__(self,maxdeg,way,time_step,a,Me):
        """
    Parameters
    ----------
    grid : object (from class GRID)
    way : !!!! A corriger !!!! 
        """

        # Load the Load Love Numbers
        self.time_step=time_step
        self.maxdeg=maxdeg
        self.time_step_number=len(time_step)
        self.h_e=np.loadtxt(way+'/h_e.dat',unpack=True)[1,:maxdeg]
        self.k_e=np.loadtxt(way+'/k_e.dat',unpack=True)[1,:maxdeg]
        self.h=np.repeat(self.h_e,np.arange(1,maxdeg+1,1))
        self.k=np.repeat(self.k_e,np.arange(1,maxdeg+1,1))
        self.k_ve=np.loadtxt(way+'/k_ve.dat',unpack=True)
        self.h_ve=np.loadtxt(way+'/h_ve.dat',unpack=True)

        self.love_time=np.loadtxt(way+'/time.dat',unpack=True)


        time_step_diff=-(time_step.reshape(-1,1)-time_step)
        time_step_diff=time_step_diff.flatten()
        time_step_diff=time_step_diff[time_step_diff>0]


        time_interval=np.where((repmat(time_step_diff,len(self.love_time),1).transpose()-repmat(self.love_time,len(time_step_diff),1))==0)[1]
        time_interval_tri=np.zeros((self.time_step_number-1,self.time_step_number-1))
        time_interval_tri[np.triu_indices(self.time_step_number-1,0)]=time_interval
        time_interval_tri=time_interval_tri.transpose().astype(int)

        k_ve=np.concatenate((np.zeros((self.k_ve.shape[0],1)),self.k_ve),axis=1)
        h_ve=np.concatenate((np.zeros((self.h_ve.shape[0],1)),self.h_ve),axis=1)
        k_e=np.concatenate((np.zeros((1,)),self.k_e))
        h_e=np.concatenate((np.zeros((1,)),self.h_e))
        self.beta_G_l=k_ve[time_interval_tri,:maxdeg].squeeze()-k_e[:maxdeg]
        self.beta_G_l[time_interval_tri<=0]=self.beta_G_l[time_interval_tri<=0]*0
        self.beta_R_l=h_ve[time_interval_tri,:maxdeg].squeeze()-h_e[:maxdeg]
        self.beta_R_l[time_interval_tri<=0]=self.beta_R_l[time_interval_tri<=0]*0
        self.beta_l=self.beta_G_l-self.beta_R_l
        self.beta_konly_l=self.k_ve[time_interval_tri,1]-self.k_e[0]-(self.h_ve[time_interval_tri,1]-self.h_e[0])


        # Load the Tide Love Numbers
        self.h_tide_e=np.loadtxt(way+'/h_e_T.dat',unpack=True)[1,:maxdeg]
        self.k_tide_e=np.loadtxt(way+'/k_e_T.dat',unpack=True)[1,:maxdeg]
        self.h_tide=np.repeat(self.h_tide_e,np.arange(1,maxdeg+1,1))
        self.k_tide=np.repeat(self.k_tide_e,np.arange(1,maxdeg+1,1))
        self.k_tide_ve=np.loadtxt(way+'/k_ve_T.dat',unpack=True)[1:-1,:]
        self.h_tide_ve=np.loadtxt(way+'/h_ve_T.dat',unpack=True)[1:-1,:]

        k_tide_ve=np.concatenate((np.zeros((self.k_ve.shape[0],1)),self.k_ve),axis=1)
        h_tide_ve=np.concatenate((np.zeros((self.h_ve.shape[0],1)),self.h_ve),axis=1)
        k_tide_e=np.concatenate((np.zeros((1,)),self.k_e))
        h_tide_e=np.concatenate((np.zeros((1,)),self.h_e))
        self.beta_tide=k_tide_ve[time_interval_tri,:maxdeg+1].squeeze()-k_tide_e[:maxdeg+1]-(h_tide_ve[time_interval_tri,:maxdeg+1].squeeze()-h_tide_e[:maxdeg+1])
        self.beta_konly_tide=self.k_tide_ve[time_interval_tri,1]-self.k_tide_e[0]-(self.h_tide_ve[time_interval_tri,1]-self.h_tide_e[0])
        
        self.E = 1 + self.k - self.h
        self.E_T = 1 + self.k_tide - self.h_tide
        self.T = sphericalobject(coeff=get_tlm(maxdeg-1,a,Me))

        calc_beta_counter(self,maxdeg)
        self.beta_G_l=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)]
        self.beta_R_l=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)]
        self.beta_l=np.array(self.beta_l)[:,:,self.beta_counter.astype(int)]
        self.beta_konly_l=np.array(self.beta_konly_l)
        self.beta_tide=np.array(self.beta_tide)[:,:,self.beta_counter.astype(int)]
        self.beta_konly_tide=np.array(self.beta_konly_tide)
    
    def calc_beta_G(self,model_p):
        self.beta_G_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
        k_e=np.loadtxt(model_p.love_way+'/k_e.dat',unpack=True)[1,:]
        for t_j in range (1,model_p.time_step_number):
            for t_n in range(1,t_j):
                self.beta_G_l[t_j-1,t_n-1,:]=self.k_ve[np.where(np.round((model_p.time_step[t_n-1]-model_p.time_step[t_j-1]),2)==self.love_time)[0],:model_p.maxdeg+1][0].squeeze()-k_e[:model_p.maxdeg+1]
        self.beta_G_l=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)[:-1]-1]

    def calc_beta_R(self,model_p):
        self.beta_R_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
        h_e=np.loadtxt(model_p.love_way+'/h_e.dat',unpack=True)[1,:]
        for t_j in range (1,model_p.time_step_number):
            for t_n in range(1,t_j):
                self.beta_R_l[t_j-1,t_n-1,:]=self.h_ve[np.where(np.round((model_p.time_step[t_n-1]-model_p.time_step[t_j-1]),2)==self.love_time)[0],:model_p.maxdeg+1][0].squeeze()-h_e[:model_p.maxdeg+1]
        self.beta_R_l=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)[:-1]-1]
