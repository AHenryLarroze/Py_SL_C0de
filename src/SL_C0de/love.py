import numpy as np
import numpy.matlib as npmlib #used to add repmat 
from .spharm import sphericalobject

def love_lm(num,grid, group='m'):
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
    h=np.concatenate((np.zeros((1,1)),num)) # a 0 coefficient is added at the beginning of the love numbers !!!! WHY ? !!!!
    h_lm=np.array([]) # initialisation of the h_lm coefficient matrix

    #   for n in range (grid.maxdeg+1): # we work over the maxdeg + 1 degree !!!! I must understand why !!!! 
    #         h_add=npmlib.repmat(h[n],1,n+1) # reproduce the love number coefficient for n+1 time
    #         if h_lm.shape[0]==0 : # if ther is nothing in h_lm
    #             h_lm=h_add #set h_lm as h_add
    #         else : 
    #             h_lm= np.concatenate((h_lm,h_add),1) # add h_add to h_lm 
    results=grid.pool.starmap(par.f_love_lm,zip(h,[n for n in range(grid.maxdeg+1)]))
    h_lm=np.concatenate(results,axis=1)
    return h_lm.transpose() # because the output is of the wrong shape we correct it by transposing the matrix

def get_tlm(maxdeg,earth, group='l'):
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
    const = 4*np.pi*earth.a**3/earth.M_e # setting the tide constant for earth !!!! CHECK !!!!
    for n in range(maxdeg+1) :
            T[n]=const/(2*n+1) # for each nth add to T the earth constant moduladed by the nth step
            T_add=npmlib.repmat(T[n],1,n+1) # prepare T_add as the repetition of T
            if n==0 :
                T_lm=T_add # if there is nothing in T_lm preset T_lm as T_add
            else :
                T_lm=np.concatenate((T_lm,T_add),1) # else add T_add to T_lm
    return np.squeeze(T_lm) # we have to squeeze the array so that the requested indices are directly on the good axis

def calc_beta_counter(self):
    '''
    define the degree of the spherical harmonic for each beta coefficient 
    
    Parameters : 
        self (object): output of the class LOVE
        
    See the documentation of the cited class object for more information on different parameters used in the function.
    
    Returns : 
    
    Added fields : 
        beta_counter : np.array ((maxdeg+1)*(maxdeg+2)/2)
    '''
    self.beta_counter = np.ones((self.h.shape[0]+1,));
    l_it = 1
    for lm_it in range(1,len(self.h)+1):
        # for each time step if the coefficient indices is one degree then , all the next values in beta_counter will have the previous of this one +1. litteraly the degree associated to there indices.
        if lm_it == l_it*(l_it+1)/2:
            self.beta_counter[lm_it] = self.beta_counter[lm_it-1]+1
            l_it = l_it+1
        else:
            self.beta_counter[lm_it] = self.beta_counter[lm_it-1]
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
    def __init__(self,model_p,way):
        """
    Parameters
    ----------
    grid : object (from class GRID)
    way : !!!! A corriger !!!! 
        """

        # Load the Load Love Numbers
        self.h_e=np.loadtxt(way+'/h_e.dat',unpack=True)[1,:]
        self.k_e=np.loadtxt(way+'/k_e.dat',unpack=True)[1,:]
        self.h=love_lm(np.expand_dims(self.h_e,axis=1),model_p)
        self.k=love_lm(np.expand_dims(self.k_e,axis=1),model_p)
        self.k_ve=np.loadtxt(way+'/k_ve.dat',unpack=True)[1:,:]
        self.h_ve=np.loadtxt(way+'/h_ve.dat',unpack=True)[1:,:]

        self.love_time=np.loadtxt(way+'/time.dat',unpack=True)
        #indices=np.linspace(0,len(love_time),1)
        # indices_love=np.zeros((len(model_p.grid.time_step),),dtype=int)
        # for i in range(0,len(model_p.grid.time_step)):
        #     indices_love[i]=np.where(model_p.grid.time_step[i]==love_time)[0]
        

        # self.k_ve=self.k_ve[indices_love,:]
        # self.h_ve=self.h_ve[indices_love,:]

        self.beta_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
        self.beta_konly_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1))


        for t_j in range (1,len(model_p.grid.time_step)):
            for t_n in range(1,t_j):
                dt_n=np.where(np.round((model_p.time_step[t_n-1]-model_p.time_step[t_j-1]),2)==self.love_time)[0]
                self.beta_l[t_j-1,t_n-1,1:]=self.k_ve[dt_n,:model_p.maxdeg][0].squeeze()-self.k_e[:model_p.maxdeg]-(self.h_ve[dt_n,:model_p.maxdeg][0].squeeze()-self.h_e[:model_p.maxdeg])
                self.beta_konly_l[t_j-1,t_n-1]=self.k_ve[dt_n,1]-self.k_e[0]-(self.h_ve[dt_n,1]-self.h_e[0])

        # Load the Tide Love Numbers
        self.h_tide_e=np.loadtxt(way+'/h_e_T.dat',unpack=True)[1,:]
        self.k_tide_e=np.loadtxt(way+'/k_e_T.dat',unpack=True)[1,:]
        self.h_tide=love_lm(np.expand_dims(self.h_tide_e,axis=1),model_p)
        self.k_tide=love_lm(np.expand_dims(self.k_tide_e,axis=1),model_p)
        self.k_tide_ve=np.loadtxt(way+'/k_ve_T.dat',unpack=True)[1:-1,:]
        self.h_tide_ve=np.loadtxt(way+'/h_ve_T.dat',unpack=True)[1:-1,:]

        self.beta_tide = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
        self.beta_konly_tide = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1))

        for t_j in range (1,len(model_p.grid.time_step)):
            for t_n in range(1,t_j):
                dt_n=np.where(np.round((model_p.time_step[t_n-1]-model_p.time_step[t_j-1]),2)==self.love_time)[0]
                self.beta_tide[t_j-1,t_n-1,1:]=self.k_tide_ve[dt_n,:model_p.maxdeg][0].squeeze()-self.k_tide_e[:model_p.maxdeg]-(self.h_tide_ve[dt_n,:model_p.maxdeg][0].squeeze()-self.h_tide_e[:model_p.maxdeg])
                self.beta_konly_tide[t_j-1,t_n-1]=self.k_tide_ve[dt_n,1]-self.k_tide_e[0]-(self.h_tide_ve[dt_n,1]-self.h_tide_e[0])
        
        self.E = 1 + self.k - self.h
        self.E_T = 1 + self.k_tide - self.h_tide
        self.T = sphericalobject(get_tlm(model_p.maxdeg,model_p),'coeff')

        calc_beta_counter(self)
        self.beta_l=np.array(self.beta_l)[:,:,self.beta_counter.astype(int)[:-1]-1]
        self.beta_konly_l=np.array(self.beta_konly_l)
        self.beta_tide=np.array(self.beta_tide)[:,:,self.beta_counter.astype(int)[:-1]-1]
        self.beta_konly_tide=np.array(self.beta_konly_tide)
        
        # initiate mapping from l to lm
        self.initiale_spherical_object()

    def initiale_spherical_object(self):
        self.h = sphericalobject(np.squeeze(self.h),'coeff')
        self.k = sphericalobject(np.squeeze(self.k),'coeff')
        self.h_tide = sphericalobject(np.squeeze(self.h_tide),'coeff')
        self.k_tide = sphericalobject(np.squeeze(self.k_tide),'coeff')
        self.E = sphericalobject(np.squeeze(self.E),'coeff')
        self.E_T = sphericalobject(np.squeeze(self.E_T),'coeff')
    
    def calc_beta_G(self,model_p):
        self.beta_G_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
        k_e=np.loadtxt(model_p.love_way+'/k_e.dat',unpack=True)[1,:]
        for t_j in range (1,model_p.time_step_number):
            for t_n in range(1,t_j):
                #print(np.where(np.round(model_p.time_step[t_n]-model_p.time_step[t_j],2)))
                #print(t_j,t_n,np.where(np.round((model_p.time_step[t_n]-model_p.time_step[t_j]),2)==self.love_time)[0])
                self.beta_G_l[t_j-1,t_n-1,:]=self.k_ve[np.where(np.round((model_p.time_step[t_n-1]-model_p.time_step[t_j-1]),2)==self.love_time)[0],:model_p.maxdeg+1][0].squeeze()-k_e[:model_p.maxdeg+1]
        self.beta_G_l=np.array(self.beta_G_l)[:,:,self.beta_counter.astype(int)[:-1]-1]

    def calc_beta_R(self,model_p):
        self.beta_R_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
        h_e=np.loadtxt(model_p.love_way+'/h_e.dat',unpack=True)[1,:]
        for t_j in range (1,model_p.time_step_number):
            for t_n in range(1,t_j):
                #print(np.where(np.round(model_p.time_step[t_n]-model_p.time_step[t_j],2)))
                #print(t_j,t_n,np.where(np.round((model_p.time_step[t_n]-model_p.time_step[t_j]),2)==self.love_time)[0])
                self.beta_R_l[t_j-1,t_n-1,:]=self.h_ve[np.where(np.round((model_p.time_step[t_n-1]-model_p.time_step[t_j-1]),2)==self.love_time)[0],:model_p.maxdeg+1][0].squeeze()-h_e[:model_p.maxdeg+1]
        self.beta_R_l=np.array(self.beta_R_l)[:,:,self.beta_counter.astype(int)[:-1]-1]
