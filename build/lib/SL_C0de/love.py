import numpy as np
import numpy.matlib as npmlib #used to add repmat 
from spharm import sphericalobject
from scipy import io
import math
import par_love as par
from multiprocessing import Pool
import itertools



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

def calc_beta(self,model_p) :
    '''
    Add to the object from the class LOVE the beta 
    coefficient used for the calculation of viscoelastic deformation 
    of earth. Two parameters are added : beta_tide and beta_konly_tide.
    
        Parameters : 
            self (object): output of the class LOVE
            grid (object): an output of the class GRID
            
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
        
        Added fields : 
            beta_l (np.array) : a matrix of the beta coefficient with shape time_step, time_step, maxdeg+2. 
                                Where time_step is the number of time step defined in the GRID object and maxdeg is the maximum 
                                degree of spherical harmonic in the GRID object.
            beta_konly_l (np.array) : a matrix of the beta with only k coefficient with shape time_step, time_step, maxdeg+2. 
                                Where time_step is the number of time step defined in the GRID object and maxdeg is the maximum 
                                degree of spherical harmonic in the GRID object.
    '''
    # i'll must reduce the number of variable used 
    # initialize the output
    self.beta_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1)) 
    self.beta_konly_l = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1))
    # for t_it in range(1,len(model_p.grid.time_step)):  # loop on the time step
    #     #initialise the temporary variable for calculation use
    #     beta_l_int=self.beta_l[t_it-1]
    #     beta_konly_l_int=self.beta_konly_l[t_it-1]
    #     for n in range(1,t_it):
        #     # initialize the beta for each time step
        #     beta = np.zeros((maxdeg,)) 
        #     for lm in range(maxdeg):
        #         # num_mod is the number of love number needed for each degree-order
        #         num_mod = self.mode_found[lm][0] 
        #         # calculate the beta coefficient for each time step and degree-order
        #         beta[lm] = np.sum((self.k_amp[lm,:num_mod] - self.h_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod]* (-model_p.grid.time_step[t_it] + model_p.grid.time_step[n]))))
        #     # add to the beta a coefficient of value 0 !!!!
        #     beta_l_int[n-1,:]=np.concatenate((np.array([0]),beta))
        #     # for rotation only needed for degree 2
        #     lm=1
        #     num_mod=self.mode_found[lm][0]
        #     # calculate the beta konly 
        #     beta_konly_l_int[n-1] = np.sum((self.k_amp[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-model_p.grid.time_step[t_it] + model_p.grid.time_step[n]))))
        # #for each time step update the beta coeffcient
        # self.beta_l[t_it-1]=beta_l_int
        # self.beta_konly_l[t_it-1]=beta_konly_l_int
    results=model_p.pool.starmap(par.f_calc_beta,zip([t_it for t_it in range(1,len(model_p.grid.time_step))],itertools.repeat(model_p.maxdeg),self.beta_l,self.beta_konly_l,itertools.repeat(self.mode_found),itertools.repeat(self.k_amp),itertools.repeat(self.h_amp),itertools.repeat(self.spoles),itertools.repeat(model_p.grid.time_step)))
    self.beta_l=np.concatenate([np.expand_dims(results[i][0],axis=0) for i in range(len(model_p.grid.time_step)-1)],axis=0)
    self.beta_konly_l=np.concatenate([np.expand_dims(results[i][1],axis=0) for i in range(len(model_p.grid.time_step)-1)],axis=0)
    return

def calc_beta_tide(self,model_p):
    '''
    Add to the object from the class LOVE the beta tide 
    coefficient used for the calculation of viscoelastic diformation 
    of earth. Two parameters are added : beta_tide and beta_konly_tide.
    
        Parameters : 
            self (object): output of the class LOVE
            grid (object): an output of the class GRID
            
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
        
        Added fields : 
            beta_tide (np.array) : a matrix of the beta tide coefficient with shape time_step, time_step, maxdeg+2. 
                                Where time_step is the number of time step defined in the GRID object and maxdeg is the maximum 
                                degree of spherical harmonic in the GRID object.
            beta_konly_tide (np.array) : a matrix of the beta tide with only k coefficient with shape time_step, time_step, maxdeg+2. 
                                Where time_step is the number of time step defined in the GRID object and maxdeg is the maximum 
                                degree of spherical harmonic in the GRID object.
    '''
    # initialize the beta coefficient of the love number
    self.beta_tide = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1,model_p.maxdeg+1))
    self.beta_konly_tide = np.zeros((len(model_p.grid.time_step)-1,len(model_p.grid.time_step)-1))
    # for t_it in range(1,len(model_p.grid.time_step)):
    #     beta_tide_int=self.beta_tide[t_it-1]
    #     beta_konly_tide_int=self.beta_konly_tide[t_it-1]
    #     for n in range(1,t_it):
    #         #initialize the beta coefficient at each time step time step for computational use
    #         beta = np.zeros((maxdeg,))
    #         for lm in range(maxdeg):
    #             num_mod = self.mode_found[lm][0]
    #             beta[lm]= sum((self.k_amp_tide[lm,:num_mod] - self.h_amp_tide[lm,:num_mod])/self.spoles[lm,:num_mod]* (1 - np.exp(- self.spoles[lm,:num_mod] * (-model_p.grid.time_step[t_it] + model_p.grid.time_step[n]))))
    #         beta_tide_int[n-1,:]=np.concatenate((np.array([0]),beta))
    #         lm=1
    #         num_mod=self.mode_found[lm][0]
    #         beta_konly_tide_int[n-1] = np.sum((self.k_amp_tide[lm,:num_mod])/self.spoles[lm,:num_mod] * (1 - np.exp(- self.spoles[lm,:num_mod] * (-model_p.grid.time_step[t_it] + model_p.grid.time_step[n]))))
    #     self.beta_tide[t_it-1]=beta_tide_int
    #     self.beta_konly_tide[t_it-1]=beta_konly_tide_int
    #     # for rotation only needed for degree 2
    results=model_p.pool.starmap(par.f_calc_beta_tide,zip([t_it for t_it in range(1,len(model_p.grid.time_step))],itertools.repeat(model_p.maxdeg),self.beta_tide,self.beta_konly_tide,itertools.repeat(self.mode_found),itertools.repeat(self.k_amp_tide),itertools.repeat(self.h_amp_tide),itertools.repeat(self.spoles),itertools.repeat(model_p.grid.time_step)))
    self.beta_tide=np.concatenate([np.expand_dims(results[i][0],axis=0) for i in range(len(model_p.grid.time_step)-1)],axis=0)
    self.beta_konly_tide=np.concatenate([np.expand_dims(results[i][1],axis=0) for i in range(len(model_p.grid.time_step)-1)],axis=0)
    return

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

def calc_rot_visc(model_p,t_it):
    '''
    Update the object from the class spherical_sea_level the fields delLa, sdelI and sdelm.
    delLa are the spherical harmonic coefficient associated to the earth rotation.
    sdelI is the deformation of ice effect due to the rotation.
    sdelm is the deformation associated to the rotation. 
    
        Parameters : 
            SL (object): output of the class spherical_sea_level
            love (object): output of the class LOVE
            t_it (int): time step iteration
            model_p (object): output of the class World_Model_Parameter
            
        See the documentation of the cited class object for more information on different parameters used in the function.
        
        Returns : 
        
        Added fields : 
    '''
    # extract degree 2 coefficient from the load
    L00 = model_p.SL.delL[0]
    L20 = model_p.SL.delL[3]
    L21 = model_p.SL.delL[4] 
        
    C = model_p.C 
    # calculate the load effect constant 
    I1=math.sqrt(32/15)*math.pi*model_p.a**4*np.real(L21)
    I2=math.sqrt(32/15)*math.pi*model_p.a**4*np.imag(L21)
    I3=8/3*math.pi*model_p.a**4*(L00-L20/math.sqrt(5))
    I=np.array([I1,I2,I3])
    if t_it==1 : #initialise the rotational coefficient
            V_lm=np.array([0,0,0])
            V_lm_T=np.array([0,0,0])
    else : #apply the visco elastic properties of the earth on the rotation using the load.
        V_lm = np.dot(model_p.love.beta_konly_l[t_it-1,:t_it-1],model_p.SL.sdelI[:t_it-1,:])
        V_lm_T = np.dot(model_p.love.beta_konly_tide[t_it-1,:t_it-1],model_p.SL.sdelm[:t_it-1,:])
    temp = 1/(1-model_p.love.k_el_tide[1]/model_p.k_f)*(1/model_p.CminA * ((1+model_p.love.k_el[1])*I + V_lm.squeeze()) + V_lm_T.squeeze()/model_p.k_f)
    # calculate the perturbation to the rotational potential from Milne 1998
    m1=temp[0]
    m2=temp[1]
    temp = -1/(1-model_p.love.k_el_tide[1]/model_p.k_f)*(1/C * ((1+model_p.love.k_el[1])*I + V_lm.squeeze()))
    m3=temp[2]
    m=np.array([m1,m2,m3])
    # update the rotational load using.
    model_p.SL.sdelI[t_it-1,:] = I - np.sum(model_p.SL.sdelI[:t_it-1,:],0)
    model_p.SL.sdelm[t_it-1,:] = m - np.sum(model_p.SL.sdelm[:t_it-1,:],0)
    # calculate the rotational perturbation of the earth potential, just for the 6 first coefficient (no use to calculate further)
    model_p.SL.delLa = np.zeros(model_p.SL.delL.shape)+1j*0
    model_p.SL.delLa[0] = model_p.a**2 * model_p.omega**2/3 * (np.sum(m**2) + 2*m3)+1j*0
    model_p.SL.delLa[3] = model_p.a**2 * model_p.omega**2/(6*5**.5) * (m1**2 + m2**2 - 2*m3**2 - 4*m3)+1j*0
    model_p.SL.delLa[4] = model_p.a**2 * model_p.omega**2/30**.5 * (m1*(1+m3) - 1j*m2*(1+m3))
    model_p.SL.delLa[5] = model_p.a**2 * model_p.omega**2/5**.5 * 24**.5 * ( (m2**2-m1**2) + 1j*2*m1*m2 )

class LOVE(object):
    """
    A class used to represent the Love numbers parameters

    ...

    Attributes
    ----------
        !!!! PREM_maxdeg insert a check to block the resolution of the model under the maximum precision of the PREM model !!!!
        mode_found : np.array (PREM_maxdeg)
            !!!! c'est quoi ? !!!! 
        k_amp : np.array(PREM_maxdeg, maximum(mode_found))
            !!!! c'est quoi ? !!!!
        h_amp : np.array(PREM_maxdeg, maximum(mode_found))
            !!!! c'est quoi ? !!!!
        k_amp_tide : np.array(PREM_maxdeg, maximum(mode_found))
            !!!! c'est quoi ? !!!!
        h_amp_tide : np.array(PREM_maxdeg, maximum(mode_found))
            !!!! c'est quoi !!!!
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
        !!!! vérifier le contenu, ça devient compliquer avec ce bazar 
        beta_l : 
        beta_k_only :
        beta_tide :
        beta_konly_tide :
        beta_counter :
        !!!!
        
    Methods
    -------
        pyspharm2shtools()
            A method to convert the spherical harmonic coefficient to the shape used in pyshtools.   

    """
    def __init__(self,grid):
        """
    Parameters
    ----------
    grid : object (from class GRID)
    way : !!!! A corriger !!!! 
        """
        love = io.loadmat('prem.l96C.umVM5.lmVM5.mat')
        self.mode_found=love['mode_found']
        self.k_amp=love['k_amp']
        self.h_amp=love['h_amp']
        self.k_amp_tide=love['k_amp_tide']
        self.h_amp_tide=love['h_amp_tide']
        self.spoles=love['spoles']     
        self.k_el=love['k_el']
        self.k_el_tide=love['k_el_tide']
        self.h = love_lm(love['h_el'],grid)
        self.k = love_lm(love['k_el'],grid)
        self.h_tide = love_lm(love['h_el_tide'],grid)
        self.k_tide = love_lm(love['k_el_tide'],grid)
        self.E = 1 + self.k - self.h
        self.T = get_tlm(grid.maxdeg,grid)
        #self.T_ml = reorder_l_to_m_primary(self.T_lm)
        self.E_T = 1 + self.k_tide - self.h_tide
        
        # calculate betas
        calc_beta(self,grid)
        
        # calculate tidal betas
        calc_beta_tide(self,grid)
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
        self.T = sphericalobject(np.squeeze(self.T),'coeff')
        self.E_T = sphericalobject(np.squeeze(self.E_T),'coeff')
    