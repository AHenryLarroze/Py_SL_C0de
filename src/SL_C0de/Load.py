from .spharm import sphericalobject
import numpy as np
import sys
import logging
import math

def calc_rot_visc(L,model_p,t_it):
    '''
    .. note::
        This function must be updated to work with the new versions of the code.

    '''
    # extract degree 2 coefficient from the load
    L00 = L.delL.coeff[0]
    L20 = L.delL.coeff[3]
    L21 = L.delL.coeff[4] 
        
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
        V_lm = np.dot(model_p.love.beta_konly_l[t_it-1,:t_it-1],L.sdelI[:t_it-1,:])
        V_lm_T = np.dot(model_p.love.beta_konly_tide[t_it-1,:t_it-1],L.sdelm[:t_it-1,:])
    temp = 1/(1-model_p.love.k_tide.coeff[1]/model_p.k_f)*(1/model_p.CminA * ((1+model_p.love.k.coeff[1])*I + V_lm.squeeze()) + V_lm_T.squeeze()/model_p.k_f)
    # calculate the perturbation to the rotational potential from Milne 1998
    m1=temp[0]
    m2=temp[1]
    temp = -1/(1-model_p.love.k_tide.coeff[1]/model_p.k_f)*(1/C * ((1+model_p.love.k.coeff[1])*I + V_lm.squeeze()))
    m3=temp[2]
    m=np.array([m1,m2,m3])
    # update the rotational load using.
    L.sdelI[t_it-1,:] = I - np.sum(L.sdelI[:t_it-1,:],0)
    L.sdelm[t_it-1,:] = m - np.sum(L.sdelm[:t_it-1,:],0)
    # calculate the rotational perturbation of the earth potential, just for the 6 first coefficient (no use to calculate further)
    L.delLa.coeff = np.zeros(L.delL.coeff.shape)+1j*0
    L.delLa.coeff[0] = model_p.a**2 * model_p.omega**2/3 * (np.sum(m**2) + 2*m3)+1j*0
    L.delLa.coeff[3] = model_p.a**2 * model_p.omega**2/(6*5**.5) * (m1**2 + m2**2 - 2*m3**2 - 4*m3)+1j*0
    L.delLa.coeff[4] = model_p.a**2 * model_p.omega**2/30**.5 * (m1*(1+m3) - 1j*m2*(1+m3))
    L.delLa.coeff[5] = model_p.a**2 * model_p.omega**2/5**.5 * 24**.5 * ( (m2**2-m1**2) + 1j*2*m1*m2 )


class LOAD(object):
    """
    The _`LOAD` class is used to estimate the deformation of the earth and geoid. this method allow you to calculate the elastic response and the vicous one separately. 

    .. note::
        This class might be included directly inside the LOAD class to reduce the amount of code. 

    ...

    Attributes
    ----------
        maxdeg : int
            The maximum spherical harmonic degree of the load used.
        time_step : np.array([time_step_number,])
            The time steps used to modelize the load.
        
    Methods
    -------
        `save_prev`_ : 
            A method used to save the actual computation in order to use it in a future computation after alteration of this computation. 
        `calc_viscuous`_ :
            A method used to compute the viscous deformation of the ground and geoïd.
        
    """
    
    def __init__(self,maxdeg,time_step):
        """
            Parameters
            ----------
        """
        time_step_number=len(time_step)
        N=int((maxdeg+1)*(maxdeg+2)/2)
        self.N=N
        self.delL=sphericalobject(coeff=np.zeros((N,)))
        # self.sdelL=np.zeros((time_step_number,N))+0j
        # self.delL_prev=sphericalobject(coeff=np.zeros((N,)))
        # self.sdelI=np.zeros((time_step_number-1,3))+0j
        # self.sdelm=np.zeros((time_step_number-1,3))+0j
        # self.delLa=sphericalobject(coeff=np.zeros((self.N,)))
        # self.sdelLa=np.zeros((time_step_number,self.N))+0j
        # self.delLa_prev= sphericalobject(coeff=np.zeros((self.N,)))
        self.V_lm=sphericalobject(coeff=np.zeros((N,)))
        self.V_lm_T=sphericalobject(coeff=np.zeros((N,)))


    # def modify(self,t_it,delL):
    #     self.delL.modify(delL,'coeff')
    #     self.sdelL[t_it-1,:]=self.delL.coeff-self.delL_prev.coeff
    #     return
    
    def save_prev(self):
        '''
        The _`save_prev` method is used to save a previously computed value of load. This can be used to calculate deformation based on a difference. 
        
        .. note::
            This method might not be required in the present code.

        Attribute
        ---------
            None
        
        Return
        ------
            None

        '''
        self.delL_prev.modify(self.delL.coeff.copy(),'coeff')
        self.delLa_prev.modify(self.delLa.coeff.copy(),'coeff')
    
    def calc_viscuous(self,sdelL,beta,t_it):
        '''
        The _`calc_viscuous` method is used to calculate the ground and geoïd deformation based on viscuous love numbers.

        Attribute
        ---------
            sdelL : np.array([t_it,(maxdeg+1)(maxdeg+20/2)])
                The load grid used to estimate the ground vertical mouvement. This include all previous loading history because of the viscous comportment of earth. 
            beta : np.array([time_step_number,time_step_number,(maxdeg+1)(maxdeg+2)/2])
                The beta love numbers used to compute the earth deformation to include the viscous part. These love numbers are particularly heavy in the memory due to the representation of the time. 
            t_it : int
                The time iteration at wich the computation is performed.
            
        '''
        if t_it==0 :
            self.V_lm.coeff=beta[0,0]*sdelL
        else :
            self.V_lm.coeff=(beta[t_it-1,:t_it-1]*sdelL).sum(0)

    
    # def calc_rotational_potential(self,model_p,t_it):
    #     calc_rot_visc(self,model_p,t_it)
    #     self.sdelLa[t_it-1]=self.delLa.coeff-self.delLa_prev.coeff
    
    # def calc_viscuous_load_T(self,model_p,t_it,sdelLa):
    #     results=model_p.pool.starmap(par.f_V_lm_T,zip(model_p.love.beta_tide.transpose()[:6],[t_it for i in range(6)],sdelLa.transpose()[:6]))
    #     self.V_lm_T.modify(np.concatenate((results,np.zeros((int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2-6),))+1j)),'coeff')