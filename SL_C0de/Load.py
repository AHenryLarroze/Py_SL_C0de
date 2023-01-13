from spharm import sphericalobject
import numpy as np
import sys
import logging
import par_Load as par
import math

def calc_rot_visc(L,model_p,t_it):
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
    A class used to represent the different load over the world

    ...

    Attributes
    ----------

        
    Methods
    -------
        create_GRID() 
    """
    
    def __init__(self,model_p):
        """
            Parameters
            ----------
        """
        N=int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2)
        self.N=N
        self.delL=sphericalobject(np.zeros((N,)),'coeff')
        self.sdelL=np.zeros((model_p.time_step_number,N))+0j
        self.delL_prev=sphericalobject(np.zeros((N,)),'coeff')
        self.sdelI=np.zeros((len(model_p.grid.time_step)-1,3))+0j
        self.sdelm=np.zeros((len(model_p.grid.time_step)-1,3))+0j
        self.delLa=sphericalobject(np.zeros((self.N,)),'coeff')
        self.sdelLa=np.zeros((model_p.time_step_number,self.N))+0j
        self.delLa_prev= sphericalobject(np.zeros((self.N,)),'coeff')
        self.V_lm=sphericalobject(np.zeros((N,)),'coeff')
        self.V_lm_T=sphericalobject(np.zeros((N,)),'coeff')


    def modify(self,t_it,delL):
        self.delL.modify(delL,'coeff')
        self.sdelL[t_it-1,:]=self.delL.coeff-self.delL_prev.coeff
        return
    
    def save_prev(self):
        self.delL_prev.modify(self.delL.coeff.copy(),'coeff')
        self.delLa_prev.modify(self.delLa.coeff.copy(),'coeff')
    
    def calc_viscuous_load(self,model_p,t_it,sdelL,beta):
        # for name in self.load_name:
        #     results=model_p.pool.starmap(par.f_V_lm,zip(model_p.love.beta_R_l.transpose(),[t_it for i in range(int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2))],np.array([self.sdelL[i][name].coeff for i in range(t_it)]).transpose()))
        #     self.V_lm_R[name].modify(np.concatenate(results),'coeff')
        #     results=model_p.pool.starmap(par.f_V_lm,zip(model_p.love.beta_G_l.transpose(),[t_it for i in range(model_p.time_step_number)],np.array([self.sdelL[i][name].coeff for i in range(t_it)]).transpose()))
        #     self.V_lm_G[name].modify(np.concatenate(results),'coeff')
        results=model_p.pool.starmap(par.f_V_lm,zip(beta.transpose(),[t_it for i in range(int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2))],sdelL.transpose()))
        self.V_lm.modify(results,'coeff')
    
    def calc_rotational_potential(self,model_p,t_it):
        # delLa,sdelI,sdelm=calc_rot_visc(self,model_p,t_it,name)
        calc_rot_visc(self,model_p,t_it)
        # self.delLa[name].modify(delLa,'coeff')
        # self.sdelI[name].modify(sdelI,'coeff')
        # self.sdelm[name].modify(sdelm,'coeff')
        self.sdelLa[t_it-1]=self.delLa.coeff-self.delLa_prev.coeff
    
    def calc_viscuous_load_T(self,model_p,t_it,sdelLa):
        # for name in self.load_name:
        #     results=model_p.pool.starmap(par.f_V_lm_T,zip(model_p.love.beta_R_tide.transpose(),[t_it for i in range(int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2))],np.array([self.sdelLa[name].coeff for i in range(t_it)]).transpose()))
        #     self.V_lm_R_T[name].modify(np.concatenate(results),'coeff')
        #     results=model_p.pool.starmap(par.f_V_lm_T,zip(model_p.love.beta_G_tide.transpose(),[t_it for i in range(model_p.time_step_number)],np.array([self.sdelLa[name].coeff for i in range(t_it)]).transpose()))
        #     self.V_lm_G_T[name].modify(np.concatenate(results),'coeff')
        results=model_p.pool.starmap(par.f_V_lm_T,zip(model_p.love.beta_tide.transpose()[:6],[t_it for i in range(6)],sdelLa.transpose()[:6]))
        self.V_lm_T.modify(np.concatenate((results,np.zeros((int((model_p.maxdeg+1)*(model_p.maxdeg+2)/2-6),))+1j)),'coeff')