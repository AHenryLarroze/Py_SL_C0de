import numpy as np
import love as Love
from ocean_function import spherical_ocean_function
from grid import GRID
from spharm import legendre, GaussQuad
from plot import world_plot
from multiprocessing import Pool
import par_SL_model as par
from spharm import sphericalobject
from Load import LOAD
import multiprocessing as mp
from time import time
import pickle
from datetime import datetime

class SL_model(world_plot):
    """
    A class used to represent the world parameters.

    ...

    Attributes
    ----------
        g : float
            gravitational constant (m/s**2)
        maxdeg : int
            Maximum degree of the spherical harmonic coefficient. (define the resolution of the model)
        M_e : float
            Earth mass (kg)
        omega : float
            Eath's angular velocity (rad/s) !!!! Donner une publication pour justifier cette valeur !!!!
        k_hydro : float 
            Hydrostatic love number of degree 2 
        G : float
            Universal gravitational constant (m**3/(kg*s**2))
        k_f : float
            !!!! trouver la définition de k_f !!!!
        a : int or float (depends on the precision needed)
            Earth Radius (km)
        CminA : float
            smallest and largest principal moments of the inertia tensor
        C : float
            !!!! trouver la définition de C !!!!
        sqrt_32_15 : float
            square root of 32/15 for computationnal velocity
        k_max : int
            maximum number of itteration for solving the sea level equation at each time step
        epsilon : float
            minimal precision for the solution of the sea level equation
        topo_it_max : int
            maximum number of iteration for the convergence of topography 
        !!!! Il devrait y avoir le maximum le critère de convergence pour la topo aussi où est t'il ? !!!!
        
    Methods
    -------
        create_GRID() 
    """
    
    def __init__(self,g=9.80616,maxdeg=64,M_e=5.9742*10**24,omega=7.292*10**(-5),k_hydro=0.934,G=6.67408*10**(-11),k_f=0.942,a = 6371000,C=8.034*10**37,k_max=10,epsilon=10**-4,topo_it_max=10,nb_workers=6,ice_way='ice6g_data',topo_way='topo_SL',sed_way='Irrawady_sedimentationgrid',love_way='VM5.l40.um316.lm5.1024'):
        """
            Parameters
            ----------
        """
        mp.set_start_method('spawn')
        self.rotation=True
        self.g = g # gravitational constant of earth m**2/s
        self.maxdeg=maxdeg # maximum degree of the spherical harmonic used
        self.M_e = M_e # Earth mass [kg]
        self.omega = omega # Earth's angular velocity rad/s
        self.k_hydro = k_hydro # degree-2 hydrostatic Love number
        self.G=G # Newton gravitational constant
        self.k_f=k_f # !!!! A trouver !!!!
        self.a = a # earth radius (m)
        self.CminA = (self.k_f*self.a**5*self.omega**2)/(3*self.G) # a, c = smallest and largest principal moments of the inertia tensor
        self.C=C
        self.k_max = k_max # maximum number of convergence iteration for the resolution of the sea level equation at 1 time
        self.epsilon = epsilon # convergence criterion for the convergence iteration of the sea level equation at 1 time
        self.topo_it_max = topo_it_max # maximum numbre of iteration for the convergence of the initial topography
        self.nb_workers=nb_workers
        self.pool = Pool(self.nb_workers)
        self.ice_way=ice_way
        self.sed_way=sed_way
        self.love_way=love_way
        self.topo_way=topo_way
        #self.init_grided_parameters()

    def set_from_file(self,way):
        with open(way+'.pkl', 'rb') as handle:
            data = pickle.load(handle)
        self.g=data['g']
        self.M_e=data['M_e']
        self.omega=data['omega']
        self.k_hydro=data['k_hydro']
        self.G=data['G']
        self.k_f=data['k_f']
        self.a=data['a']
        self.CminA=data['CminA']
        self.C=data['C']
        self.maxdeg=data['maxdeg']
        self.rotation=data['rotation']
        self.k_max=data['k_max']
        self.epsilon=data['epsilon']
        self.topo_it_max=data['topo_it_max']
        self.nb_workers=data['nb_workers']
        self.time_step=data['time_step']
        self.ice_way=data['ice_way']
        self.sed_way=data['sed_way']
        self.love_way=data['love_way']
        self.topo_way=data['topo_way']
        self.init_grided_parameters()
        self.sed.rho=data['rho_sed']
        self.oc.rho=data['rho sea']
        self.ice.rho=data['rho_ice']
        self.grid.time_step=self.time_step
        self.ice.ice_corrected=data['ice_corrected']
        self.SL.delS.saved=data['delS']
        self.SL.sdelS.saved=data['sdelS']
        self.topo.topo=data['topo']


    def create_GRID(self):
        self.grid=GRID(self)

    def init_grided_parameters(self): 

        print('Grid creation')
        self.P_lm=legendre(self.maxdeg,self.pool) #Calculate the Legendre associated functions for the Gaussian grid
        self.x, self.w= GaussQuad(self.maxdeg) # calculate the Gaussian grid parameters
        self.grid=GRID(self)

        print('Ice loading')
        # create a better object for the input and prepare the input for the functions
        self.ice=self.grid.create_ICE() # adapt the load method of the spherical_ice class to the data you are using
        self.ice.quick_load(self.grid,self.pool,self.ice_way) # change to quick_load for parallelised version ! ne marche pas. 
        self.time_step_number=len(self.grid.time_step)-1

        print('Sed loading')
        self.sed = self.grid.create_SED()
        self.sed.load(self,self.sed_way)

        print('Sea level initialisation')
        self.SL=self.grid.create_SL()

        print('Love number load')
        self.love=self.grid.create_LOVE(self,self.love_way)

        print('Topography load')
        self.topo = self.grid.create_TOPO() # adapt the load method of the spherical_topo class to the data you are using
        self.topo.load(self,self.topo_way)

        print('Ocean initialisation')
        self.oc=spherical_ocean_function().evaluate_ocean(self.topo.topo_pres).grdtocoeff(self)
        self.oc.update_oc_0()# initialize the grid at t=t_it(t_it=0)
        self.oc.area=self.oc.coeff[0]

        print('Initialize Load')
        self.load=LOAD(self)

        

    def run(self):
        for self.topo_it in range(self.topo_it_max): # Il y a un pb sur le chargement de topo à corriger !
            self.topo_loop()
            self.topo.topo_initial=self.ice.topo_pres_ice_corrected - (self.topo.topo[-1,:,:]-self.topo.topo[0,:,:]) 
        return

    def topo_loop(self):
        #regarder si il n'y a pas une boucle while plutot avec un critère de convergence pour stopper la boucle.
        # initialization of the prev variable
        self.reset() #reset need to be before the loop otherwise we would loose the data. FOr future release maybe put this at the end when we would have a good way to save data.

        self.ice.deli_00_prev = 0 #modifier cette entrès quand cette variable seras dans l'objet
        
        #sdelI et sdelm doivent être défini autrement car ils ont une forme particulière !!!
        
        #set topography
        self.topo.topo[0,:,:]=self.topo.topo_initial.copy()
        #for i in range(len(grid.time_step)):
        #    topo.topo[i,:,:]=topo.topo[i,:,:]-ice.ice_corrected[i,:,:]+ice.ice[i,:,:]         
        self.topo.ice_correction(self)
        self.topo.update_topo_0()
        
        
        self.oc.evaluate_ocean(self.topo.topo_0).grdtocoeff(self)
        self.oc.update_oc_0()
        self.oc.save_prev()
        self.SL.saved=np.array([])
        self.oc.saved=np.array([])

        self.reset_save_data()
        
        for t_it in range(1,len(self.grid.time_step)):
            t1=time()
            self.time_loop(t_it)
            t2=time()
            print('topoiteration : ',self.topo_it,' at time : ',self.grid.time_step[t_it],' in ',t2-t1)
        self.ice.topo_pres_ice_corrected = self.topo.topo_pres - self.ice.ice[-1,:,:] + self.ice.ice_corrected[-1,:,:]
        

    def time_loop(self,t_it):
        # initial topo with the ice correction : 
        self.topo.modify(self.topo.topo[t_it,:,:])
        # grd correspond donc au topo_j défini dans le code de kendal et al.
        self.oc.evaluate_ocean(self.topo.grd).grdtocoeff(self)
        self.SL.TO.modify(self.topo.topo_0*(self.oc.grd-self.oc.oc_0_grd)).grdtocoeff(self)
        self.sed.calc_del_sed(t_it).grdtocoeff(self)
        self.ice.calc_del_ice(t_it).grdtocoeff(self)
        
        self.ice.sdeli_00=self.ice.coeff[0]-self.ice.deli_00_prev # peut être à passer dans l'initialisation ou le modify de ice object après leurs création
        
        #apply the loop solving the see level equation
        k = 0
        chi = self.epsilon * 2
        while (k < self.k_max) and (chi >= self.epsilon):
            chi=self.SL_conv_loop(t_it,k) 
            k+=1
        #update the prev variable for the next loop
        self.save_data()
        self.save_prev()
        self.SL.ESL[t_it]=np.real(self.ice.coeff[0]/self.oc.area*self.ice.rho/self.oc.rho)

        #SL.ESL[t_it]=ice.deli_00_prev/oc.area*ice.rho/oc.rho
        self.SL.delSL.modify(self.SL.delSLcurl.grd+self.SL.delPhi_g)
        self.topo.topo[t_it,:,:]=-self.SL.delSL.grd+self.topo.topo_0#update cette ligne avec le nouvel objet spherical_top

    def SL_conv_loop(self,t_it,k):
        if k == 0 and self.topo_it==0:
            #print('oc',oc.coeff.shape,'TO',SL.TO.coeff.shape,'TO_prev',SL.TO.prev.shape,'ICE',ice.coeff.shape,'SED',sed.coeff.shape)
            self.SL.sdelS.modify(self.oc.prev/self.oc.prev[0]*(-self.ice.rho/self.oc.rho*self.ice.sdeli_00 + self.SL.TO.coeff[0]-self.SL.TO.prev[0])-self.SL.TO.coeff-self.SL.TO.prev,'coeff')
        self.SL.delS.modify(self.SL.delS.prev + self.SL.sdelS.coeff,'coeff')
        self.load.modify(t_it,self.oc.rho*self.SL.delS.coeff+self.ice.rho*self.ice.coeff+self.sed.rho*self.sed.coeff)
        # calculate viscous contribution
        # beta contains the viscous love numbers for time t_it,
        # row index goes over the time increments, column
        # index goes over lm
        self.load.calc_viscuous_load(self,t_it,self.load.sdelL,self.love.beta_l)
        if self.rotation : 
            self.load.calc_rotational_potential(self,t_it)
        self.load.calc_viscuous_load_T(self,t_it,self.load.sdelLa)
        # calculate sea level perturbation
        # add ice and sea level and multiply with love numbers
        self.SL.delSLcurl_fl.modify(self.love.E.coeff * self.love.T.coeff *self.load.delL.coeff+self.love.T.coeff*self.load.V_lm.coeff+1/self.g*self.love.E_T.coeff*self.load.delLa.coeff+1/self.g*self.load.V_lm_T.coeff,'coeff')
        self.SL.delSLcurl.modify(self.SL.delSLcurl_fl.coeff - self.ice.coeff+ -self.sed.coeff,'coeff').coefftogrd(self)
        self.SL.RO.modify(self.SL.delSLcurl.grd*self.oc.grd).grdtocoeff(self)
        self.SL.delPhi_g = np.real(1/self.oc.coeff[0] * (- self.ice.rho/self.oc.rho*self.ice.coeff[0] - self.SL.RO.coeff[0] + self.SL.TO.coeff[0]))
        self.SL.sdelS_new=self.SL.RO.coeff + self.SL.delPhi_g*self.oc.coeff -  self.SL.TO.coeff - self.SL.delS.prev
        chi = np.abs((np.sum(np.abs(self.SL.sdelS_new)) - np.sum(np.abs(self.SL.delS.coeff))) / np.sum(np.abs(self.SL.delS.coeff)))
        self.SL.sdelS.modify(self.SL.sdelS_new.copy(),'coeff')
        return chi
    
    def reset(self):
        N=int((self.maxdeg+1)*(self.maxdeg+2)/2)
        n=self.SL.TO.coeff.shape[0]
        self.SL.TO.prev=np.zeros(n)+1j*0
        self.load.delLa_prev=sphericalobject(np.zeros((N,)),'coeff')
        self.load.delL_prev=sphericalobject(np.zeros((N,)),'coeff')
        self.SL.delS.prev=np.zeros(n)+1j*0
        self.load.sdelI=np.zeros((len(self.grid.time_step)-1,3))+1j
        self.load.sdelm=np.zeros((len(self.grid.time_step)-1,3))+1j
        self.load.sdelLa=np.zeros((self.time_step_number,N))+1j
        self.load.sdelL=np.zeros((self.time_step_number,N))+1j
        self.SL.ESL = np.zeros((len(self.grid.time_step),))

    def calc_RSL(self):
        self.SL.RSL=np.zeros(self.topo.topo.shape)
        for t_it in range(len(self.grid.time_step)):
            self.SL.RSL[t_it,:,:]=self.topo.topo[-1,:,:]-self.ice.ice_corrected[-1,:,:]+self.sed.sed[-1,:,:]-(self.topo.topo[t_it,:,:]-self.ice.ice_corrected[t_it,:,:]+self.sed.sed[t_it,:,:])


    def save_prev(self):
        self.SL.delS.save_prev()
        self.SL.TO.save_prev()
        self.load.save_prev()
        self.ice.deli_00_prev=self.ice.coeff[0]
        self.oc.save_prev()

    def reset_save_data(self):
        self.SL.delS.saved=np.array([])
        self.SL.sdelS.saved=np.array([])

    def save_data(self):
        self.SL.sdelS.save()
        self.SL.delS.save()

    def calc_delR(self,delL,sdelL):
        if not('beta_R_l' in self.love.__dict__):
            self.love.calc_beta_R(self)
        self.delR_e=np.zeros((self.time_step_number,int((self.maxdeg+1)*(self.maxdeg+2)/2)))+1j
        self.delR_v=np.zeros((self.time_step_number,int((self.maxdeg+1)*(self.maxdeg+2)/2)))+1j
        for t_it in range (1,self.time_step_number):
            print(self.grid.time_step[t_it])
            self.delR_e[t_it,:]=self.love.T.coeff*self.love.h.coeff*delL[t_it,:]
            self.load.calc_viscuous_load(self,t_it,sdelL,self.love.beta_R_l)
            self.delR_v[t_it,:]=self.love.T.coeff*self.load.V_lm.coeff
    

    def save_model(self) :
        
        Data_dict={'g' : self.g,'M_e' : self.M_e,'omega' : self.omega,'k_hydro' : self.k_hydro, 'G' : self.G,'k_f' : self.k_f, 'a' : self.a, 'CminA' : self.CminA, 'C' : self.C, 'maxdeg' : self.maxdeg, 'rotation' : self.rotation, 'k_max' : self.k_max, 'epsilon' : self.epsilon, 'topo_it_max': self.topo_it_max, 'nb_workers' : self.nb_workers, 'time_step' : self.grid.time_step, 'sed_way' : 'Irrawady_sedimentationgrid', 'rho_sed' : self.sed.rho, 'ice_way' : 'ice6g_data', 'ice_corrected' : self.ice.ice_corrected, 'rho_ice' : self.ice.rho, 'delS' : self.SL.delS.saved, 'sdelS' : self.SL.sdelS.saved, 'rho_sea' : self.oc.rho,  'love_way' : 'VM5.l40.um316.lm5.1024/', 'topo_way' : 'topo_SL', 'topo' : self.topo.topo}

        # datetime object containing current date and time
        now = datetime.now()

        # dd/mm/YY H:M:S
        dt_string = now.strftime("%d.%m.%Y_%H.%M.%S")

        with open('result_'+dt_string+'.pkl', 'wb') as f:
            pickle.dump(Data_dict, f)