import numpy as np
import love as Love
from ocean_function import spherical_ocean_function
from grid import GRID
from spharm import legendre, GaussQuad
from plot import world_plot
from multiprocessing import Pool
import par_SL_model as par
import itertools

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
    
    def __init__(self,g=9.80616,maxdeg=64,M_e=5.9742*10**24,omega=7.292*10**(-5),k_hydro=0.934,G=6.67408*10**(-11),k_f=0.942,a = 6371000,C=8.034*10**37,k_max=10,epsilon=10**-4,topo_it_max=10,nb_workers=6):
        """
            Parameters
            ----------
        """
        
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
        self.init_grided_parameters()

    def create_GRID(self):
        self.grid=GRID(self)

    def init_grided_parameters(self): 

        self.P_lm=legendre(self.maxdeg,self.nb_workers) #Calculate the Legendre associated functions for the Gaussian grid
        self.x, self.w= GaussQuad(self.maxdeg) # calculate the Gaussian grid parameters
        self.grid=GRID(self)

        # create a better object for the input and prepare the input for the functions
        self.ice=self.grid.create_ICE() # adapt the load method of the spherical_ice class to the data you are using
        self.ice.quick_load(self.grid,self.nb_workers) # change to quick_load for parallelised version ! ne marche pas. 

        self.sed = self.grid.create_SED()
        self.sed.load(self)

        self.SL=self.grid.create_SL()
        self.love=self.grid.create_LOVE(self)

        self.topo = self.grid.create_TOPO() # adapt the load method of the spherical_topo class to the data you are using
        self.topo.load(self)


        self.oc=spherical_ocean_function().evaluate_ocean(self.topo.topo_pres).grdtocoeff(self)
        self.oc.update_oc_0()# initialize the grid at t=t_it(t_it=0)
        self.oc.area=self.oc.coeff[0]

        

    def run(self):
        for self.topo_it in range(self.topo_it_max): # Il y a un pb sur le chargement de topo à corriger !
            self.topo_loop()
            self.topo.topo_initial=self.ice.topo_pres_ice_corrected - (self.topo.topo[-1,:,:]-self.topo.topo[0,:,:]) 
        return

    def topo_loop(self):
        #regarder si il n'y a pas une boucle while plutot avec un critère de convergence pour stopper la boucle.
        # initialization of the prev variable
        self.SL.reset()

        self.delL_prev=np.zeros(self.love.T.coeff.shape)+1j*0
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
        self.SL.ESL = np.zeros((len(self.grid.time_step),))
        for t_it in range(1,len(self.grid.time_step)):
            print(self.topo_it,t_it)
            self.time_loop(t_it)
        
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
        
        self.SL.delS_prev=self.SL.delS
        self.SL.TO.save_prev()
        self.SL.delL_prev=self.SL.delL
        self.SL.delLa_prev=self.SL.delLa 
        self.ice.deli_00_prev=self.ice.coeff[0]
        self.SL.ESL[t_it]=self.ice.coeff[0]/self.oc.area*self.ice.rho/self.oc.rho
        self.oc.save_prev()

        #SL.ESL[t_it]=ice.deli_00_prev/oc.area*ice.rho/oc.rho
        self.SL.delSL.modify(self.SL.delSLcurl.grd+self.SL.delPhi_g)
        self.topo.topo[t_it,:,:]=-self.SL.delSL.grd+self.topo.topo_0#update cette ligne avec le nouvel objet spherical_top

    def SL_conv_loop(self,t_it,k):
        if k == 0 and self.topo_it==0:
            #print('oc',oc.coeff.shape,'TO',SL.TO.coeff.shape,'TO_prev',SL.TO.prev.shape,'ICE',ice.coeff.shape,'SED',sed.coeff.shape)
            self.SL.sdelS=self.oc.prev/self.oc.prev[0]*(-self.ice.rho/self.oc.rho*self.ice.sdeli_00 + self.SL.TO.coeff[0]-self.SL.TO.prev[0])-self.SL.TO.coeff-self.SL.TO.prev
        self.SL.delS=self.SL.delS_prev + self.SL.sdelS
        self.SL.delL=self.ice.rho*self.ice.coeff + self.oc.rho*self.SL.delS+ self.sed.rho*self.sed.coeff
        self.SL.sdelL[t_it-1]=self.SL.delL-self.SL.delL_prev
        # calculate viscous contribution
        # beta contains the viscous love numbers for time t_it,
        # row index goes over the time increments, column
        # index goes over lm
        #SL.V_lm = np.zeros((int((grid.maxdeg)*(grid.maxdeg+1)/2),))+1j*0
        #SL.V_lm_T = np.zeros((int((grid.maxdeg)*(grid.maxdeg+1)/2),))+1j*0
        if t_it == 1:
            self.SL.V_lm = np.zeros((int((self.grid.maxdeg+1)*(self.grid.maxdeg+2)/2),))+1j*0
        else :
            # for lm_it in range(int((self.grid.maxdeg+1)*(self.grid.maxdeg+2)/2)):
            #     self.SL.V_lm[lm_it] = np.dot(self.love.beta_l[t_it-1,:t_it-1,int(self.love.beta_counter[lm_it]-1)].transpose(), self.SL.sdelL[:t_it-1,lm_it])# regarder comment gérer cet partie !!!            
            results=self.pool.starmap(par.f_V_lm,zip(self.love.beta_l.transpose(),[t_it for i in range(len(self.SL.delL))],self.SL.sdelL.transpose()))
            self.SL.V_lm=np.array(results)
        Love.calc_rot_visc(self,t_it)
        self.SL.sdelLa[t_it-1,:]=self.SL.delLa - self.SL.delLa_prev
        if t_it == 1 :
            self.SL.V_lm_T = np.zeros((int((self.grid.maxdeg+1)*(self.grid.maxdeg+2)/2),))+1j*0
        else :
            # for lm_it  in range(6) : # don't need to loop over all degrees
            #     self.SL.V_lm_T[lm_it] = np.dot(self.love.beta_tide[t_it-1,:t_it-1,int(self.love.beta_counter[lm_it]-1)].transpose(), self.SL.sdelLa[:t_it-1,lm_it]); # pareil réfléchir à comment gérer cette indexation
            results=self.pool.starmap(par.f_V_lm,zip(self.love.beta_tide.transpose(),[t_it for i in range(len(self.SL.delL))],self.SL.sdelLa.transpose()))
            self.SL.V_lm_T=np.array(results)
        # calculate sea level perturbation
        # add ice and sea level and multiply with love numbers
        
        self.SL.delSLcurl_fl=self.love.E.coeff * self.love.T.coeff *self.SL.delL+self.love.T.coeff*self.SL.V_lm+1/self.g*self.love.E_T.coeff*self.SL.delLa+1/self.g*self.SL.V_lm_T
        self.SL.delSLcurl.modify(self.SL.delSLcurl_fl - self.ice.coeff -self.sed.coeff,'coeff').coefftogrd(self)
        self.SL.RO.modify(self.SL.delSLcurl.grd*self.oc.grd).grdtocoeff(self)
        self.SL.delPhi_g = np.real(1/self.oc.coeff[0] * (- self.ice.rho/self.oc.rho*self.ice.coeff[0] - self.SL.RO.coeff[0] + self.SL.TO.coeff[0]))
        self.SL.sdelS_new=self.SL.RO.coeff + self.SL.delPhi_g*self.oc.coeff -  self.SL.TO.coeff - self.SL.delS_prev
        chi = np.abs((np.sum(np.abs(self.SL.sdelS_new)) - np.sum(np.abs(self.SL.delS))) / np.sum(np.abs(self.SL.delS)))
        self.SL.sdelS=self.SL.sdelS_new.copy()
        return chi

    def calc_RSL(self):
        self.SL.RSL=np.zeros(self.topo.topo.shape)
        for t_it in range(len(self.grid.time_step)):
            self.SL.RSL[t_it,:,:]=self.topo.topo[-1,:,:]-self.ice.ice_corrected[-1,:,:]-self.sed.sed[-1,:,:]-(self.topo.topo[t_it,:,:]-self.ice.ice_corrected[t_it,:,:]-self.sed.sed[t_it,:,:])

