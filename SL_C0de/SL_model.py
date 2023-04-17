import numpy as np
from grid import GRID
from spharm import legendre, GaussQuad
from plot import world_plot
from multiprocessing import Pool
from spharm import sphericalobject
from Load import LOAD
import pickle
from datetime import datetime
import tomli
import warnings
import joblib

class SL_MODEL(world_plot):
    """
    A class used to manage all parameters of the Sea level equation

    Attributes
    ----------

    Parameters
    ----------
        parameter_file : str
            The way to the parameters of the model. To have details about these parameters you can refer to the example toml file.

        empty : bool 
            This parameter give you the possibility of not loading data from a TOML file. The model will be empty but you can use the set_from_file() method to fill it based on a .pkl file. 
        
    Methods
    -------
        create_GRID() 
    """
    def __init__(self,parameter_file='Constant.toml',empty=False):
        if empty :
            warnings.warn('You have created an empty file. To fill it, please use the methods set_from_file.')
        else :
            self.set_from_parameter_file(parameter_file)
            self.P_lm=legendre(self.maxdeg,self.pool) #Calculate the Legendre associated functions for the Gaussian grid
            self.x, self.w= GaussQuad(self.maxdeg) # calculate the Gaussian grid parameters
            self.grid=GRID(self)
            self.create_parameters()

    def set_from_parameter_file(self,parameter_file='Constant.toml'):
        '''
        This method is used to load parameters from a parameters_file. This file is based on a template presented in the example_constant.toml.

        Attributes
        ----------

        Parameters
        ----------
            parameter_file : str
                The way to the parameters of the model. To have details about these parameters you can refer to the example toml file.
        
        '''
        with open(parameter_file, mode="rb") as fp:
            config = tomli.load(fp)
        self.rotation=True
        self.g = config['Earth_param']['g'] # gravitational constant of earth m**2/s
        self.maxdeg=config['Model_param']['maxdeg'] # maximum degree of the spherical harmonic used
        self.M_e = config['Earth_param']['M_e'] # Earth mass [kg]
        self.omega = config['Earth_param']['omega'] # Earth's angular velocity rad/s
        self.k_hydro = config['Earth_param']['k_hydro'] # degree-2 hydrostatic Love number
        self.G=config['Earth_param']['G'] # Newton gravitational constant
        self.k_f=config['Earth_param']['k_f'] # !!!! A trouver !!!!
        self.a = config['Earth_param']['a'] # earth radius (m)
        self.CminA = (self.k_f*self.a**5*self.omega**2)/(3*self.G) # a, c = smallest and largest principal moments of the inertia tensor
        self.C=config['Earth_param']['C']
        self.k_max = config['Model_param']['k_max'] # maximum number of convergence iteration for the resolution of the sea level equation at 1 time
        self.epsilon = config['Model_param']['epsilon'] # convergence criterion for the convergence iteration of the sea level equation at 1 time
        self.topo_it_max = config['Model_param']['topo_it_max'] # maximum numbre of iteration for the convergence of the initial topography
        self.del_topo=config['Model_param']['del_topo']
        data_way=config['Input_data']['data_way']
        self.ice_way=data_way+config['Input_data']['ice_way']
        self.sed_way=data_way+config['Input_data']['sed_way']
        self.love_way=data_way+config['Input_data']['love_way']
        self.topo_way=data_way+config['Input_data']['topo_way']
        self.time_step=config['Time']['time_step']
        self.time_step_number=len(self.time_step)-1

    def set_from_output_file(self,output_file,way_to_way='',large=False):
        """
        This function is used to load previous computation of the model. This can be used mainly to calculate other parameters of the model for the post process part. This method use the pickle module to load the pkl file if it's a small file of joblib if it's a big one. 

        Attributes
        ----------
            output_file : str
                The way to the .pkl file where the output of a previous computation of the model are stored.
            way_to_way : str
                You can use this parameter to change the way to the stored data (sedimentation grid, ice, topography). It'll be a  $toml file like the one for the parametrisation but only with the Input_data part (example): 

            [Input_data]
            ice_way='your_ice_data'
            topo_way='your_topographic_data'
            sed_way='your_sedimentation_grid'
            love_way=your_love_numbers'

            large : bool
                This parameter let you define if it's a large or a small .pkl file. This is nessessary because of the different modules used to save these files. 
        """
        self.load_output_file(output_file,large,'model')

        with open(way_to_way, mode="rb") as fp:
            config = tomli.load(fp)
        data_way=config['Input_data']['data_way']
        self.ice_way=data_way+config['Input_data']['ice_way']
        self.sed_way=data_way+config['Input_data']['sed_way']
        self.love_way=data_way+config['Input_data']['love_way']
        self.topo_way=data_way+config['Input_data']['topo_way']

        self.P_lm=legendre(self.maxdeg,self.pool) #Calculate the Legendre associated functions for the Gaussian grid
        self.x, self.w= GaussQuad(self.maxdeg) # calculate the Gaussian grid parameters

        self.create_GRID() #create the grid framework

        self.create_parameters() #Create all grid parameters


        self.ice.quick_load(self.grid,self.pool,self.ice_way) #set the ice self
        self.sed.load(self,self.sed_way) # set the sedimentation self
        self.topo.load(self,self.topo_way) # set the topography self, this function need the sediment and ice self ! 
        #Disk_load=self.grid.disk(self,0,180,1,100)
        # we need to calculate the ocean function after setting the topography due to it's dependencies
        self.create_ocean(self.topo.topo_pres) # create the ocean function (see theoria (TABOO), spada et al., 2003) from the actual topography.
        self.load_output_file(output_file,large,'module')

    def load_output_file(self,output_file,large=False,type='model'):
        """
        This method is the base method to load parameters from a toml file. It's called by set_from_file and the init function. 

        Attributes
        ----------
            output_file : str
                The way to the .pkl file where the output of a previous computation of the model are stored.

            large : bool
                This parameter let you define if it's a large or a small .pkl file. This is nessessary because of the different modules used to save these files. 

            type : (str)
                Can be 'model' or 'module'. Define the type of data you want to load, there is the model parameters where we have the basic toml file. And then we have the module_parameters with the density of the ice, sediment etc. 
        """
        if type=='model' :
            if large :
                with open(output_file+'.pkl', 'rb') as handle:
                    data = joblib.load(handle)
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
                self.love_way=data['love_way'][:-1]
                self.topo_way=data['topo_way']
            else :
                with open(output_file+'.pkl', 'rb') as handle:
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
                self.love_way=data['love_way'][:-1]
                self.topo_way=data['topo_way']
        elif type=='module': 
            if large :
                with open(output_file+'.pkl', 'rb') as handle:
                    data = joblib.load(handle)
                self.sed.rho=data['rho_sed']
                self.oc.rho=data['rho_sea']
                self.ice.rho=data['rho_ice']
                self.grid.time_step=self.time_step
                self.ice.ice_corrected=data['ice_corrected']
                self.SL.delS.saved=data['delS']
                self.SL.sdelS.saved=data['sdelS']
                self.topo.topo=data['topo']
            else :
                with open(output_file+'.pkl', 'rb') as handle:
                    data = pickle.load(handle)
                self.sed.rho=data['rho_sed']
                self.oc.rho=data['rho_sea']
                self.ice.rho=data['rho_ice']
                self.grid.time_step=self.time_step
                self.ice.ice_corrected=data['ice_corrected']
                self.SL.delS.saved=data['delS']
                self.SL.sdelS.saved=data['sdelS']
                self.topo.topo=data['topo']
        else :
            warnings.warn('Wrong argument for the type parameter : can be either "model" or "module"')
                


    def create_GRID(self):
        """
        This function is used to create the grids parameters. refer to the grid module to have more informations. 
        """
        self.grid=GRID(self)
    
    def create_ocean(self,topo):
        """
        This function is used to create the ocean object wich contains all the paramters of the ocean plus method to calculate. This module more than creating the ocean, it preset it with the actual state of the ocean based on the actual topography. It calcultae the actual ocean area wich is used for the calculation of the ESL curve. For more information refer to the module ocean_function.
        """
        self.oc=spherical_ocean_function().evaluate_ocean(topo).grdtocoeff(self)
        self.oc.update_oc_0()# initialize the grid at t=t_it(t_it=0)
        self.oc.area=self.oc.coeff[0]
    
    def create_LOAD(self):
        """
        This function is used to create the loads parameters. For more information refer to the Load module.
        """
        return LOAD(self)

    def create_parameters(self):
        """
        This function create all the model parameters based on the grids parameters. It interpolate the ice topographic and sediment grid on the model grid. For more information refer to the init method of the class of the following modules : ice, SeaLevel, sediment, love, topography and load. 
        """
        self.ice=self.grid.create_ICE() # adapt the load method of the spherical_ice class to the data you are using
        self.sed = self.grid.create_SED()
        self.SL=self.grid.create_SL()
        self.grid.time_step=self.time_step
        self.love=self.grid.create_LOVE(self,self.love_way)
        self.topo = self.grid.create_TOPO() # adapt the load method of the spherical_topo class to the data you are using
        self.load=self.create_LOAD()
    
    def reset(self):
        """
        This method can be use to reset to zero some parameters used in the code. This function is called each time we restart the computation after topographic correction of the previous itteration. For more information refere to the model workflow !!!A faire !!! .
        """
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

    def save_prev(self):
        """
        This method is used to save a parameter of the previous interation at each time step. We use this function to simplify the loop when we are resolving the SL equation.
        """
        self.SL.delS.save_prev()
        self.SL.TO.save_prev()
        self.load.save_prev()
        self.ice.deli_00_prev=self.ice.coeff[0]
        self.oc.save_prev()



    def save_model(self,output_way=r'C:/Users/ahenry01/Desktop/Python_code/SL_C0de_output/',subname='') :

        '''
        This method is used to save the model output and parameters. Among these outputs you have the Sea level variation, the corrected ice and the dynamic topography.
        
        Attributes
        ----------
            output_way : str
                The way to the output file.
            subname : str
                A suplementary string to precise the output.
        '''
        
        Data_dict={'g' : self.g,'M_e' : self.M_e,'omega' : self.omega,'k_hydro' : self.k_hydro, 'G' : self.G,'k_f' : self.k_f, 'a' : self.a, 'CminA' : self.CminA, 'C' : self.C, 'maxdeg' : self.maxdeg, 'rotation' : self.rotation, 'k_max' : self.k_max, 'epsilon' : self.epsilon, 'topo_it_max': self.topo_it_max, 'nb_workers' : self.nb_workers, 'time_step' : self.grid.time_step, 'sed_way' : self.sed_way, 'rho_sed' : self.sed.rho, 'ice_way' : self.ice_way, 'ice_corrected' : self.ice.ice_corrected, 'rho_ice' : self.ice.rho, 'delS' : self.SL.delS.saved, 'sdelS' : self.SL.sdelS.saved, 'rho_sea' : self.oc.rho,  'love_way' : self.love_way, 'topo_way' : self.topo_way, 'topo' : self.topo.topo}

        # datetime object containing current date and time
        now = datetime.now()

        # dd/mm/YY H:M:S
        dt_string = now.strftime("%d.%m.%Y_%H.%M.%S")

        with open(output_way+subname+'result_'+dt_string+'.pkl', 'wb') as f:
            pickle.dump(Data_dict, f)