from .grid import ICE_TIME_GRID
from .grid import SEDIMENT_TIME_GRID
from .grid import TOPOGRAPHIC_TIME_GRID
from .grid import OCEAN_TIME_GRID
from .grid import LOAD_TIME_GRID
from .spharm import sphericalobject
from .love import LOVE
import numpy as np
import copy
import os
import sys

def Precomputation(ice_grid,sed_grid,topo_grid,Output_way,stop=26,step=0.5,maxdeg=512,irregular_time_step=None,backend=False,plot=False):
    '''
    The _`Precomputation` method prepare the different data to match temporal and spatial resolution of the modelisation. spatial resolution (m,n) and temporal resolution (step). The input data format is described in :ref:`Mass grid format<grid_format>`.

    Attribute :
    ----------- 
        ice_grid : dict(name = str, grid = np.array((n_i,m_i,t_i)), time = np.array((t_i,)), lon = np.array((n_i)), lat = np.array((m_i)))
            ice thickness grid
        sed_grid : dict(name = str, grid = np.array((n_j,m_j,t_j)), time = np.array((t_j,)), lon = np.array((n_j)), lat = np.array((m_j)),frame = (lon1,lon2,lat1,lat2))
            sediment thickness
        topo_grid : dict(name = str, grid = np.array((n_k,m_k)), lon = np.array((n_j)), lat = np.array((m_j)))
            present topography 
        Output_way : str
            The filepath toward wich the different grid will be saved
        stop : float 
            age at wich the computation stop
        step : float
            time step used for the temporal resolution of the model
        maxdeg : int 
            The maximum degree of spherical harmonics, that define the spatial resolution of the model. (m = 2*maxdeg, n = maxdeg)
        irregular_time_step : np.array((time_step_number,))
            An irregular time array to compute the model over non regular time step
        backend : bool
            If required you can ask a backend during the run of each computation. True, will generate a backend, False will not. Defaul is False
        Plot : bool
            If needed you can ask plot of each output compared to the initial plot. You will need the cartopy doc.
    Returns : 
    --------- 
        None 

    The data output of the function are saved to the output way.

    '''

    time_step=np.arange(start=stop,stop=-step,step=-step)
    if time_step[-1]<0:
        time_step[-1]=0
    if not(irregular_time_step is None):
        time_step=irregular_time_step

    ice_time_grid=ICE_TIME_GRID(time_step,maxdeg,grid_name=ice_grid['name'])
    ice_time_grid.interp_on_time_and_space(ice_grid['grid'],ice_grid['time'],ice_grid['lon'],ice_grid['lat'],grid_type='global',backend=backend)
    ice_time_grid.save(save_way=Output_way)

    sed_time_grid=SEDIMENT_TIME_GRID(time_step,maxdeg,grid_name=sed_grid['name'])
    sed_time_grid.interp_on_time_and_space(sed_grid['grid'],sed_grid['time'],sed_grid['lon'],sed_grid['lat'],backend=backend,grid_type='local')
    sed_time_grid.height_time_grid=sed_time_grid.height_time_grid=sed_time_grid.height_time_grid[:,::-1,:]# latitude are inverted ! 
    sed_time_grid.save(save_way=Output_way)


    topo_time_grid=TOPOGRAPHIC_TIME_GRID(time_step,maxdeg,grid_name=topo_grid['name'])
    topo_time_grid.topo_pres=topo_time_grid.interp_on(topo_grid['grid'],topo_grid['lon'],topo_grid['lat'],grid_type='global')+ice_time_grid.height_time_grid.sum(0)

    # there is no ice in topo_sl

    for i in range(topo_time_grid.time_step_number-1):
        if i==0:
            topo_time_grid.height_time_grid[i,:,:]=topo_time_grid.topo_pres-sed_time_grid.height_time_grid.sum(0)-ice_time_grid.height_time_grid[1:,:,:].sum(0)
        else :
            topo_time_grid.height_time_grid[i,:,:]=topo_time_grid.topo_pres-sed_time_grid.height_time_grid[i:,:,:].sum(0)-ice_time_grid.height_time_grid[i:,:,:].sum(0)
    
    topo_time_grid.save(save_way=Output_way)

    if plot :
        plot_precomputation(ice_grid['time'],ice_grid['grid'],ice_time_grid,save_way=Output_way)
        plot_precomputation(sed_grid['time'],sed_grid['grid'],sed_time_grid,area=sed_grid['frame'],lon_init=sed_grid['lon'],lat_init=sed_grid['lat'],save_way=Output_way)
        colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 256))
        colors_land = plt.cm.terrain(np.linspace(0.25, 1, 256))
        all_colors = np.vstack((colors_undersea, colors_land))
        terrain_map = colors.LinearSegmentedColormap.from_list('terrain_map',all_colors)
        plot_precomputation(np.array([1]),topo_grid['grid'],topo_time_grid,topo=True,save_way=Output_way,vcenter=0,cmap=terrain_map)

def SLE_forward_modeling(Input_way,ice_name,sed_name,topo_name,ocean_name,love_way,love_file,conv_lim,Output_way):
    '''
    The _`SLE_forward_modeling` method solve the SLE with no constrain on final topography. This result in a forward modeling of the SL. It can be used for test and exploration of models. For informations on iterations see :ref:`description of iteration <iteration_desc>`.

    Attribute : 
    ------------
        Input_way : str
            File location of the input data (ice, sediment and topography)
        ice_name : str
            The name of the ice data file
        sed_name : str
            The name of the sediment data file
        topo_name : str
            The name of the topographic data file
        ocean_name : str
            The name of the output file containig the data on the ocean
        love_way : str 
            way of the love numbers file used to compute earth visco elastic deformation
        conv_lim : float
            Limit of precision required for the convergence of the SLE resolution (10^-3)
        output_way : str 
            filepath for saving results of the modelisation

    Returns :
    --------- 
        None 

    The resulting ocean class object (:ref:`OCEAN_TIME_GRID <OCEAN_TIME_GRID>`) containing the ocean thickness is saved in the outputway under the name ocean. 
    '''
    ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+'/'+ice_name))
    sed_time_grid=SEDIMENT_TIME_GRID(from_file=(True,Input_way+'/'+sed_name))
    topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/'+topo_name))

    ocean_time_grid=OCEAN_TIME_GRID(ice_time_grid.time_step,ice_time_grid.maxdeg,grid_name=ocean_name)
    ocean_time_grid.time_step_number=ocean_time_grid.time_step_number

    love_number=LOVE(ice_time_grid.maxdeg,love_way+'/'+love_file,ice_time_grid.time_step,6371000,5.9742e24,type='time')
    love_number.dev_beta()

    TO=sphericalobject(coeff=np.zeros(ice_time_grid.height_time_coeff[0,:].shape))
    TO.prev=np.zeros(ice_time_grid.height_time_coeff[0,:].shape)

    

    ice_time_grid.ice_correction(topo_time_grid,ocean_time_grid)
    ice_time_grid.timegrdtotimecoeff()

    load=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff[1:t_it,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[1:t_it,:]*sed_time_grid.rho+ocean_time_grid.height_time_coeff[1:t_it,:]*ocean_time_grid.rho,beta_l=love_number.beta_l,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg)

    topo_time_grid.update_0()
    ocean_time_grid.evaluate_ocean(topo_time_grid.grd_0).grdtocoeff()
    ocean_time_grid.update_0()
    ocean_time_grid.save_prev()
    topo_time_grid.grid_from_step(0)
    # grd correspond donc au topo_j défini dans le code de kendal et al.
    ocean_time_grid.evaluate_ocean(topo_time_grid.grd).grdtocoeff()
    TO.grd=topo_time_grid.grd_0*(ocean_time_grid.grd-ocean_time_grid.grd_0)
    TO.grdtocoeff()

    track_conv=np.array([])

    for t_it in range (1,ice_time_grid.time_step_number-1):
        topo_time_grid.grid_from_step(t_it)
        # grd correspond donc au topo_j défini dans le code de kendal et al.
        ocean_time_grid.evaluate_ocean(topo_time_grid.grd).grdtocoeff()
        TO.grd=topo_time_grid.grd_0*(ocean_time_grid.grd-ocean_time_grid.grd_0)
        TO.grdtocoeff()
        sed_time_grid.coeff_from_step(t_it)
        ice_time_grid.coeff_from_step(t_it)
        
        conv_it=0
        conv_it=ocean_time_grid.sea_level_solver(load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,conv_lim)

        track_conv=np.append(track_conv,conv_it)     
        
        TO.prev=TO.coeff.copy()

        ocean_time_grid.save_prev()
        
        topo_time_grid.height_time_grid[t_it,:,:]=topo_time_grid.grd_0-(ocean_time_grid.delSLcurl.grd+ocean_time_grid.delPhi_g)
    
    if not(os.path.exists(Output_way+'/model_output')):
            os.mkdir(Output_way+'/model_output')
    Output_way=Output_way+'/model_output'

    if not(os.path.exists(Output_way+'/'+love_file)):
        os.mkdir(Output_way+'/'+love_file)
    ocean_time_grid.save(Output_way+'/'+love_file)
    ice_time_grid.save(Output_way+'/'+love_file)
    topo_time_grid.save(Output_way+'/'+love_file)
    print(sed_time_grid.time_grid_name)
    sed_time_grid.save(Output_way+'/'+love_file)

def SLE_solver(Input_way,ice_name,sed_name,topo_name,ocean_name,love_way,love_file,topo_lim,conv_lim,Output_way):
    '''
    The _`SLE_solver` solve the SLE and converge toward the actual topography. The computation of the SLE resolution is describe in :ref:`Sea level equation resolution <SLE_res>`. For informations on iterations see :ref:`description of iteration <iteration_desc>`.

    Attribute :
    ----------- 
        Input_way : str
            File location of the input data (ice, sediment and topography)
        ice_name : str
            The name of the ice data file
        sed_name : str
            The name of the sediment data file
        topo_name : str
            The name of the topographic data file
        ocean_name : str
            The name of the output file containig the data on the ocean
        love_way : str 
            way of the love numbers file used to compute earth visco elastic deformation
        conv_lim : float
            Limit of precision required for the convergence of the SLE resolution (10^-3)
        topo_lim : float
            Topography convergence limit toward wich the iteration will converge (1 m)
        output_way : str 
            filepath for saving results of the modelisation

    Returns : 
    --------- 
        None 

    The resulting ocean class object (:ref:`OCEAN_TIME_GRID <OCEAN_TIME_GRID>`), updated ice thickness (:ref:`ICE_TIME_GRID <ICE_TIME_GRID>`) and updated  topography (:ref:`TOPOGRAPHIC_TIME_GRID <TOPOGRAPHIC_TIME_GRID>`) are saved to output_way. 
    '''
    
    ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+'/'+ice_name))
    sed_time_grid=SEDIMENT_TIME_GRID(from_file=(True,Input_way+'/'+sed_name))
    topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/'+topo_name))

    maxdeg=ice_time_grid.maxdeg

    time_step=ice_time_grid.time_step.copy()

    # Initiate the base elements
    from SL_C0de.spharm import sphericalobject
    from SL_C0de.Load import LOAD

    ocean_time_grid=OCEAN_TIME_GRID(time_step,maxdeg,grid_name=ocean_name)
    ocean_time_grid.time_step_number=ocean_time_grid.time_step_number

    ice_time_grid.timegrdtotimecoeff()
    sed_time_grid.timegrdtotimecoeff()

    love_number=LOVE(maxdeg,love_way+'/'+love_file,time_step,6371000,5.9742e24)
    love_number.dev_beta()
    love_number.dev_beta_tide()
    TO=sphericalobject(coeff=np.zeros(ice_time_grid.height_time_coeff[0,:].shape))


    topo_diff_median=np.inf
    #sdel_topo_diff=np.inf
    topo_it=0
    while topo_diff_median>topo_lim : #and sdel_topo_diff>10**(-1):
        delPhi_g_time=np.array([])
        TO.prev=np.zeros(ice_time_grid.height_time_coeff[0,:].shape)

        if topo_diff_median != np.inf :
            topo_time_grid.height_time_grid[0,:,:]=topo_initial.copy()
        # topo_time_grid.height_time_grid = topo_time_grid.height_time_grid - ice_time_grid.height_time_grid + ice_time_grid.ice # for resetting the corrected ice.
        ice_time_grid.ice_correction(topo_time_grid,ocean_time_grid)
        ice_time_grid.timegrdtotimecoeff()

        t_it=-1

        load=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff[1:t_it,:]*ice_time_grid.rho+sed_time_grid.height_time_coeff[1:t_it,:]*sed_time_grid.rho+ocean_time_grid.height_time_coeff[1:t_it,:]*ocean_time_grid.rho,beta_l=love_number.beta_l,E=love_number.E,E_T=love_number.E_T,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg)
        
        topo_time_grid.update_0()
        ocean_time_grid.evaluate_ocean(topo_time_grid.grd_0).grdtocoeff()
        ocean_time_grid.update_0()
        ocean_time_grid.save_prev()
        topo_time_grid.grid_from_step(0)
        # grd correspond donc au topo_j défini dans le code de kendal et al.
        ocean_time_grid.evaluate_ocean(topo_time_grid.grd).grdtocoeff()
        TO.grd=topo_time_grid.grd_0*(ocean_time_grid.grd-ocean_time_grid.grd_0)
        TO.grdtocoeff()

        track_conv=np.array([])

        for t_it in range (1,ice_time_grid.time_step_number-1):
            topo_time_grid.grid_from_step(t_it)
            # grd correspond donc au topo_j défini dans le code de kendal et al.
            ocean_time_grid.evaluate_ocean(topo_time_grid.grd).grdtocoeff()
            TO.grd=topo_time_grid.grd_0*(ocean_time_grid.grd-ocean_time_grid.grd_0)
            TO.grdtocoeff()
            sed_time_grid.coeff_from_step(t_it)
            ice_time_grid.coeff_from_step(t_it)
            


            print('time_iteration : ',ice_time_grid.time_step[t_it])
            if topo_it==0 : 
                conv_it=0
            else :
                conv_it=1
            conv_it=ocean_time_grid.sea_level_solver(load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,conv_lim)

            print(conv_it)

            track_conv=np.append(track_conv,conv_it)     
            
            TO.prev=TO.coeff.copy()

            ocean_time_grid.save_prev()
            
            topo_time_grid.height_time_grid[t_it,:,:]=topo_time_grid.grd_0-(ocean_time_grid.delSLcurl.grd+ocean_time_grid.delPhi_g)
            

            #delPhi_g_time=np.append(delPhi_g_time,ocean_time_grid.delPhi_g)

        topo_it+=1
        topo_pres_ice_corrected=topo_time_grid.topo_pres-ice_time_grid.ice.sum(0)+ice_time_grid.height_time_grid.sum(0)
        topo_diff=np.abs(topo_time_grid.height_time_grid[-1,:,:]-topo_pres_ice_corrected).max().max()
        sdel_topo_diff=np.abs(topo_diff-np.abs(topo_time_grid.height_time_grid[-1,:,:]-topo_pres_ice_corrected).max().max())
        topo_diff_mean=np.abs(topo_time_grid.height_time_grid[-1,:,:]-topo_pres_ice_corrected).mean().mean()
        topo_diff_median=np.median(np.median(np.abs(topo_time_grid.height_time_grid[-1,:,:]-topo_pres_ice_corrected)))
        print(topo_it,' : ',topo_diff, topo_diff_mean, topo_diff_median, sdel_topo_diff,track_conv-(topo_it>0))
        topo_initial=topo_pres_ice_corrected - (topo_time_grid.height_time_grid[-1,:,:]-topo_time_grid.height_time_grid[0,:,:])

    if not(os.path.exists(Output_way+'/model_output')):
        os.mkdir(Output_way+'/model_output')
    Output_way=Output_way+'/model_output'

    if not(os.path.exists(Output_way+'/'+love_file)):
        os.mkdir(Output_way+'/'+love_file)
    ocean_time_grid.save(Output_way+'/'+love_file)
    ice_time_grid.save(Output_way+'/'+love_file)
    topo_time_grid.save(Output_way+'/'+love_file)
    sed_time_grid.save(Output_way+'/'+love_file)
    return load

def find_files(filename, search_path):
    result = []

    # Wlaking top-down from the root
    for root, dir, files in os.walk(search_path):
        if filename in dir:
            result.append(os.path.join(root, filename))
    return result

def calculate_deformation(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,Output_way,backend=False,plot=dict(plot=False)) :

    '''
    The _`calculate_deformation` method calculate the earth viscoelastic response resulting from the different masses caculated with the :ref:`SLE_solver <SLE_solver>`. This computation follows the decomposition of the SLE described in :ref:`Computation of ground and geoid subsidence from different load source <G_R_comp>`

    Attribute :
    ----------- 
        love_number : :ref:`LOVE <LOVE>` 
            The loaded love number in the form of a LOVE class object
        ice_time_grid : :ref:`grid.ICE_TIME_GRID <ICE_TIME_GRID>`
            ice thickness in the form of ICE_TIME_GRID class object
        sed_time_grid : :ref:`grid.SEDIMENT_TIME_GRID <SEDIMENT_TIME_GRID>` 
            sediment thickness in the form of SEDIMENT_TIME_GRID class object
        ocean_time_grid : :ref:`grid.OCEAN_TIME_GRID <OCEAN_TIME_GRID>` 
            ocean thickness calculated using the SLE solver in the form of OCEAN_TIME_GRID class object
        a : float
            radius of earth in meter
        Me : float
            mass of earth (kg)
        Output_way : str
            way of the output where the load is calculated
        backend : bool
            set if the function is writing its state of computation

    Returns :  
    ---------
        None 

    The resulting LOAD and GEOID are stored in a LOAD file in the same file then the input files. The structure is based on the :ref:`LOAD_TIME_GRID <LOAD_TIME_GRID>`  class object.
    '''
    beta_l=love_number.beta_R_l

    ice_load_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='ICE_LOAD')
    ice_load_time_grid.calc_elastic_time()
    ice_load_time_grid.calc_viscuous_time(backend=backend)
    ice_load_time_grid.save(save_way=Output_way)    
    ice_load_time_grid=0

    sediment_load_time_grid=LOAD_TIME_GRID(sdelL=sed_time_grid.height_time_coeff*sed_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='SEDIMENT_LOAD')
    sediment_load_time_grid.calc_elastic_time()
    sediment_load_time_grid.calc_viscuous_time(backend=backend)
    sediment_load_time_grid.save(save_way=Output_way)
    sediment_load_time_grid=0

    ocean_load_time_grid=LOAD_TIME_GRID(sdelL=ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEAN_LOAD')
    ocean_load_time_grid.calc_elastic_time()
    ocean_load_time_grid.calc_viscuous_time(backend=backend)
    ocean_load_time_grid.save(save_way=Output_way)
    ocean_load_time_grid=0

    total_load_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho+sed_time_grid.height_time_coeff*sed_time_grid.rho+ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='TOTAL_LOAD')
    total_load_time_grid.calc_elastic_time()
    total_load_time_grid.calc_viscuous_time(backend=backend)
    total_load_time_grid.save(save_way=Output_way)
    total_load_time_grid=0

    #Calculating the geoïd deformation

    beta_l=love_number.beta_G_l

    ice_geoid_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho,beta_l=beta_l,E=1+love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='ICE_GEOID')
    ice_geoid_time_grid.calc_elastic_time()
    ice_geoid_time_grid.calc_viscuous_time(backend=backend)
    ice_geoid_time_grid.save(save_way=Output_way)
    ice_geoid_time_grid=0

    sediment_geoid_time_grid=LOAD_TIME_GRID(sdelL=sed_time_grid.height_time_coeff*sed_time_grid.rho,beta_l=beta_l,E=1+love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='SEDIMENT_GEOID')
    sediment_geoid_time_grid.calc_elastic_time()
    sediment_geoid_time_grid.calc_viscuous_time(backend=backend)
    sediment_geoid_time_grid.save(save_way=Output_way)
    sediment_geoid_time_grid=0

    ocean_geoid_time_grid=LOAD_TIME_GRID(sdelL=ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEAN_GEOID')
    ocean_geoid_time_grid.calc_elastic_time()
    ocean_geoid_time_grid.calc_viscuous_time(backend=backend)
    ocean_geoid_time_grid.save(save_way=Output_way)
    ocean_geoid_time_grid=0

    total_geoid_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho+sed_time_grid.height_time_coeff*sed_time_grid.rho+ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='TOTAL_GEOID')
    total_geoid_time_grid.calc_elastic_time()
    total_geoid_time_grid.calc_viscuous_time(backend=backend)
    total_geoid_time_grid.save(save_way=Output_way)
    total_geoid_time_grid=0



def calculate_sediment_ocean_interaction(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,topo_time_grid,Output_way,backend=False,plot=dict(plot=False)) :

    '''
    The _`calculate_sediment_ocean_interaction` method calculate the earth viscoelastic response resulting from sediment under sea surface. This is used to retrieve the effect of ocean replacement by sediment on the sediment load. The resulting subsidence is the true subsidence induced by sediment. The problems and resolution in equation is described in :ref:`True sediment subsidence <sed_subs>`.

    Attribute :
    ----------- 
        love_number : :ref:`LOVE <LOVE>` 
            The loaded love number in the form of a LOVE class object
        ice_time_grid : :ref:`grid.ICE_TIME_GRID <ICE_TIME_GRID>`
            ice thickness in the form of ICE_TIME_GRID class object
        sed_time_grid : :ref:`grid.SEDIMENT_TIME_GRID <SEDIMENT_TIME_GRID>` 
            sediment thickness in the form of SED_TIME_GRID class object
        ocean_time_grid : :ref:`grid.OCEAN_TIME_GRID <OCEAN_TIME_GRID>` 
            ocean thickness calculated using the SLE solver in the form of OCEAN_TIME_GRID class object
        a : float
            radius of earth in meter
        Me : float
            mass of earth (kg)
        Output_way : str
            way of the output where the load is calculated
        backend : bool
            set if the function is writing its state of computation

    Returns : 
    --------- 
        None 

    The resulting LOAD and GEOID are stored in a LOAD file in the same file then the input files. The structure is based on the :ref:`LOAD_TIME_GRID <LOAD_TIME_GRID>` class object.
    '''

    beta_l=love_number.beta_R_l # Load the earth love numbers
    #Preparing a new grid that compute the load of the substracted sediment volume to the ocean.
    #We create a new object to compute the oceanic sediment
    oceanic_sediment_time_grid=copy.copy(sed_time_grid)
    oceanic_sediment_time_grid.rho=ocean_time_grid.rho
    

    for t_it in range (oceanic_sediment_time_grid.time_step_number-1): # At each time step we apply the ocean function to the sediment height grid
        oceanic_sediment_time_grid.height_time_grid[t_it,:,:]=oceanic_sediment_time_grid.height_time_grid[t_it,:,:]*ocean_time_grid.evaluate_ocean(topo_time_grid.height_time_grid[t_it,:,:]).grd
    oceanic_sediment_time_grid.timegrdtotimecoeff()


    oceanic_sediment_load_time_grid=LOAD_TIME_GRID(sdelL=oceanic_sediment_time_grid.height_time_coeff*oceanic_sediment_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=oceanic_sediment_time_grid.time_step,maxdeg=oceanic_sediment_time_grid.maxdeg,grid_name='OCEANIC_SEDIMENT_LOAD')    
    #computing the earth deformation
    oceanic_sediment_load_time_grid.calc_elastic_time()
    oceanic_sediment_load_time_grid.calc_viscuous_time(backend=backend)
    oceanic_sediment_load_time_grid.save(save_way=Output_way)
    oceanic_sediment_load_time_grid.clean_memory()

    beta_l=love_number.beta_G_l # Load the ocean love numbers

    #computing the geoid deformation
    oceanic_sediment_geoid_time_grid=LOAD_TIME_GRID(sdelL=oceanic_sediment_time_grid.height_time_coeff*oceanic_sediment_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=oceanic_sediment_time_grid.time_step,maxdeg=oceanic_sediment_time_grid.maxdeg,grid_name='OCEANIC_SEDIMENT_GEOID')
    oceanic_sediment_geoid_time_grid.calc_elastic_time()
    oceanic_sediment_geoid_time_grid.calc_viscuous_time(backend=backend)
    oceanic_sediment_geoid_time_grid.save(save_way=Output_way)
    oceanic_sediment_geoid_time_grid.clean_memory()



def Post_process(Input_way,sed_name,ice_name,ocean_name,topo_name,love_way):
    '''
    The _`Post_process` calculate the earth viscoelastic response resulting from the different masses caculated with the SLE_solver. This computation is made for all models available in the files. 

    Attribute :
    ----------- 
        input_way_sed : str
            way where the sediment precomputed are stored
        input_way_model_output : str
            filepath to the different masses calculated by the SLE solver
        love_way : str
            way to the love number output

    Returns : 
    --------- 
        None 

    The resulting LOAD and GEOID are stored in a LOAD file in the same file then the input files. The structure is based on the :ref:`LOAD_TIME_GRID <LOAD_TIME_GRID>` class object.
    '''
    


    earth_model_name_list=os.listdir(Input_way)

    for earth_model_name in earth_model_name_list :

        print('calculation for : ' + earth_model_name)

        Output_way=Input_way+'/'+earth_model_name+'/LOAD/'

        if not(os.path.exists(Output_way)):
            os.mkdir(Output_way)

        sed_time_grid=SEDIMENT_TIME_GRID(from_file=(True,Input_way+'/'+earth_model_name+'/'+sed_name))
        sed_time_grid.timegrdtotimecoeff()

        ocean_time_grid=OCEAN_TIME_GRID(from_file=(True,Input_way+'/'+earth_model_name+'/'+ocean_name))
        #print(ocean_time_grid.height_time_coeff.shape)
        ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+'/'+earth_model_name+'/'+ice_name))
        #print(ice_time_grid.time_step_number)

        np.set_printoptions(threshold=sys.maxsize)
        a=6371000
        Me=5.9742e24
        love_way_found=find_files(earth_model_name,love_way)[0]
        love_number=LOVE(ice_time_grid.maxdeg,love_way_found,ice_time_grid.time_step,a,Me)

        ice_time_grid.timegrdtotimecoeff()

        calculate_deformation(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,Output_way)

        #Loading the topography to compute the ocean function at each time step
        topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/'+earth_model_name+'/'+topo_name))

        print('Oceanic sediment calculation')
        calculate_sediment_ocean_interaction(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,topo_time_grid,Output_way)

def plot_model_RSL(plot,Input_way,Output_way):
    load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}/LOAD/TOTAL_LOAD'))
    geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}/LOAD/TOTAL_GEOID'))
    dESL=calculate_dESL(Input_way)
    load_time_grid.height_time_coeff=(geoid_time_grid.viscuous_deformation+geoid_time_grid.elastic_deformation-load_time_grid.viscuous_deformation-load_time_grid.elastic_deformation)*dESL[2:,np.newaxis].repeat(load_time_grid.viscuous_deformation.shape[1],axis=1)
    time_grid=load_time_grid
    if plot['plot']:
        if len(plot['frames'])>0 :
            fig,ax=plt.subplots(len(plot['frames'])*len(plot['times']),1,figsize=(29.7,10*len(plot['frames'])*len(plot['times'])),subplot_kw={'projection': ccrs.PlateCarree(), "aspect": 1})
            k=0
            for i_time in range(len(plot['times'])):
                for i_frame in range(len(plot['frames'])):
                    axes=ax[k]
                    time=plot['times'][i_time]
                    central_longitude=0
                    #ax_total  = plt.subplot(subplot_size+i_time*3+1, projection=ccrs.PlateCarree(central_longitude=central_longitude))
                    ax_total=plot_load_frame(axes,time_grid,plot['frames'][i_frame],time,plot['frames_resolution'][i_frame],name=f'RSL',vmin=plot['frames_min_max'][i_frame][0],vmax=plot['frames_min_max'][i_frame][1])
                    k+=1
        plt.close(fig)
        fig.savefig(f'{Output_way}/RSL.pdf')

def plot_model_result_map(Input_way,plot):
    '''
    The _`plot_model_result_map` function run the plot functions for all type of load. This function use the `plot_frame`_ function, see this function for more details. 

    Attribute :
    ----------- 
        input_way_sed : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        plot : dict(plot=True,times=[n_t],frames=[n_f*(lon1,lon1,lat1,lat2)],frames_resolution=[n_f],frames_min_max=np.array([[[min_load],max_load,min_geoid,max_geoid]]),contours_v=[n_f*2*[contours...]],transects=[n_t*(lat1,lon1,lat2,lon2)],point_density=[n_t], transect_min_max=[n_t(min_load,max_load,min_geoid,max_geoid)], points=np.array([n_p[lat,lon]]))
            The plot dictionnary used to define all the plot parameters of the output. In this dictionary, there is three main group of parameters, the frame plot parameters, the transect plot parameters and the points plot parameters. The frame plot parameters contain the frame locations in frames, the frames resolution in frames_resolution, the frames min and max value to plot in frames_min_max and the contours value to plot in contours_v. The transect plot parameters contains the location of the begining and end of the transects in transects, the point density along the transect in point_density and the min max value of the transect in transect_min_max. The points plot parameters contains the locations of the differents plots in points. 

    Returns : 
    --------- 
        None 

    '''
    Input_way=Input_way+'/LOAD'
    plot_frame(Input_way,plot,Input_way,'ICE')
    plot_frame(Input_way,plot,Input_way,'SEDIMENT')
    plot_frame(Input_way,plot,Input_way,'OCEAN')
    plot_frame(Input_way,plot,Input_way,'NO_SEDIMENT')
    plot_frame(Input_way,plot,Input_way,'TOTAL')

def calc_transect(Input_way,time,transect,point_density,type,type_bis):
    '''
    The _`calc_transect` function evaluate the spherical harmonics coefficient function run the plot functions for all type of load. This function use the `plot_frame`_ function, see this function for more details. 

    Attribute :
    ----------- 
        input_way_sed : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        time : float
            The selected time to plot the transect.
        transect : (lat1,lon1,lat2,lon2)
            The coordinate of the beginning and the end of the transect.
        point_dentity : int
            Number of points along the transect, this will define the resolution
        type : str
            This define the type of load considered. Must be "SEDIMENT", "NO_SEDIMENT", "OCEAN", "ICE" or "TOTAL".
        type_bis : str
            This define if the plot concern the geoid of the ground deformation. Must be "LOAD" or "GEOID".

    Returns : 
    --------- 
        None 
        
    '''
    
    if type is 'OCEAN' :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_{type_bis}'))
        t_it=np.where(load_time_grid.time_step==time)[0][0]
        load_time_grid.coeff=load_time_grid.viscuous_deformation[t_it-2,:]-load_time_grid.viscuous_deformation[t_it-3,:]+load_time_grid.elastic_deformation[t_it-2,:]-load_time_grid.elastic_deformation[t_it-3,:]
        load_time_grid.coeff=load_time_grid.coeff/(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
        tr_ocean=load_time_grid.along_transect(coord=transect,point_density=point_density)
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_{type_bis}'))
        load_time_grid.coeff=load_time_grid.viscuous_deformation[t_it-2,:]-load_time_grid.viscuous_deformation[t_it-3,:]+load_time_grid.elastic_deformation[t_it-2,:]-load_time_grid.elastic_deformation[t_it-3,:]
        load_time_grid.coeff=load_time_grid.coeff/(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
        tr_oceanic_sediment=load_time_grid.along_transect(coord=transect,point_density=point_density)
        tr=tr_ocean+tr_oceanic_sediment
    elif type is 'SEDIMENT' :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_{type_bis}'))
        t_it=np.where(load_time_grid.time_step==time)[0][0]
        load_time_grid.coeff=load_time_grid.viscuous_deformation[t_it-2,:]-load_time_grid.viscuous_deformation[t_it-3,:]+load_time_grid.elastic_deformation[t_it-2,:]-load_time_grid.elastic_deformation[t_it-3,:]
        load_time_grid.coeff=load_time_grid.coeff/(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
        tr_sediment=load_time_grid.along_transect(coord=transect,point_density=point_density)
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_{type_bis}'))
        load_time_grid.coeff=load_time_grid.viscuous_deformation[t_it-2,:]-load_time_grid.viscuous_deformation[t_it-3,:]+load_time_grid.elastic_deformation[t_it-2,:]-load_time_grid.elastic_deformation[t_it-3,:]
        load_time_grid.coeff=load_time_grid.coeff/(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
        tr_oceanic_sediment=load_time_grid.along_transect(coord=transect,point_density=point_density)
        tr=tr_sediment-tr_oceanic_sediment
    elif type is 'NO_SEDIMENT' :
        ice_load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\ICE_{type_bis}'))
        ocean_load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEAN_{type_bis}'))
        t_it=np.where(ice_load_time_grid.time_step==time)[0][0]

        ice_load_time_grid.coeff=ice_load_time_grid.viscuous_deformation[t_it-2,:]-ice_load_time_grid.viscuous_deformation[t_it-3,:]+ice_load_time_grid.elastic_deformation[t_it-2,:]-ice_load_time_grid.elastic_deformation[t_it-3,:]
        ice_load_time_grid.coeff=ice_load_time_grid.coeff/(ice_load_time_grid.time_step[t_it-1]-ice_load_time_grid.time_step[t_it])

        ocean_load_time_grid.coeff=ocean_load_time_grid.viscuous_deformation[t_it-2,:]-ocean_load_time_grid.viscuous_deformation[t_it-3,:]+ocean_load_time_grid.elastic_deformation[t_it-2,:]-ocean_load_time_grid.elastic_deformation[t_it-3,:]
        ocean_load_time_grid.coeff=ocean_load_time_grid.coeff/(ocean_load_time_grid.time_step[t_it-1]-ocean_load_time_grid.time_step[t_it])

        tr_ice=ice_load_time_grid.along_transect(coord=transect,point_density=point_density)
        tr_ocean=ocean_load_time_grid.along_transect(coord=transect,point_density=point_density)

        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_{type_bis}'))
        load_time_grid.coeff=load_time_grid.viscuous_deformation[t_it-2,:]-load_time_grid.viscuous_deformation[t_it-3,:]+load_time_grid.elastic_deformation[t_it-2,:]-load_time_grid.elastic_deformation[t_it-3,:]
        load_time_grid.coeff=load_time_grid.coeff/(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
        tr_oceanic_sediment=load_time_grid.along_transect(coord=transect,point_density=point_density)
        tr=tr_ice+tr_ocean+tr_oceanic_sediment
    else :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_{type_bis}'))
        t_it=np.where(load_time_grid.time_step==time)[0][0]
        load_time_grid.coeff=load_time_grid.viscuous_deformation[t_it-2,:]-load_time_grid.viscuous_deformation[t_it-3,:]+load_time_grid.elastic_deformation[t_it-2,:]-load_time_grid.elastic_deformation[t_it-3,:]
        load_time_grid.coeff=load_time_grid.coeff/(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
        tr=load_time_grid.along_transect(coord=transect,point_density=point_density)
    dt=(load_time_grid.time_step[t_it-1]-load_time_grid.time_step[t_it])
    return tr,dt


def plot_model_result_cross_section(Input_way,plot):
    '''
    The _`plot_model_result_cross_section` function evaluate the spherical harmonics coefficient function and plot the results along a defined transect in plot.

    Attribute :
    ----------- 
        Input_way_sed : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        plot : dict(plot=True,times=[n_t],frames=[n_f*(lon1,lon1,lat1,lat2)],frames_resolution=[n_f],frames_min_max=np.array([[[min_load],max_load,min_geoid,max_geoid]]),contours_v=[n_f*2*[contours...]],transects=[n_t*(lat1,lon1,lat2,lon2)],point_density=[n_t], transect_min_max=[n_t(min_load,max_load,min_geoid,max_geoid)], points=np.array([n_p[lat,lon]]))
            The plot dictionnary used to define all the plot parameters of the output. In this dictionary, there is three main group of parameters, the frame plot parameters, the transect plot parameters and the points plot parameters. The frame plot parameters contain the frame locations in frames, the frames resolution in frames_resolution, the frames min and max value to plot in frames_min_max and the contours value to plot in contours_v. The transect plot parameters contains the location of the begining and end of the transects in transects, the point density along the transect in point_density and the min max value of the transect in transect_min_max. The points plot parameters contains the locations of the differents plots in points. 

    Returns : 
    --------- 
        None 
        
    '''
    Input_way=Input_way+'/LOAD'
    sediment_color=(0.4,0.7,0.5)
    sediment_color_dark=(0.2,0.4,0.25)
    ice_color=(0.1,0.8,0.8)
    ocean_color=(0.2,0.2,0.6)
    glaciohydrostatic_color=(0.1,0.1,0.4)
    if len(plot['transects'])>0 :
        fig,ax=plt.subplots(len(plot['transects'])*len(plot['times']),2,figsize=(29.7,5*len(plot['transects'])*len(plot['times'])))
        fig_loc,ax_map=plt.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree(), "aspect": 1})
        alpha_ocean=0
        coast_line_width=0.2
        ax_map.set_global()
        cartopy.mpl.geoaxes.GeoAxes.gridlines(ax_map,crs=ccrs.PlateCarree(),draw_labels=True)
        ax_map.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
        ax_map.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)
        for i_transect in range(len(plot['transects'])):
            ax_map.plot((plot['transects'][i_transect][1],plot['transects'][i_transect][3]),(plot['transects'][i_transect][0],plot['transects'][i_transect][2]),color='r',linewidth=0.5)
            ax_map.annotate(f'{i_transect})', # this is the text
                 (plot['transects'][i_transect][1],plot['transects'][i_transect][0]), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,1), # distance from text to points (x,y)
                 ha='center',
                 fontsize=3)
        plt.close(fig_loc)
        fig_loc.savefig(f'{Input_way}/transect_localisation.pdf')
        k=0
        for i_time in range(len(plot['times'])):
            for i_transect in range(len(plot['transects'])):
                coord=plot['transects'][i_transect]
                theta=np.linspace(coord[0],coord[2],plot['point_density'][i_transect])/180*2*np.pi
                phi=np.linspace(coord[1],coord[3],plot['point_density'][i_transect])/180*2*np.pi
                phiA=np.append(phi,0)
                phiB=np.append(0,phi)
                D= np.append(0,6371*np.arccos(np.sin(phiA)*np.sin(phiB)+np.cos(phiA)*np.cos(phiB)*np.cos(np.append(0,np.append(np.diff(theta),0))))[1:-1].cumsum())/2
                axes=ax[k,:]
                type_bis='LOAD'
                tr_ice_load,dt=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'ICE',type_bis)
                tr_sed_load,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'SEDIMENT',type_bis)
                tr_oc_load,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'OCEAN',type_bis)
                tr_nosed_load,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'NO_SEDIMENT',type_bis)
                tr_tot_load,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'TOTAL',type_bis)
                axes[0].plot(D,tr_ice_load,label='ice_load',color=ice_color)
                axes[0].plot(D,tr_sed_load,label='sediment_load',color=sediment_color)
                axes[0].plot(D,tr_oc_load,label='ocean_load',color=ocean_color)
                axes[0].plot(D,tr_nosed_load,label='glaciohydrostatic_load',color=glaciohydrostatic_color)
                axes[0].plot(D,tr_tot_load,label='total load',color='k')
                axes[0].set_ylim(plot['transect_min_max'][i_transect][0],plot['transect_min_max'][i_transect][1])
                axes[0].set_xlabel(f'distance (km)')
                axes[0].set_ylabel(f'Vertical Land motion averaged over {dt*1000} year (mm/yr)')
                t=plot['times'][i_time]
                axes[0].grid()
                axes[0].set_title(f'VLM estimated at {t} kyr along transect {i_transect})')
                axes[0].legend()
                type_bis='GEOID'
                tr_ice_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'ICE',type_bis)
                tr_sed_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'SEDIMENT',type_bis)
                tr_oc_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'OCEAN',type_bis)
                tr_nosed_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'NO_SEDIMENT',type_bis)
                tr_tot_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'TOTAL',type_bis)
                axes[1].plot(D,tr_ice_geoid,label='ice_geoid',color=ice_color)
                axes[1].plot(D,tr_sed_geoid,label='sediment_geoid',color=sediment_color)
                axes[1].plot(D,tr_oc_geoid,label='ocean_geoid',color=ocean_color)
                axes[1].plot(D,tr_nosed_geoid,label='glaciohydrostatic_geoid',color=glaciohydrostatic_color)
                axes[1].plot(D,tr_tot_geoid,label='total geoid',color='k')
                axes[1].set_xlabel(f'distance (km)')
                axes[1].set_ylabel(f'Vertical Geoid motion averaged over {dt*1000} year (mm/yr)')
                t=plot['times'][i_time]
                axes[1].set_title(f'VGM estimated at {t} kyr along transect {i_transect})')
                axes[1].grid()
                axes[1].set_ylim(plot['transect_min_max'][i_transect][2],plot['transect_min_max'][i_transect][3])
                axes[1].legend()
                k+=1
        plt.close(fig)
        fig.savefig(f'{Input_way}/transect_deformation.pdf')

    

def calc_point(Input_way,points,type,type_bis):
    '''
    The _`calc_point` function evaluate the spherical harmonics coefficient for a type of model output at precise points locations.

    Attribute :
    ----------- 
        Input_way : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        time : float
            The selected time to plot the transect.
        points : np.array(n_p*[lat,lon])
            The locations of the different points you wan't to extract data.
        type : str
            This define the type of load considered. Must be "SEDIMENT", "NO_SEDIMENT", "OCEAN", "ICE" or "TOTAL".
        type_bis : str
            This define if the plot concern the geoid of the ground deformation. Must be "LOAD" or "GEOID".

    Returns : 
    --------- 
        None 
        
    '''
    
    if type is 'OCEAN' :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_ocean=load_time_grid.point_time(points)

        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_oceanic_sediment=load_time_grid.point_time(points)

        points=points_ocean+points_oceanic_sediment
    elif type is 'SEDIMENT' :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_sediment=load_time_grid.point_time(points)

        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_oceanic_sediment=load_time_grid.point_time(points)

        points=points_sediment-points_oceanic_sediment

    elif type is 'NO_SEDIMENT' :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEAN_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_ocean=load_time_grid.point_time(points)

        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\ICE_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_ice=load_time_grid.point_time(points)

        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points_oceanic_sediment=load_time_grid.point_time(points)

        points=points_ocean+points_oceanic_sediment+points_ice
    else :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_{type_bis}'))
        load_time_grid.height_time_coeff=load_time_grid.viscuous_deformation+load_time_grid.elastic_deformation
        points=load_time_grid.point_time(points)

    return points,load_time_grid.time_step

def plot_model_output_points(Input_way,plot):
    '''
    The _`plot_model_output_points` function plot the output of the model at different points.
    
    Attribute :
    ----------- 
        Input_way : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        plot : dict(plot=True,times=[n_t],frames=[n_f*(lon1,lon1,lat1,lat2)],frames_resolution=[n_f],frames_min_max=np.array([[[min_load],max_load,min_geoid,max_geoid]]),contours_v=[n_f*2*[contours...]],transects=[n_t*(lat1,lon1,lat2,lon2)],point_density=[n_t], transect_min_max=[n_t(min_load,max_load,min_geoid,max_geoid)], points=np.array([n_p[lat,lon]]))
            The plot dictionnary used to define all the plot parameters of the output. In this dictionary, there is three main group of parameters, the frame plot parameters, the transect plot parameters and the points plot parameters. The frame plot parameters contain the frame locations in frames, the frames resolution in frames_resolution, the frames min and max value to plot in frames_min_max and the contours value to plot in contours_v. The transect plot parameters contains the location of the begining and end of the transects in transects, the point density along the transect in point_density and the min max value of the transect in transect_min_max. The points plot parameters contains the locations of the differents plots in points. 
        

    Returns : 
    --------- 
        None 
        
    '''
    sediment_color=(0.4,0.7,0.5)
    sediment_color_dark=(0.2,0.4,0.25)
    ice_color=(0.1,0.8,0.8)
    ocean_color=(0.2,0.2,0.6)
    glaciohydrostatic_color=(0.7,0.2,0.2)
    if len(plot['transects'])>0 :
        fig,axes=plt.subplots(len(plot['points']),3,figsize=(29.7,4*len(plot['points'])))
        fig_loc,ax_map=plt.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree(), "aspect": 1})
        alpha_ocean=0
        coast_line_width=0.2
        ax_map.set_global()
        cartopy.mpl.geoaxes.GeoAxes.gridlines(ax_map,crs=ccrs.PlateCarree(),draw_labels=True)
        ax_map.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
        ax_map.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)
        for i_point in range(len(plot['points'])):
            ax_map.scatter((plot['points'][i_point,1],plot['points'][i_point,1]),(plot['points'][i_point,0],plot['points'][i_point,0]),color='r',s=0.05)
            ax_map.annotate(f'{i_point})', # this is the text
                 (plot['points'][i_point,1],plot['points'][i_point,0]), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,0.2), # distance from text to points (x,y)
                 ha='center',
                 fontsize=1)
        plt.close(fig_loc)
        fig_loc.savefig(f'{Input_way}/LOAD/Points_localisation.pdf')
        coord=plot['points']
        type_bis='LOAD'
        points_ice_load,time_step=calc_point(Input_way+'/LOAD',plot['points'],'ICE',type_bis)
        points_sed_load,_=calc_point(Input_way+'/LOAD',plot['points'],'SEDIMENT',type_bis)
        points_oc_load,_=calc_point(Input_way+'/LOAD',plot['points'],'OCEAN',type_bis)
        points_nosed_load,_=calc_point(Input_way+'/LOAD',plot['points'],'NO_SEDIMENT',type_bis)
        points_tot_load,_=calc_point(Input_way+'/LOAD',plot['points'],'TOTAL',type_bis)
        for i_point in range(len(plot['points'])):
            axes[i_point,0].plot(time_step,points_ice_load[i_point,:],label='ice_load',color=ice_color)
            axes[i_point,0].plot(time_step,points_sed_load[i_point,:],label='sediment_load',color=sediment_color)
            axes[i_point,0].plot(time_step,points_oc_load[i_point,:],label='ocean_load',color=ocean_color)
            axes[i_point,0].plot(time_step,points_nosed_load[i_point,:],label='glaciohydrostatic_load',color=glaciohydrostatic_color)
            axes[i_point,0].plot(time_step,points_tot_load[i_point,:],label='total load',color='k')
            #axes[i_point,0].set_ylim(plot['transect_min_max'][i_transect][0],plot['transect_min_max'][i_transect][1])
            axes[i_point,0].set_xlabel(f'time (kyr)')
            axes[i_point,0].set_ylabel(f'Vertical Land motion (mm/yr)')
            axes[i_point,0].grid()
            axes[i_point,0].set_title(f'VLM at point {i_point})')
            axes[i_point,0].legend()
        type_bis='GEOID'
        points_ice_geoid,time_step=calc_point(Input_way+'/LOAD',plot['points'],'ICE',type_bis)
        points_sed_geoid,_=calc_point(Input_way+'/LOAD',plot['points'],'SEDIMENT',type_bis)
        points_oc_geoid,_=calc_point(Input_way+'/LOAD',plot['points'],'OCEAN',type_bis)
        points_nosed_geoid,_=calc_point(Input_way+'/LOAD',plot['points'],'NO_SEDIMENT',type_bis)
        points_tot_geoid,_=calc_point(Input_way+'/LOAD',plot['points'],'TOTAL',type_bis)
        for i_point in range(len(plot['points'])):
            axes[i_point,1].plot(time_step,points_ice_geoid[i_point,:],label='ice_geoid',color=ice_color)
            axes[i_point,1].plot(time_step,points_sed_geoid[i_point,:],label='sediment_geoid',color=sediment_color)
            axes[i_point,1].plot(time_step,points_oc_geoid[i_point,:],label='ocean_geoid',color=ocean_color)
            axes[i_point,1].plot(time_step,points_nosed_geoid[i_point,:],label='glaciohydrostatic_geoid',color=glaciohydrostatic_color)
            axes[i_point,1].plot(time_step,points_tot_geoid[i_point,:],label='total_geoid',color='k')
            #axes[i_point,0].set_ylim(plot['transect_min_max'][i_transect][0],plot['transect_min_max'][i_transect][1])
            axes[i_point,1].set_xlabel(f'time (kyr)')
            axes[i_point,1].set_ylabel(f'Vertical Geoid motion (mm/yr)')
            axes[i_point,1].grid()
            axes[i_point,1].set_title(f'VGM at point {i_point})')
            axes[i_point,1].legend()
        dESL=calculate_dESL(Input_way)
        for i_point in range(len(plot['points'])):
            axes[i_point,2].plot(time_step,dESL[:].cumsum()-dESL[:].cumsum()[-1],label='ESL',linestyle='dashed',color='k')
            axes[i_point,2].plot(time_step,(dESL[:]+points_tot_geoid[i_point,:]-points_tot_load[i_point,:]).cumsum()-(dESL[:]+points_tot_geoid[i_point,:]-points_tot_load[i_point,:]).cumsum()[-1],label='RSL',color='k')
            axes[i_point,2].legend()
            axes[i_point,2].grid()
            axes[i_point,2].set_xlabel(f'time (kyr)')
            axes[i_point,2].set_ylabel(f'Depth below present sea level (mbpsl)')
            axes[i_point,2].set_title(f'RSL and ESL at point {i_point})')

        # type_bis='GEOID'
        # tr_ice_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'ICE',type_bis)
        # tr_sed_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'SEDIMENT',type_bis)
        # tr_oc_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'OCEAN',type_bis)
        # tr_tot_geoid,_=calc_transect(Input_way,plot['times'][i_time],plot['transects'][i_transect],plot['point_density'][i_transect],'TOTAL',type_bis)
        # axes[1].plot(D,tr_ice_geoid,label='ice_load',color=ice_color)
        # axes[1].plot(D,tr_sed_geoid,label='sediment_load',color=sediment_color)
        # axes[1].plot(D,tr_oc_geoid,label='ocean_load',color=ocean_color)
        # axes[1].plot(D,tr_tot_geoid,label='total load',color='k')
        # axes[1].set_xlabel(f'distance (km)')
        # axes[1].set_ylabel(f'Vertical Geoid motion averaged over {dt*1000} year (mm/yr)')
        # t=plot['times'][i_time]
        # axes[1].set_title(f'VGM estimated at {t} kyr along transect {i_transect})')
        # axes[1].grid()
        # axes[1].set_ylim(plot['transect_min_max'][i_transect][2],plot['transect_min_max'][i_transect][3])
        # axes[1].legend()
        plt.close(fig)
        fig.savefig(f'{Input_way}/LOAD/Points_deformation.pdf')

def plot_frame(Input_way,plot,Output_way,type):
    '''
    The _`plot_frame` function define the spherical harmonic coefficient to be expanded into a grid and plot in the frame defined in plot. This function make the difference between sediment subsidence and true sediment subsidence where water replaced by sediment is included. 
    
    Attribute :
    ----------- 
        Input_way : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        plot : dict(plot=True,times=[n_t],frames=[n_f*(lon1,lon1,lat1,lat2)],frames_resolution=[n_f],frames_min_max=np.array([[[min_load],max_load,min_geoid,max_geoid]]),contours_v=[n_f*2*[contours...]],transects=[n_t*(lat1,lon1,lat2,lon2)],point_density=[n_t], transect_min_max=[n_t(min_load,max_load,min_geoid,max_geoid)], points=np.array([n_p[lat,lon]]))
            The plot dictionnary used to define all the plot parameters of the output. In this dictionary, there is three main group of parameters, the frame plot parameters, the transect plot parameters and the points plot parameters. The frame plot parameters contain the frame locations in frames, the frames resolution in frames_resolution, the frames min and max value to plot in frames_min_max and the contours value to plot in contours_v. The transect plot parameters contains the location of the begining and end of the transects in transects, the point density along the transect in point_density and the min max value of the transect in transect_min_max. The points plot parameters contains the locations of the differents plots in points. 
        Output_way : str
            way where the plot will be saved.
        type : str
            This define the type of load considered. Must be "SEDIMENT", "NO_SEDIMENT", "OCEAN", "ICE" or "TOTAL".

        

    Returns : 
    --------- 
        None 
        
    '''
    if type is "SEDIMENT" :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_LOAD'))
        sediment_load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_LOAD'))
        load_time_grid.viscuous_deformation=sediment_load_time_grid.viscuous_deformation-load_time_grid.viscuous_deformation
        load_time_grid.elastic_deformation=sediment_load_time_grid.elastic_deformation-load_time_grid.elastic_deformation
        geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_GEOID'))
        sediment_geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_GEOID'))
        geoid_time_grid.viscuous_deformation=sediment_geoid_time_grid.viscuous_deformation-geoid_time_grid.viscuous_deformation
        geoid_time_grid.elastic_deformation=sediment_geoid_time_grid.elastic_deformation-geoid_time_grid.elastic_deformation
    if type is "OCEAN" :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_LOAD'))
        sediment_load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_LOAD'))
        load_time_grid.viscuous_deformation=sediment_load_time_grid.viscuous_deformation+load_time_grid.viscuous_deformation
        load_time_grid.elastic_deformation=sediment_load_time_grid.elastic_deformation+load_time_grid.elastic_deformation
        geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_GEOID'))
        sediment_geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_GEOID'))
        geoid_time_grid.viscuous_deformation=sediment_geoid_time_grid.viscuous_deformation+geoid_time_grid.viscuous_deformation
        geoid_time_grid.elastic_deformation=sediment_geoid_time_grid.elastic_deformation+geoid_time_grid.elastic_deformation
    elif type is "NO_SEDIMENT":
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\TOTAL_LOAD'))
        sediment_load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\SEDIMENT_LOAD'))
        oceanic_sediment_load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_LOAD'))
        load_time_grid.viscuous_deformation=-sediment_load_time_grid.viscuous_deformation+load_time_grid.viscuous_deformation+oceanic_sediment_load_time_grid.viscuous_deformation
        load_time_grid.elastic_deformation=-sediment_load_time_grid.elastic_deformation+load_time_grid.elastic_deformation+oceanic_sediment_load_time_grid.elastic_deformation
        
        ocean_geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEAN_GEOID'))
        geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\ICE_GEOID'))
        oceanic_sediment_geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\OCEANIC_SEDIMENT_GEOID'))
        geoid_time_grid.viscuous_deformation=ocean_geoid_time_grid.viscuous_deformation+geoid_time_grid.viscuous_deformation+oceanic_sediment_geoid_time_grid.viscuous_deformation
        geoid_time_grid.elastic_deformation=ocean_geoid_time_grid.elastic_deformation+geoid_time_grid.elastic_deformation+oceanic_sediment_geoid_time_grid.elastic_deformation
    else :
        load_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_LOAD'))
        geoid_time_grid=LOAD_TIME_GRID(from_file=(True,f'{Input_way}\{type}_GEOID'))
    plot_model_output(load_time_grid,plot,Output_way,f'load_{type}',vmin=(plot['frames_min_max'][:,0]),vmax=plot['frames_min_max'][:,1],contours_v=plot['contours_v'][0])
    plot_model_output(geoid_time_grid,plot,Output_way,f'geoid_{type}',vmin=(plot['frames_min_max'][:,2]),vmax=plot['frames_min_max'][:,3],contours_v=plot['contours_v'][1])

def plot_model_output(time_grid,plot,Output_way,type=None,vmin=None,vmax=None,contours_v=None):
    '''
    The _`plot_model_output` function plot the frames for each frames and at each time steps. 
    
    Attribute :
    ----------- 
        time_grid : TIME_GRID from `ref`:{TIME_GRID <TIME_GRID>}
        plot : dict(plot=True,times=[n_t],frames=[n_f*(lon1,lon1,lat1,lat2)],frames_resolution=[n_f],frames_min_max=np.array([[[min_load],max_load,min_geoid,max_geoid]]),contours_v=[n_f*2*[contours...]],transects=[n_t*(lat1,lon1,lat2,lon2)],point_density=[n_t], transect_min_max=[n_t(min_load,max_load,min_geoid,max_geoid)], points=np.array([n_p[lat,lon]]))
            The plot dictionnary used to define all the plot parameters of the output. In this dictionary, there is three main group of parameters, the frame plot parameters, the transect plot parameters and the points plot parameters. The frame plot parameters contain the frame locations in frames, the frames resolution in frames_resolution, the frames min and max value to plot in frames_min_max and the contours value to plot in contours_v. The transect plot parameters contains the location of the begining and end of the transects in transects, the point density along the transect in point_density and the min max value of the transect in transect_min_max. The points plot parameters contains the locations of the differents plots in points. 
        Output_way : str
            way where the plot will be saved.
        type : str
            This define the type of load considered. Must be "SEDIMENT", "NO_SEDIMENT", "OCEAN", "ICE" or "TOTAL".
        vmin : float 
            The minimal value for the colormap of your plot
        vmax : float
            The maximal value for the colormap of your plot
        contours_v : list
            List of the contours values to be plot

    Returns : 
    --------- 
        None 
        
    '''
    if plot['plot']:
        if len(plot['frames'])>0 :
            fig,ax=plt.subplots(len(plot['frames'])*len(plot['times']),3,figsize=(29.7,10*len(plot['frames'])*len(plot['times'])),subplot_kw={'projection': ccrs.PlateCarree(), "aspect": 1})
            k=0
            for i_time in range(len(plot['times'])):
                for i_frame in range(len(plot['frames'])):
                    axes=ax[k,:]
                    time=plot['times'][i_time]
                    central_longitude=0
                    #ax_total  = plt.subplot(subplot_size+i_time*3+1, projection=ccrs.PlateCarree(central_longitude=central_longitude))
                    time_grid.height_time_coeff=time_grid.viscuous_deformation+time_grid.elastic_deformation
                    ax_total=plot_load_frame(axes[0],time_grid,plot['frames'][i_frame],time,plot['frames_resolution'][i_frame],name=f'total {type}',vmin=vmin[i_frame],vmax=vmax[i_frame],contour=contours_v[i_frame])

                    #ax_viscous  = plt.subplot(subplot_size+i_time*3+2, projection=ccrs.PlateCarree(central_longitude=central_longitude))
                    time_grid.height_time_coeff=time_grid.viscuous_deformation
                    ax_viscous=plot_load_frame(axes[1],time_grid,plot['frames'][i_frame],time,plot['frames_resolution'][i_frame],name=f'viscous {type}',vmin=vmin[i_frame],vmax=vmax[i_frame],contour=contours_v[i_frame])

                    #ax_elastic = plt.subplot(subplot_size+i_time*3+3, projection=ccrs.PlateCarree(central_longitude=central_longitude))
                    time_grid.height_time_coeff=time_grid.elastic_deformation
                    ax_elastic=plot_load_frame(axes[2],time_grid,plot['frames'][i_frame],time,plot['frames_resolution'][i_frame],name=f'elastic {type}',vmin=vmin[i_frame],vmax=vmax[i_frame],contour=contours_v[i_frame]) 
                    k+=1
        plt.close(fig)
        fig.savefig(f'{Output_way}/{type}_deformation.pdf')

import cartopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib import cm
import matplotlib.colors as colors

def plot_precomputation(initial_time_step,initial_grid,precomputed_grid,selected_time=1,area=None,save_way='',lon_init=None,lat_init=None,topo=False,vmin=None,vmax=None,vcenter=0,cmap=cm.get_cmap('bwr', 255)):
    '''
    The _`plot_precomputation` function plot the output of the precomputation as the ice, sediment and totpography. These plot allow you to compare the entry and the output to be sure there is no errors.
    
    Attribute :
    ----------- 
        Input_way : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        plot : dict(plot=True,times=[n_t],frames=[n_f*(lon1,lon1,lat1,lat2)],frames_resolution=[n_f],frames_min_max=np.array([[[min_load],max_load,min_geoid,max_geoid]]),contours_v=[n_f*2*[contours...]],transects=[n_t*(lat1,lon1,lat2,lon2)],point_density=[n_t], transect_min_max=[n_t(min_load,max_load,min_geoid,max_geoid)], points=np.array([n_p[lat,lon]]))
            The plot dictionnary used to define all the plot parameters of the output. In this dictionary, there is three main group of parameters, the frame plot parameters, the transect plot parameters and the points plot parameters. The frame plot parameters contain the frame locations in frames, the frames resolution in frames_resolution, the frames min and max value to plot in frames_min_max and the contours value to plot in contours_v. The transect plot parameters contains the location of the begining and end of the transects in transects, the point density along the transect in point_density and the min max value of the transect in transect_min_max. The points plot parameters contains the locations of the differents plots in points. 
        Output_way : str
            way where the plot will be saved.
        type : str
            This define the type of load considered. Must be "SEDIMENT", "NO_SEDIMENT", "OCEAN", "ICE" or "TOTAL".
        vmin : float 
            The minimal value for the colormap of your plot
        vmax : float
            The maximal value for the colormap of your plot

    Returns : 
    --------- 
        None 
        
    '''
    our_ind=np.where((precomputed_grid.time_step-selected_time)==0)[0][0]
    in_ind=np.where((initial_time_step-selected_time)==0)[0][0]
    alpha_ocean=0
    coast_line_width=0.5
    norm1 = colors.TwoSlopeNorm(vmin=-0.1,vmax=None,vcenter=0)
    fig = plt.figure(figsize=(12, 6), facecolor="none")

    # plot the initial data

    norm = colors.TwoSlopeNorm(vmin=vmin,vmax=vmax,vcenter=vcenter)

    if area is None :
        ax_init  = plt.subplot(121, projection=ccrs.Mollweide())
        ax_init.set_global()
    else :
        ax_init  = plt.subplot(121, projection=ccrs.PlateCarree())
        ax_init.set_extent(area)

    if not(lon_init is None) :
        m_init = ax_init.scatter(lon_init,lat_init,c=initial_grid[:in_ind+1,:].sum(0),transform=ccrs.PlateCarree(),zorder=0,norm=norm,cmap=cmap)
        cbar_init=plt.colorbar(mappable=m_init, orientation="horizontal")
        cbar_init.set_label(precomputed_grid.time_grid_name + '(m)')
        #cbar_init.set_ticks([0, initial_grid[:in_ind+1].sum(0).max()])
        ax_init.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
        ax_init.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)
        cartopy.mpl.geoaxes.GeoAxes.gridlines(ax_init,crs=ccrs.PlateCarree(),draw_labels=True)
    else : 
        if topo : 
            m_init = ax_init.imshow(initial_grid, origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, cmap=cmap, interpolation="gaussian",norm=norm)
            cbar_init=plt.colorbar(mappable=m_init, orientation="horizontal")
            cbar_init.set_label(precomputed_grid.time_grid_name + '(m)')
            #cbar_init.set_ticks([0, initial_grid[:in_ind+1].sum(0).max()])
            ax_init.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
            ax_init.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)
        else : 
            m_init = ax_init.imshow(initial_grid[:in_ind+1].sum(0), origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, cmap=cmap, interpolation="gaussian",norm=norm)
            cbar_init=plt.colorbar(mappable=m_init, orientation="horizontal")
            cbar_init.set_label(precomputed_grid.time_grid_name + '(m)')
            #cbar_init.set_ticks([0, initial_grid[:in_ind+1].sum(0).max()])
            ax_init.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
            ax_init.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)
        cartopy.mpl.geoaxes.GeoAxes.gridlines(ax_init,crs=ccrs.PlateCarree(),draw_labels=True)

    if area is None :
        ax_precomp  = plt.subplot(122, projection=ccrs.Mollweide())
        ax_precomp.set_global()
    else :
        ax_precomp  = plt.subplot(122, projection=ccrs.PlateCarree())
        ax_precomp.set_extent(area)

    colormap = cmap
    if topo :
        m2 = ax_precomp.imshow(precomputed_grid.height_time_grid[-1,:,:], origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, cmap=colormap, interpolation="gaussian",norm=norm)
        cbar2=plt.colorbar(mappable=m2, orientation="horizontal")
        cbar2.set_label(precomputed_grid.time_grid_name + ' interpolated (m)')
        #cbar2.set_ticks([0, precomputed_grid.height_time_grid[:our_ind+1,:,:].sum(0).max()])
        ax_precomp.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
        ax_precomp.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)

    else : 
        m2 = ax_precomp.imshow(precomputed_grid.height_time_grid[:our_ind+1,:,:].sum(0), origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, cmap=colormap, interpolation="gaussian",norm=norm)
        cbar2=plt.colorbar(mappable=m2, orientation="horizontal")
        cbar2.set_label(precomputed_grid.time_grid_name + ' interpolated (m)')
        #cbar2.set_ticks([0, precomputed_grid.height_time_grid[:our_ind+1,:,:].sum(0).max()])
        ax_precomp.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
        ax_precomp.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)
    cartopy.mpl.geoaxes.GeoAxes.gridlines(ax_precomp,crs=ccrs.PlateCarree(),draw_labels=True)

    plt.close(fig)

    fig.savefig(save_way+'/'+precomputed_grid.time_grid_name+'.pdf')

def calculate_dESL(Input_way):
    '''
    The _`calculate_dESL` function compute the ESL variation over time including a varaible ocean surface.
    
    Attribute :
    ----------- 
        Input_way : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.

    Returns : 
    --------- 
        dESL : np.array([time_step_number,])
            The variation of the ESL at each time step.
        
    '''
    ocean_grid=OCEAN_TIME_GRID(from_file=(True,Input_way+"/OCE"))

    ice=ICE_TIME_GRID(from_file=(True,Input_way+"/ICE_ICE6G"))

    topo=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+"/topo_SL"))
    dESL=np.zeros(ice.time_step_number)
    for t_it in range(ice.time_step_number-1):
        ocean_grid.evaluate_ocean(topo.height_time_grid[t_it,:,:]).grdtocoeff()
        dESL[t_it]=np.real(1/ocean_grid.coeff[0] * (- ice.rho/ocean_grid.rho*ice.height_time_coeff[t_it,0]))
    return dESL

import cartopy
import cartopy.crs as ccrs
import matplotlib
from matplotlib import cm
import matplotlib as mpl
import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt

def plot_load_frame(ax,data,frame,time,resolution,name,add_contour=None,vmin=None,vmax=None,vcenter=0,contour=[None],color=None):
    '''
    The _`plot_load_frame` function plot the data in entry into a defined frame following different arguments fro the plot.
    
    Attribute :
    ----------- 
        ax : matplotlib.pyplot.axes
            The axis on wich the data will be plot.
        data : LOAD_TIME_GRID
            A LOAD_TIME_GRID object with the LOAD_TIME_GRID.coeff to be ploted.
        frame : (lon1,lon2,lat1,lat2)
            The coordinate of the frame to plot.
        time : float
            The time of the plot that will be plot
        resolution : int
            The resolution of the plot, the resulting resolution will be resolution x resolution*2
        name : str
            The name of the plot, will be used in the title and legend of the plot
        vmin : float
            The minimal value of the plot colormap
        vmax : float
            The maximal value of the plot colormap
        vcenter : float
            The central value of the plot colormap
        contour : np.array([n_c])
            The contour that will be ploted on the map
        color : np.array([n_c])
            The contour fill color to be ploted
        

    Returns : 
    --------- 
        ax : matplotlib.pyplot.axes
            The axis updated from the plots
        
    '''
    t_it=np.where(data.time_step==time)[0][0]
    data.coeff=(data.coeff_from_step(t_it-2).coeff-data.coeff_from_step(t_it-3).coeff)/(data.time_step[t_it-1]-data.time_step[t_it])
    grid,lon_hd,lat_hd=data.coefftogrdhd(resolution)
    lon_hd-=180
    grid=np.concatenate((grid[:,resolution-1:],grid[:,:resolution]),axis=1)
    lon_lim_min=np.abs(lon_hd-frame[0]).argmin()
    lon_lim_max=np.abs(lon_hd-frame[1]).argmin()
    lat_lim_min=np.abs(lat_hd-frame[2]).argmin()
    lat_lim_max=np.abs(lat_hd-frame[3]).argmin()
    lon=lon_hd[lon_lim_min:lon_lim_max+1]
    lat=lat_hd[lat_lim_min:lat_lim_max+1]
    grid=grid[lat_lim_min:lat_lim_max+1,lon_lim_min:lon_lim_max+1]


    alpha_ocean=0
    coast_line_width=0.5
    cmap=cm.get_cmap('bwr', 100)
    #cmap=mpl.colors.LinearSegmentedColormap.from_list('hsv_2',mpl.colormaps['hsv_r'].resampled(255)(range(255))[20:-20],gamma=1.0)
    
    if not(color is None):
        cmap=None

    norm = colors.TwoSlopeNorm(vmin=vmin,vmax=vmax,vcenter=vcenter)
    if not(contour[0] is None):
        if contour[0]<0 and contour[-1]>0 :
            norm = colors.TwoSlopeNorm(vmin=contour[0],vmax=contour[-1],vcenter=0)
        elif contour[0]<0 and contour[-1]<0 :
            norm  = colors.TwoSlopeNorm(vmin=contour[0],vmax=contour[-1],vcenter=contour[-1])
        elif contour[0]>0 and contour[-1]>0 :
            norm  = colors.TwoSlopeNorm(vmin=contour[0],vmax=contour[-1],vcenter=contour[0])
    
    #print(frame_corr)
    ax.set_extent(frame)
    #ax.set_global()
    if contour[0] is None :
        m2 = ax.contourf(lon,lat,grid,levels=50,origin='lower', transform=ccrs.PlateCarree(central_longitude=0),extent=frame, zorder=0, cmap=cmap,colors=color,norm=norm)
    else :
        m2 = ax.contourf(lon,lat,grid,levels=contour,origin='lower', transform=ccrs.PlateCarree(central_longitude=0),extent=frame, zorder=0, cmap=cmap,colors=color,norm=norm)
    if not(add_contour is None):
        CS=ax.contour(lat,lon,add_contour,5,colors='k',linewidths=0.5,linestyles='solid')
        ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    
    cbar2=plt.colorbar(mappable=m2, orientation="horizontal")
    cbar2.set_label(f'{name} VGM at {time} kyr, averaged over {(data.time_step[t_it-1]-data.time_step[t_it])*1000} year (mm/yr); maximum subsidenc : {np.min(np.min(grid))}')
    # if not(contour is None):
    #     cbar2.set_ticks(contour)
    # elif not(vmax is None) and not(vmin is None):
    #     cbar2.set_ticks([vmin, 0, vmax])
    # elif not(vmin is None):
    #     cbar2.set_ticks([vmin, 0, grid.max()])
    # elif not(vmax is None):
    #     cbar2.set_ticks([grid.min(), 0, vmax])
    # else:
    #     cbar2.set_ticks([grid.min(), 0, grid.max()])
    cartopy.mpl.geoaxes.GeoAxes.gridlines(ax,crs=ccrs.PlateCarree(),draw_labels=True)
    ax.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
    ax.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)

    return ax

from netCDF4 import Dataset



def export_to_netcdf(Input_way,time,resolution,frame,save_way,type,type_bis):
    '''
    The _`export_to_netcdf` function export the data from the input_way, type and type_bis, of the frame shape into a netcdf file. This file can be open in any gis software.

    Attribute : 
    -----------
        Input_way : str
            way where the load data are located. If you are using the function from SL_C0de.SOLVER library, this way should be xxx/model_output/earth_model_name.
        time : float
            The time of the plot in kyr
        resolution : int
            the spatial resolution of the plot. ! This resolution concern the entire globe and not juste the frame.
        frame : (lon1,lon2,lat1,lat2)
            The coordinate of the frame to plot.
        save_way : str
            The filepath to save the netcdf file
        type : str 
            The type of loading taken into account. Can be : "SEDIMENT", "OCEAN", "ICE", "NO_SEDIMENT".
        type_bis : str
            The type of affected surface by the loading. Can be : "GEOID", "LOAD".
            
    ''' 
    data=LOAD_TIME_GRID(from_file=(True,f'{Input_way}/LOAD/{type}_{type_bis}'))
    t_it=np.where(data.time_step==time)[0][0]
    data.height_time_coeff=data.viscuous_deformation+data.elastic_deformation
    data.coeff=(data.coeff_from_step(t_it-2).coeff-data.coeff_from_step(t_it-3).coeff)/(data.time_step[t_it-1]-data.time_step[t_it])
    grid,lon_hd,lat_hd=data.coefftogrdhd(resolution)
    lon_hd-=180
    grid=np.concatenate((grid[:,resolution-1:],grid[:,:resolution]),axis=1)
    lon_lim_min=np.abs(lon_hd-frame[0]).argmin()
    lon_lim_max=np.abs(lon_hd-frame[1]).argmin()
    lat_lim_min=np.abs(lat_hd-frame[2]).argmin()
    lat_lim_max=np.abs(lat_hd-frame[3]).argmin()
    lon=lon_hd[lon_lim_min:lon_lim_max+1]
    lat=lat_hd[lat_lim_min:lat_lim_max+1]
    grid=grid[lat_lim_min:lat_lim_max+1,lon_lim_min:lon_lim_max+1]

    data.ncgrid=Dataset(f'{save_way}/{data.time_grid_name}{frame}.nc','w','NETCDF4_CLASSIC')
    data.ncgrid.title=data.time_grid_name

    data.ncgrid.createDimension('lon',len(lon))
    data.ncgrid.createDimension('lat',len(lat))


    lati=data.ncgrid.createVariable('lat', np.float32, ('lat',))
    lati.units = 'degrees_north'
    lati.long_name = 'latitude'
    lati[:]=lat

    long=data.ncgrid.createVariable('lon', np.float32, ('lon',))
    long.units = 'degrees_east'
    long.long_name = 'longitude'
    long[:]=lon

    grided=data.ncgrid.createVariable('grid', np.float32, ('lat','lon'))
    grided.units = 'mm/yr'
    grided.long_name = 'grid'
    grided[:,:]=grid

    data.ncgrid.close()
