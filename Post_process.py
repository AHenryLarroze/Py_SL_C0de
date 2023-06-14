from SL_C0de.grid import OCEAN_TIME_GRID
from SL_C0de.grid import SEDIMENT_TIME_GRID
from SL_C0de.grid import ICE_TIME_GRID
from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID

import os

def find_files(filename, search_path):
    result = []

    # Wlaking top-down from the root
    for root, dir, files in os.walk(search_path):
        if filename in dir:
            result.append(os.path.join(root, filename))
    return result

def calculate_deformation(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,Output_way,backend=False) :

    beta_l=love_number.beta_R_l

    ice_load_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='ICE_LOAD_122_512')
    ice_load_time_grid.calc_elastic_time()
    ice_load_time_grid.calc_viscuous_time(backend=backend)
    ice_load_time_grid.save(save_way=Output_way)
    ice_load_time_grid=0

    sediment_load_time_grid=LOAD_TIME_GRID(sdelL=sed_time_grid.height_time_coeff*sed_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='SEDIMENT_LOAD_122_512')
    sediment_load_time_grid.calc_elastic_time()
    sediment_load_time_grid.calc_viscuous_time(backend=backend)
    sediment_load_time_grid.save(save_way=Output_way)
    sediment_load_time_grid=0

    ocean_load_time_grid=LOAD_TIME_GRID(sdelL=ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEAN_LOAD_122_512')
    ocean_load_time_grid.calc_elastic_time()
    ocean_load_time_grid.calc_viscuous_time(backend=backend)
    ocean_load_time_grid.save(save_way=Output_way)
    ocean_load_time_grid=0

    total_load_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho+sed_time_grid.height_time_coeff*sed_time_grid.rho+ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='TOTAL_LOAD_122_512')
    total_load_time_grid.calc_elastic_time()
    total_load_time_grid.calc_viscuous_time(backend=backend)
    total_load_time_grid.save(save_way=Output_way)
    total_load_time_grid=0

    #Calculating the geoïd deformation

    beta_l=love_number.beta_G_l

    ice_geoid_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='ICE_GEOID_122_512')
    ice_geoid_time_grid.calc_elastic_time()
    ice_geoid_time_grid.calc_viscuous_time(backend=backend)
    ice_geoid_time_grid.save(save_way=Output_way)
    ice_geoid_time_grid=0

    sediment_geoid_time_grid=LOAD_TIME_GRID(sdelL=sed_time_grid.height_time_coeff*sed_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='SEDIMENT_GEOID_122_512')
    sediment_geoid_time_grid.calc_elastic_time()
    sediment_geoid_time_grid.calc_viscuous_time(backend=backend)
    sediment_geoid_time_grid.save(save_way=Output_way)
    sediment_geoid_time_grid=0

    ocean_geoid_time_grid=LOAD_TIME_GRID(sdelL=ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEAN_GEOID_122_512')
    ocean_geoid_time_grid.calc_elastic_time()
    ocean_geoid_time_grid.calc_viscuous_time(backend=backend)
    ocean_geoid_time_grid.save(save_way=Output_way)
    ocean_geoid_time_grid=0

    total_geoid_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho+sed_time_grid.height_time_coeff*sed_time_grid.rho+ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='TOTAL_GEOID_122_512')
    total_geoid_time_grid.calc_elastic_time()
    total_geoid_time_grid.calc_viscuous_time(backend=backend)
    total_geoid_time_grid.save(save_way=Output_way)
    total_geoid_time_grid=0

import copy

def calculate_sediment_ocean_interaction(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,topo_time_grid,Output_way,backend=False) :
    beta_l=love_number.beta_R_l # Load the earth love numbers
    #Preparing a new grid that compute the load of the substracted sediment volume to the ocean.
    #We create a new object to compute the oceanic sediment
    oceanic_sediment_time_grid=copy.copy(sed_time_grid)
    oceanic_sediment_time_grid.rho=ocean_time_grid.rho
    

    for t_it in range (oceanic_sediment_time_grid.time_step_number-1): # At each time step we apply the ocean function to the sediment height grid
        oceanic_sediment_time_grid.height_time_grid[t_it,:,:]=oceanic_sediment_time_grid.height_time_grid[t_it,:,:]*ocean_time_grid.evaluate_ocean(topo_time_grid.height_time_grid[t_it,:,:]).grd
    oceanic_sediment_time_grid.timegrdtotimecoeff()


    oceanic_sediment_load_time_grid=LOAD_TIME_GRID(sdelL=oceanic_sediment_time_grid.height_time_coeff*oceanic_sediment_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=oceanic_sediment_time_grid.time_step,maxdeg=oceanic_sediment_time_grid.maxdeg,grid_name='OCEANIC_SEDIMENT_LOAD_122_512')    
    #computing the earth deformation
    oceanic_sediment_load_time_grid.calc_elastic_time()
    oceanic_sediment_load_time_grid.calc_viscuous_time(backend=backend)
    oceanic_sediment_load_time_grid.save(save_way=Output_way)
    oceanic_sediment_load_time_grid.clean_memory()

    beta_l=love_number.beta_G_l # Load the ocean love numbers

    #computing the geoid deformation
    oceanic_sediment_geoid_time_grid=LOAD_TIME_GRID(sdelL=oceanic_sediment_time_grid.height_time_coeff*oceanic_sediment_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=oceanic_sediment_time_grid.time_step,maxdeg=oceanic_sediment_time_grid.maxdeg,grid_name='OCEANIC_SEDIMENT_GEOID_122_512')
    oceanic_sediment_geoid_time_grid.calc_elastic_time()
    oceanic_sediment_geoid_time_grid.calc_viscuous_time(backend=backend)
    oceanic_sediment_geoid_time_grid.save(save_way=Output_way)
    oceanic_sediment_geoid_time_grid.clean_memory()




Input_way='C:/Users/ahenry01/Desktop/Python_code/Interpolated_grid'
sed_time_grid=SEDIMENT_TIME_GRID(from_file=(True,Input_way+'/SED/sed_AYS2_26_512'))
sed_time_grid.timegrdtotimecoeff()


Input_way='F:/SL_C0de/AYEYARWADY/26kyr/'


earth_model_name_list=os.listdir(Input_way)

for earth_model_name in earth_model_name_list :

    print('calculation for : ' + earth_model_name)

    Output_way=Input_way+earth_model_name+'/LOAD/'

    if not(os.path.exists(Output_way)):
        os.mkdir(Output_way)

    ocean_time_grid=OCEAN_TIME_GRID(from_file=(True,Input_way+earth_model_name+'/OCEAN_AYS2_26_512'))
    ocean_time_grid.height_time_coeff=ocean_time_grid.height_time_coeff[:52,:]
    #print(ocean_time_grid.height_time_coeff.shape)
    ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+earth_model_name+'/ice_ICE6G_26_512'))
    #print(ice_time_grid.time_step_number)


    from SL_C0de.love import LOVE
    import numpy as np
    import sys
    np.set_printoptions(threshold=sys.maxsize)
    a=6371000
    Me=5.9742e24
    love_way=find_files(earth_model_name,'C:/Users/ahenry01/Desktop/Python_code/SL_C0de_data')[0]
    love_number=LOVE(ice_time_grid.maxdeg,love_way,ice_time_grid.time_step,a,Me)

    from SL_C0de.grid import LOAD_TIME_GRID
    ice_time_grid.timegrdtotimecoeff()

    calculate_deformation(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,Output_way)

    #Loading the topography to compute the ocean function at each time step
    topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+earth_model_name+'/topo_topo_SL_AYS2_26_512'))

    print('Oceanic sediment calculation')
    calculate_sediment_ocean_interaction(love_number,ice_time_grid,sed_time_grid,ocean_time_grid,a,Me,topo_time_grid,Output_way)

