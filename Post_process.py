from SL_C0de.grid import OCEAN_TIME_GRID
from SL_C0de.grid import SEDIMENT_TIME_GRID
from SL_C0de.grid import ICE_TIME_GRID
from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID

Input_way='C:/Users/ahenry01/Desktop/Python_code/Output_grid/VM5a_122'

Output_way=Input_way+'/'

ocean_time_grid=OCEAN_TIME_GRID(from_file=(True,Input_way+'/OCEAN_122_512'))
ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+'/ice_ICE6G_122_512'))

Input_way='C:/Users/ahenry01/Desktop/Python_code/Interpolated_grid'
sed_time_grid=SEDIMENT_TIME_GRID(from_file=(True,Input_way+'/SED/sed_AYS1_122_512'))

from SL_C0de.love import LOVE
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
a=6371000
Me=5.9742e24
love_number=LOVE(ice_time_grid.maxdeg,'C:/Users/ahenry01/Desktop/Python_code/SL_C0de_data/VM5a_122',ice_time_grid.time_step,a,Me)

from SL_C0de.grid import LOAD_TIME_GRID
sed_time_grid.timegrdtotimecoeff()
ice_time_grid.timegrdtotimecoeff()
sdelL_total=ice_time_grid.height_time_coeff*ice_time_grid.rho+sed_time_grid.height_time_coeff*sed_time_grid.rho+ocean_time_grid.height_time_coeff*ocean_time_grid.rho

#Calculating the load deformation

beta_l=love_number.beta_R_l

ice_load_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='ICE_LOAD_122_512')
ice_load_time_grid.calc_elastic_time()
ice_load_time_grid.calc_viscuous_time(backend=True)
ice_load_time_grid.save(save_way=Output_way)
ice_load_time_grid.clean_memory()

sediment_load_time_grid=LOAD_TIME_GRID(sdelL=sed_time_grid.height_time_coeff*sed_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='SEDIMENT_LOAD_122_512')
sediment_load_time_grid.calc_elastic_time()
sediment_load_time_grid.calc_viscuous_time(backend=True)
sediment_load_time_grid.save(save_way=Output_way)
sediment_load_time_grid.clean_memory()

ocean_load_time_grid=LOAD_TIME_GRID(sdelL=ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEAN_LOAD_122_512')
ocean_load_time_grid.calc_elastic_time()
ocean_load_time_grid.calc_viscuous_time(backend=True)
ocean_load_time_grid.save(save_way=Output_way)
ocean_load_time_grid.clean_memory()

total_load_time_grid=LOAD_TIME_GRID(sdelL=sdelL_total,beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='TOTAL_LOAD_122_512')
total_load_time_grid.calc_elastic_time()
total_load_time_grid.calc_viscuous_time(backend=True)
total_load_time_grid.save(save_way=Output_way)
total_load_time_grid.clean_memory()

#Calculating the geoïd deformation

beta_l=love_number.beta_G_l

ice_geoid_time_grid=LOAD_TIME_GRID(sdelL=ice_time_grid.height_time_coeff*ice_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='ICE_GEOID_122_512')
ice_geoid_time_grid.calc_elastic_time()
ice_geoid_time_grid.calc_viscuous_time(backend=True)
ice_geoid_time_grid.save(save_way=Output_way)
ice_geoid_time_grid.clean_memory()

sediment_geoid_time_grid=LOAD_TIME_GRID(sdelL=sed_time_grid.height_time_coeff*sed_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='SEDIMENT_GEOID_122_512')
sediment_geoid_time_grid.calc_elastic_time()
sediment_geoid_time_grid.calc_viscuous_time(backend=True)
sediment_geoid_time_grid.save(save_way=Output_way)
sediment_geoid_time_grid.clean_memory()

ocean_geoid_time_grid=LOAD_TIME_GRID(sdelL=ocean_time_grid.height_time_coeff*ocean_time_grid.rho,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEAN_GEOID_122_512')
ocean_geoid_time_grid.calc_elastic_time()
ocean_geoid_time_grid.calc_viscuous_time(backend=True)
ocean_geoid_time_grid.save(save_way=Output_way)
ocean_geoid_time_grid.clean_memory()

total_geoid_time_grid=LOAD_TIME_GRID(sdelL=sdelL_total,beta_l=beta_l,E=love_number.k_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='TOTAL_GEOID_122_512')
total_geoid_time_grid.calc_elastic_time()
total_geoid_time_grid.calc_viscuous_time(backend=True)
total_geoid_time_grid.save(save_way=Output_way)
total_geoid_time_grid.clean_memory()

# Calculate the subsidence induced by the ocean replaced by sediment :
# I have to crop the sediment grid by the ocean grid. 

from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID

beta_l=love_number.beta_R_l # Load the earth love numbers
#Preparing a new grid that compute the load of the substracted sediment volume to the ocean.
oceanic_sediment_load_time_grid=LOAD_TIME_GRID(beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEANIC_SEDIMENT_LOAD_122_512')

Input_way='C:/Users/ahenry01/Desktop/Python_code/Output_grid/VM5a_122'
#Loading the topography to compute the ocean function at each time step
topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/topo_topo_SL_122_512'))

#We create a new object to compute the oceanic sediment
oceanic_sediment_time_grid=sed_time_grid
oceanic_sediment_time_grid.rho=ocean_time_grid.rho

for t_it in range (ocean_load_time_grid.time_step_number-1): # At each time step we apply the ocean function to the sediment height grid
    oceanic_sediment_time_grid.height_time_grid[t_it,:,:]=oceanic_sediment_time_grid.height_time_grid[t_it,:,:]*ocean_time_grid.evaluate_ocean(topo_time_grid.height_time_grid[t_it,:,:]).grd
oceanic_sediment_time_grid.timegrdtotimecoeff()

#computing the earth deformation
oceanic_sediment_load_time_grid.load=oceanic_sediment_time_grid.height_time_coeff*oceanic_sediment_time_grid.rho
oceanic_sediment_load_time_grid.calc_elastic_time()
oceanic_sediment_load_time_grid.calc_viscuous_time(backend=True)
oceanic_sediment_load_time_grid.save(save_way=Output_way)
oceanic_sediment_load_time_grid.clean_memory()

beta_l=love_number.beta_G_l # Load the ocean love numbers

#computing the geoid deformation
oceanic_sediment_geoid_time_grid=LOAD_TIME_GRID(beta_l=beta_l,E=love_number.h_e,a=a,Me=Me,time_step=ice_time_grid.time_step,maxdeg=ice_time_grid.maxdeg,grid_name='OCEANIC_SEDIMENT_GEOID_122_512')
oceanic_sediment_geoid_time_grid.load=oceanic_sediment_time_grid.height_time_coeff*oceanic_sediment_time_grid.rho
oceanic_sediment_geoid_time_grid.calc_elastic_time()
oceanic_sediment_geoid_time_grid.calc_viscuous_time(backend=True)
oceanic_sediment_geoid_time_grid.save(save_way=Output_way)
oceanic_sediment_geoid_time_grid.clean_memory()