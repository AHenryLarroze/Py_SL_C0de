import matplotlib.pyplot as plt
import math 
from SL_C0de.grid import OCEAN_TIME_GRID
from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID
from SL_C0de.grid import ICE_TIME_GRID

Input_way='C:/Users/ahenry01/Desktop/Python_code/Output_grid/VM5a_122'

ocean_time_grid_26=OCEAN_TIME_GRID(from_file=(True,Input_way+'/OCEAN_26_512'))
topo_time_grid_26=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/topo_topo_SL_26_512'))
ice_time_grid_26=ICE_TIME_GRID(from_file=(True,Input_way+'/ice_ICE6G_26_512'))
topo_time_grid_26.topo_pres=topo_time_grid_26.height_time_grid[-1,:,:].copy()

ocean_time_grid_26.timegrdtotimecoeff()
ice_time_grid_26.timegrdtotimecoeff()
ESL_26=ice_time_grid_26.height_time_coeff[:,0]/ocean_time_grid_26.evaluate_ocean(topo_time_grid_26.topo_pres).grdtocoeff().coeff[0]*(ice_time_grid_26.rho/ocean_time_grid_26.rho)

ocean_time_grid_122=OCEAN_TIME_GRID(from_file=(True,Input_way+'/OCEAN_122_512'))
topo_time_grid_122=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/topo_topo_SL_122_512'))
ice_time_grid_122=ICE_TIME_GRID(from_file=(True,Input_way+'/ice_ICE6G_122_512'))
topo_time_grid_122.topo_pres=topo_time_grid_122.height_time_grid[-1,:,:].copy()

ocean_time_grid_122.timegrdtotimecoeff()
ice_time_grid_122.timegrdtotimecoeff()
ESL_122=ice_time_grid_122.height_time_coeff[:,0]/ocean_time_grid_122.evaluate_ocean(topo_time_grid_122.topo_pres).grdtocoeff().coeff[0]*(ice_time_grid_122.rho/ocean_time_grid_122.rho)


fig=plt.figure(facecolor="none")
plt.plot(ocean_time_grid_26.time_step,(ocean_time_grid_26.height_time_coeff[:,0].cumsum()-ocean_time_grid_26.height_time_coeff[:,0].sum(0)),label='ESL for a variable ocean over 26 ka')
plt.plot(ocean_time_grid_26.time_step,ESL_26[::-1].cumsum()[::-1],label='ESL over 26 ka')
plt.plot(ocean_time_grid_122.time_step,(ocean_time_grid_122.height_time_coeff[:,0].cumsum()-ocean_time_grid_122.height_time_coeff[:,0].sum(0)),label='ESL for a variable ocean over 122 ka')
plt.plot(ocean_time_grid_122.time_step,ESL_122[::-1].cumsum()[::-1],label='ESL over 122 ka')
plt.legend()

fig.savefig('Sea_level_curve_comparison.pdf')

plt.close()

del_v_ocean=(ocean_time_grid_26.height_time_coeff[:,0].cumsum()-ocean_time_grid_26.height_time_coeff[:,0].sum(0))[::2]-(ocean_time_grid_122.height_time_coeff[:,0].cumsum()-ocean_time_grid_122.height_time_coeff[:,0].sum(0))[-26:]

fig=plt.figure(facecolor="none")
plt.plot(ocean_time_grid_26.time_step[::2],del_v_ocean,label='difference ESL for a variable ocean over 26 ka and 122 ka')
plt.plot(ocean_time_grid_26.time_step[::2],ESL_26[::-1].cumsum()[::-1][::2]-ESL_122[::-1].cumsum()[::-1][-26:],label='difference ESL over 26 ka and 122 ka')
plt.legend()

fig.savefig('Sea_level_curve_difference.pdf')

plt.close()