import matplotlib.pyplot as plt
import math 
from SL_C0de.grid import OCEAN_TIME_GRID
from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID
from SL_C0de.grid import ICE_TIME_GRID

Input_way='C:/Users/ahenry01/Desktop/Python_code/Output_grid/VM5a_122'

ocean_time_grid=OCEAN_TIME_GRID(from_file=(True,Input_way+'/OCEAN_26_512'))
topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/topo_topo_SL_26_512'))
ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+'/ice_ICE6G_26_512'))
topo_time_grid.topo_pres=topo_time_grid.height_time_grid[-1,:,:].copy()

ocean_time_grid.timegrdtotimecoeff()
ice_time_grid.timegrdtotimecoeff()
ESL=ice_time_grid.height_time_coeff[:,0]/ocean_time_grid.evaluate_ocean(topo_time_grid.topo_pres).grdtocoeff().coeff[0]*(ice_time_grid.rho/ocean_time_grid.rho)
fig=plt.figure(facecolor="none")
plt.plot(ocean_time_grid.time_step,(ocean_time_grid.height_time_coeff[:,0].cumsum()-ocean_time_grid.height_time_coeff[:,0].sum(0)),label='ESL for a variable ocean')
plt.plot(ocean_time_grid.time_step,ESL[::-1].cumsum()[::-1],label='ESL')
plt.legend()

fig.savefig('Sea_level_curve.pdf')

plt.close()

