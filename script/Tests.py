import numpy as np

from grid import TIME_GRID
import matplotlib.pyplot as plt

time_grid=TIME_GRID(6,64)
nlon=8
nlat=4
data_grid=np.random.rand(4,nlat,nlon).round()

lon=np.linspace(1,359,nlon)
lat=np.linspace(-89,89,nlat)



time_grid.interp_on_time_and_space(data_grid,np.array([3,2.6,2,1,0]),np.array([3,2.5,2,1.5,1,0.5,0]),6,lon,lat)

time_grid.plot_step_on_sphere(0)
time_grid.scatter_step_on_sphere(0)
time_grid.plot_step_on_Mollweide_projection(0,alpha_ocean=0)
plt.show()