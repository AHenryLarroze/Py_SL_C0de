from scipy import io
import numpy as np
from slcode import grid
import matplotlib.pyplot as plt


data=io.loadmat('C:/Users/ahenry01/Desktop/Python_code/SL_C0de_data/ice6G122k') #load the file, can be modified to load an other file.
ice_in=data['ice6g']
ice_time=data['ice_time']
ice_time=ice_time.squeeze()
ice_lon =data['ice_long'].squeeze()
ice_lat =data['ice_lat'].squeeze()[::-1]


# lon,lat=np.meshgrid(ice_lon,ice_lat)
# plt.pcolor(ice_in.T[0,:,:])
# plt.show()

time_step=np.linspace(120,0,120,-1)
maxdeg=64

print(len(ice_lon)*len(ice_lat))
print(ice_in.T.shape)

ice_time_grid=TIME_GRID(len(time_step),maxdeg)


ice_time_grid.interp_on_time_and_space(ice_in.T,ice_time,time_step,len(time_step),ice_lon,ice_lat,backend=True)
ice_time_grid.scatter_step_on_sphere(20)
# ice_time_grid.plot_step_on_Mollweide_projection(119,alpha_ocean=0)
plt.show()