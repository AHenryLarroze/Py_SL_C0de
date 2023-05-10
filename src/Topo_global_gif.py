import imageio
import cartopy
import cartopy.crs as ccrs
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID

import sys


maxdeg=512


Input_way='C:/Users/ahenry01/Desktop/Python_code/Output_grid/VM5a_122'

topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/topo_topo_SL_122_512'))

grd=topo_time_grid

alpha_ocean=0
coast_line_width=0.5

cmap=cm.get_cmap('bwr', 100)


for t_it in range (grd.time_step_number):

    #grid=grd.ice[:t_it+1,:,:].sum(0)

    
    grid=grd.height_time_grid[t_it,:,:]


    fig = plt.figure(figsize=(12, 12), facecolor="none")
    ax  = plt.subplot(111, projection=ccrs.Mollweide())
    ax.set_global()
    norm = colors.TwoSlopeNorm(vmin=None,vmax=None,vcenter=0)
    m = ax.imshow(grid, origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, cmap=cmap, interpolation="gaussian",norm=norm)
    ax.contour(topo_time_grid.elons,-topo_time_grid.lats,grid,[0], transform=ccrs.PlateCarree())
    cbar=plt.colorbar(mappable=m, orientation="horizontal", shrink=0.8)
    cbar.set_label(f'Topography (m) at {grd.time_step[t_it]}')
    ax.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
    ax.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)

    fig.savefig(f'temp_gif/at_{grd.time_step[t_it]}.png',transparent = False, facecolor= 'white')

    plt.close()

frames=[]
for t_it in range (grd.time_step_number):
    image = imageio.v2.imread(f'temp_gif/at_{grd.time_step[t_it]}.png')
    frames.append(image)

imageio.mimsave('Topography_122.gif',frames,duration = 100) 