import imageio
import cartopy
import cartopy.crs as ccrs
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID

import sys

np.set_printoptions(threshold=sys.maxsize)

maxdeg=512
stop=26 # define the number of time steps
step=0.5
time_step=np.arange(start=stop,stop=-step,step=-step)


Input_way='F:/SL_C0de/AYEYARWADY/Test/VM5a.l32.um21.lm22.699'

topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/topo_topo_SL_AYS2_26_512'))


grd=topo_time_grid

alpha_ocean=0
coast_line_width=0.5
cmap=cm.get_cmap('bwr', 100)


area=(93,98,13,20)
#find the closest position :
lon_lim_min=np.abs(grd.elons-area[0]).argmin()-1
lon_lim_max=np.abs(grd.elons-area[1]).argmin()+1
lat_lim_min=np.abs(grd.lats[::-1]-area[2]).argmin()-1
lat_lim_max=np.abs(grd.lats[::-1]-area[3]).argmin()+2
lon=grd.elons[lon_lim_min:lon_lim_max]
lat=grd.lats[::-1][lat_lim_min:lat_lim_max]

#define colors :

color_topo=[(68,79,137),(68,79,137),(68,101,137),(68,101,137),(38,117,207),(38,117,207),(102,153,205),(102,153,205),(122,182,245),(102,205,171),(158,215,194),(137,205,102),(205,205,102),(114,137,68),(114,137,68),(205,170,102),(205,137,102),(205,102,102),(205,102,102),(137,68,68),(137,68,68)]
color_topo=[(c[0]/255,c[1]/255,c[2]/255) for c in color_topo]
levels_topo=[-8000,-7000,-6000,-5000,-4000,-3000,-2000,-1000,-200,-100,-10,0,10,50,100,200,1000,2000,3000,4000,5000]

for t_it in range (grd.time_step_number-1):


    
    grid=grd.height_time_grid[t_it,:,:]#+sed_time_grid.height_time_grid[:t_it+1,:,:].sum(0)
    grid=grid[lat_lim_min:lat_lim_max,lon_lim_min:lon_lim_max]

    

    fig = plt.figure(figsize=(12, 12), facecolor="none")
    ax  = plt.subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent(area)
    m = ax.contourf(lon,lat,grid,levels=levels_topo, origin='lower', transform=ccrs.PlateCarree(),extent=[0,360, -89, 89], zorder=0, colors=color_topo)
    ax.contour(lon,lat,grid,[0], transform=ccrs.PlateCarree(),zorder=1)
    cbar=plt.colorbar(mappable=m, orientation="horizontal", shrink=0.8)
    cbar.set_label(f'Topography (m) at {grd.time_step[t_it]}')
    ax.add_feature(cartopy.feature.OCEAN, alpha=alpha_ocean, zorder=99, facecolor="#BBBBBB")
    ax.coastlines(resolution="50m", zorder=100, linewidth=coast_line_width)

    fig.savefig(f'temp_gif/at_{grd.time_step[t_it]}.png',transparent = False, facecolor= 'white')

    plt.close()

frames=[]
for t_it in range (grd.time_step_number-1):
    image = imageio.v2.imread(f'temp_gif/at_{grd.time_step[t_it]}.png')
    frames.append(image)

imageio.mimsave(Input_way+'Topography_local_122.gif',frames,duration = 100) 

