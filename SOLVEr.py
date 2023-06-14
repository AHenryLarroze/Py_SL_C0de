from SL_C0de.grid import ICE_TIME_GRID
from SL_C0de.grid import SEDIMENT_TIME_GRID
from SL_C0de.grid import TOPOGRAPHIC_TIME_GRID
from SL_C0de.grid import OCEAN_TIME_GRID
from SL_C0de.love import LOVE
import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)

maxdeg=512
stop=26 # define the number of time steps
step=0.5
time_step=np.arange(start=stop,stop=-step,step=-step)

Input_way='C:/Users/ahenry01/Desktop/Python_code/Interpolated_grid'
love_way='C:/Users/ahenry01/Desktop/Python_code/SL_C0de_data'
love_cat='/variable_lm'
love_file='/VM5a.l60.um21.lm22.949'


ice_time_grid=ICE_TIME_GRID(from_file=(True,Input_way+'/ICE/ice_ICE6G_26_512'))
sed_time_grid=SEDIMENT_TIME_GRID(from_file=(True,Input_way+'/SED/sed_AYS1_26_512'))
topo_time_grid=TOPOGRAPHIC_TIME_GRID(from_file=(True,Input_way+'/TOPO/topo_topo_SL_AYS1_26_512'))

Output_way='F:/SL_C0de/AYEYARWADY/26kyr_AYS1'

# Initiate the base elements
from SL_C0de.spharm import sphericalobject
from SL_C0de.Load import LOAD

ocean_time_grid=OCEAN_TIME_GRID(time_step,maxdeg,grid_name='OCEAN_AYS1_26_512')
ocean_time_grid.time_step_number=ocean_time_grid.time_step_number

ice_time_grid.timegrdtotimecoeff()
sed_time_grid.timegrdtotimecoeff()

love_number=LOVE(maxdeg,love_way+love_cat+love_file,time_step,6371000,5.9742e24)
love_number.calc_beta()
TO=sphericalobject(coeff=np.zeros(ice_time_grid.height_time_coeff[0,:].shape))

load=LOAD(maxdeg,ice_time_grid.time_step)

topo_diff_median=np.inf
#sdel_topo_diff=np.inf
topo_it=0
while topo_diff_median>10**(0) : #and sdel_topo_diff>10**(-1):
    delPhi_g_time=np.array([])
    TO.prev=np.zeros(ice_time_grid.height_time_coeff[0,:].shape)

    if topo_diff_median != np.inf :
        topo_time_grid.height_time_grid[0,:,:]=topo_initial.copy()
    # topo_time_grid.height_time_grid = topo_time_grid.height_time_grid - ice_time_grid.height_time_grid + ice_time_grid.ice # for resetting the corrected ice.
    ice_time_grid.ice_correction(topo_time_grid,ocean_time_grid)
    ice_time_grid.timegrdtotimecoeff()
    
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
        


        conv_lim=10**(-10)
        print('time_iteration : ',ice_time_grid.time_step[t_it])
        if topo_it==0 : 
            conv_it=0
        else :
            conv_it=1
        conv_it=ocean_time_grid.sea_level_solver(load,ice_time_grid,sed_time_grid,love_number,TO,t_it,conv_it,conv_lim)

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

    #Saving the code :
ocean_time_grid.timecoefftotimegrd()
import os 
if not(os.path.exists(Output_way+love_file)):
    os.mkdir(Output_way+love_file)
ocean_time_grid.save(Output_way+love_file)
#ice_time_grid.time_grid_name="ice_ICE6G_122_512"
ice_time_grid.save(Output_way+love_file)
#topo_time_grid.time_grid_name="topo_topo_SL_122_512"
topo_time_grid.save(Output_way+love_file)