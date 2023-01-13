import numpy as np
from scipy import interpolate
"""
def f_ice_quick_load(ice_time,ice_extended,ice_extended_2,ice_in,t_it,ice_long,ice_lat,inlat,inlon) :
    ice_nointerp=ice_in[t_it]
    ice_extended=np.concatenate((np.zeros((1,len(ice_long)-2)),ice_nointerp,ice_nointerp[0,-1]*np.ones((1,len(ice_long)-2))),axis=0)
    ice_extended_2=np.concatenate((np.expand_dims(ice_extended[:,-1],axis=1),ice_extended,np.expand_dims(ice_extended[:,0],axis=1)),axis=1)

    ice_interp = grid.interp_on(ice_extended_2,ice_long,ice_lat)
    
    time.sleep(0.01)
    self.ice[len(grid.time_step)-t_it-1,:,:] = ice_interp
    grid.time_step[t_it] = ice_time[len(ice_time)-t_it-1]
"""
    
def f_ice_quick_load(ice_time,ice_in,t_it,ice_long,ice_lat,inlat,inlon) :
    ice_nointerp=ice_in.squeeze()
    ice_extended=np.concatenate((np.zeros((1,len(ice_long)-2)),ice_nointerp,ice_nointerp[0,-1]*np.ones((1,len(ice_long)-2))),axis=0)
    ice_extended_2=np.concatenate((np.expand_dims(ice_extended[:,-1],axis=1),ice_extended,np.expand_dims(ice_extended[:,0],axis=1)),axis=1)
    ice_interp = interp_on(ice_extended_2,ice_long,ice_lat,inlon,inlat)
    return ice_interp, ice_time[len(ice_time)-t_it-1]

def interp_on(grd,lon,lat,inlon,inlat):
    #if grd.shape!=lon.shape :
    #    lon,lat=np.meshgrid(lon,lat)
    elons,lats=np.meshgrid(inlon,inlat)
    f = interpolate.RegularGridInterpolator((lon.transpose().squeeze(),lat.transpose().squeeze()[-1::-1]),grd.transpose()[:,-1::-1],method='linear')
    grd=f((elons.flatten(),lats.flatten()))
    grd=np.reshape(grd,(inlat.shape[0],inlon.shape[0]))
    return grd