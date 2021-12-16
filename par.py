import numpy as np
from scipy import interpolate

def f_ice_quick_load(ice_time,ice_extended,ice_extended_2,ice_in,t_it,ice_long,ice_lat,inlat,inlon) :
    ice_extended[0,:]= np.zeros((1,len(ice_long)-3))
    ice_extended[1:-2,:]=ice_in
    ice_extended[-1,:]=ice_in[0,-1]*np.ones((1,len(ice_long)-3))
    
    ice_extended_2[:,0]=ice_extended[:,-1]
    ice_extended_2[:,1:-2]=ice_extended
    ice_extended_2[:,-1] = ice_extended[:,0]

    ice_interp = interp_on(ice_extended_2,ice_long,ice_lat,inlon,inlat)
    return ice_interp, ice_time[len(ice_time)-t_it-1]
    
def f_print(ice_time,t_it,ice_in):
    return ice_time[t_it],ice_in

def interp_on(grd,lon,lat,inlon,inlat):
    if grd.shape!=lon.shape :
        lon,lat=np.meshgrid(lon,lat)
    elons,lats=np.meshgrid(inlon,inlat)
    grd = interpolate.griddata((lon.flatten(),lat.flatten()),grd.flatten(),(elons.flatten(),lats.flatten()),method='nearest')
    grd=np.reshape(grd,(inlat.shape[0],inlon.shape[0]))
    return grd