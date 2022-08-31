import numpy as np
from scipy import interpolate
def interp_on(grd,lon,lat,inlon,inlat):
    #if grd.shape!=lon.shape :
    #    lon,lat=np.meshgrid(lon,lat)
    elons,lats=np.meshgrid(inlon,inlat)
    f = interpolate.RegularGridInterpolator((lon[1,:],lat[:,1][-1::-1]),grd.transpose()[:,-1::-1],method='linear')
    grd=f((elons.flatten(),lats.flatten()))
    grd=np.reshape(grd,(inlat.shape[0],inlon.shape[0]))
    return grd

def f_sed_quick_load(sed_in,sed_long,sed_lat,inlat,inlon) :
    sed_interp = interp_on(sed_in.squeeze(),sed_long,sed_lat,inlon,inlat)
    return sed_interp