import sys
import logging
from spharm import legendre, GaussQuad
import numpy as np
import math
from spharm import sphericalobject
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.offline import plot
from scipy import interpolate
from spharm import calc_at_point
from matplotlib import colors

def degree2radians(degree):
  # convert degrees to radians
  return degree*np.pi/180

def mapping_map_to_sphere(lon, lat, radius=1):
  # this function maps the points of coords (lon, lat) to points onto the sphere of radius radius
  lon=np.array(lon, dtype=np.float64)
  lat=np.array(lat, dtype=np.float64)
  lon=degree2radians(lon)
  lat=degree2radians(lat)
  xs=radius*np.cos(lon)*np.cos(lat)
  ys=radius*np.sin(lon)*np.cos(lat)
  zs=radius*np.sin(lat)
  return xs, ys, zs

def interp_on(grd,lon,lat,inlon,inlat):
    #if grd.shape!=lon.shape :
    #    lon,lat=np.meshgrid(lon,lat)
    elons,lats=np.meshgrid(inlon,inlat)
    f = interpolate.RegularGridInterpolator((lon.transpose().squeeze(),lat.transpose().squeeze()[-1::-1]),grd.transpose()[:,-1::-1],method='linear',bounds_error=False,fill_value=-1000)
    grd=f((elons.flatten(),lats.flatten()))
    grd=np.reshape(grd,(inlat.shape[0],inlon.shape[0]))
    return grd

def write_saved_to_csv(self_p,saved,way,name):
    lat,lon=np.meshgrid(self_p.lats_res,self_p.elons_res)

    time=0.5
    ind=np.where(self_p.grid.time_step==time)

    z=saved
    z=z[ind,:]
    z=sphericalobject(z.transpose().squeeze(),'coeff').coefftogrd_hd(self_p)

    z=z.transpose().reshape((z.size,1),order='A')
    x=lon.reshape((lon.size,1))
    y=lat.reshape((lat.size,1))

    data=np.concatenate((x-360,y,z),axis=1)
    np.savetxt(way+name+'.csv',data,header='x y z',delimiter=' ',comments='')

class world_plot(object):
    def __init__(self):
        ploting='on'

    def init_resolution_plot(self,N=64) :
        self.maxdeg_hd=N
        if N>self.maxdeg :
            self.P_lm_res=legendre(N,self.pool)
            self.x_res, self.w_res= GaussQuad(N)
            x_GL = np.arccos(self.x_res)*180/math.pi - 90
            lon_GL = np.linspace(0,360,2*N+1)
            lon_GL = lon_GL[:-1]
            self.lats_res=x_GL
            self.elons_res=lon_GL
        else :
            self.P_lm_res=self.P_lm
            self.x_res=self.x
            self.w_res=self.w
            self.lats_res=self.grid.lats
            self.elons_res=self.grid.elons

    def check_self_loaded(self):
        if not('grid' in self.__dict__) : 
            logging.error('error: ', "no self created in this object please use SL_self and SL_self.init_grided_parameters()")
            sys.exit(1)
        
    def plot_RSL(self,time_plot,vmin=-1000,vmax=1000) :
        self.check_self_loaded()
        if not('RSL' in self.SL.__dict__):
            self.calc_RSL()
        if not('P_lm_res' in self.__dict__):
            self.init_resolution_plot()
        ind=np.where(self.grid.time_step==time_plot)
        RSL_at_time=sphericalobject(self.SL.RSL[ind,:,:].squeeze()).grdtocoeff(self)
        RSL_grid=RSL_at_time.coefftogrd_hd(self)
        plt.pcolor(self.elons_res,self.lats_res,RSL_grid,cmap='seismic',vmin=vmin,vmax=vmax)
        plt.colorbar()
        plt.axis('equal')
        plt.contour(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),levels=[0])
        plt.show()
    
    def plot_map(self,c,time_plot,vmin=None,vmax=None,cmap='seismic',lat=[-90,90],lon=[0,360],ax=None,fig=None,title='',colorbar_title=''):
        self.check_self_loaded()
        ind=np.where(self.grid.time_step==time_plot)
        ind=ind[0][0]
        if not('local_topo' in self.__dict__):
            lat_min=np.absolute(self.grid.lats-lat[0]).argmin()
            lat_max=np.absolute(self.grid.lats-lat[1]).argmin()
            if lat_min>lat_max :
                temps=lat_max
                lat_max=lat_min
                lat_min=temps
            lon_min=np.absolute(self.grid.elons-lon[0]).argmin()
            lon_max=np.absolute(self.grid.elons-lon[1]).argmin()
            self.local_topo=self.topo.topo_pres[lat_min:lat_max,lon_min:lon_max]
            self.lati_topo=self.grid.lats[lat_min:lat_max]
            self.long_topo=self.grid.elons[lon_min:lon_max]
        if ax==None :
            fig,ax=plt.subplots()
        if len(c.shape)==2 :
            D_tot=sphericalobject(c[ind,:],'coeff').coefftogrd_hd(self)
            lat_min=np.absolute(self.lats_res-lat[0]).argmin()
            lat_max=np.absolute(self.lats_res-lat[1]).argmin()
            if lat_min>lat_max :
                temps=lat_max
                lat_max=lat_min
                lat_min=temps
            lon_min=np.absolute(self.elons_res-lon[0]).argmin()
            lon_max=np.absolute(self.elons_res-lon[1]).argmin()
            D=D_tot[lat_min:lat_max,lon_min:lon_max]
            lati=self.lats_res[lat_min:lat_max]
            long=self.elons_res[lon_min:lon_max]
            long,lati=np.meshgrid(long,lati)
            divnorm=colors.TwoSlopeNorm(vcenter=0.)
            if vmin==None :
                vmin=np.min(D)
                vmax=np.max(D)
                v=np.max(np.abs(np.array([vmin,vmax])))
                sp=ax.pcolor(long,lati,D,cmap=cmap,norm=divnorm)
            else :
                v=np.max(np.abs(np.array([vmin,vmax])))
                sp=ax.pcolor(long,lati,D,cmap=cmap,norm=divnorm)
            fig.colorbar(sp,label=colorbar_title)
            if 'topo' in self.topo.__dict__:
                ax.contour(self.long_topo,self.lati_topo,self.local_topo,levels=[0])
        elif len(c.shape)==3:
            D=c[ind,:,:].squeeze()
            if vmin==None :
                vmin=np.min(D)
                vmax=np.max(D)
                v=np.max(np.abs(np.array([vmin,vmax])))
                ax.pcolor(self.grid.elons,self.grid.lats,D,cmap=cmap,vmin=-v,vmax=v)
            else :
                ax.pcolor(self.grid.elons,self.grid.lats,D,cmap=cmap,vmin=vmin,vmax=vmax)
            fig.colorbar()
            plt.axis('equal')
            if 'topo' in self.topo.__dict__:
                plt.contour(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),levels=[0])
        plt.title(title)
        plt.grid()
        plt.axis('equal')
        return ax,D_tot,D,fig

    def plot_map_derivation(self,c,time_plot,vmin=None,vmax=None,cmap='seismic',lat=[-90,90],lon=[0,360],ax=None,fig=None,title='',colorbar_title=''):
        self.check_self_loaded()
        ind=np.where(self.grid.time_step==time_plot)
        ind=ind[0][0]
        if not('local_topo' in self.__dict__):
            lat_min=np.absolute(self.grid.lats-lat[0]).argmin()
            lat_max=np.absolute(self.grid.lats-lat[1]).argmin()
            if lat_min>lat_max :
                temps=lat_max
                lat_max=lat_min
                lat_min=temps
            lon_min=np.absolute(self.grid.elons-lon[0]).argmin()
            lon_max=np.absolute(self.grid.elons-lon[1]).argmin()
            self.local_topo=self.topo.topo_pres[lat_min:lat_max,lon_min:lon_max]
            self.lati_topo=self.grid.lats[lat_min:lat_max]
            self.long_topo=self.grid.elons[lon_min:lon_max]
        if ax==None :
            fig,ax=plt.subplots()
        if len(c.shape)==2 :
            D_tot=sphericalobject((c[ind,:]-c[ind-1,:])/(self.time_step[ind]-self.time_step[ind-1])/1000,'coeff').coefftogrd_hd(self)
            lat_min=np.absolute(self.lats_res-lat[0]).argmin()
            lat_max=np.absolute(self.lats_res-lat[1]).argmin()
            if lat_min>lat_max :
                temps=lat_max
                lat_max=lat_min
                lat_min=temps
            lon_min=np.absolute(self.elons_res-lon[0]).argmin()
            lon_max=np.absolute(self.elons_res-lon[1]).argmin()
            D=D_tot[lat_min:lat_max,lon_min:lon_max]
            lati=self.lats_res[lat_min:lat_max]
            long=self.elons_res[lon_min:lon_max]
            long,lati=np.meshgrid(long,lati)
            divnorm=colors.TwoSlopeNorm(vcenter=0.)
            if vmin==None :
                vmin=np.min(D)
                vmax=np.max(D)
                v=np.max(np.abs(np.array([vmin,vmax])))
                sp=ax.pcolor(long,lati,D,cmap=cmap,norm=divnorm)
            else :
                v=np.max(np.abs(np.array([vmin,vmax])))
                sp=ax.pcolor(long,lati,D,cmap=cmap,norm=divnorm)
            fig.colorbar(sp,label=colorbar_title)
            if 'topo' in self.topo.__dict__:
                ax.contour(self.long_topo,self.lati_topo,self.local_topo,levels=[0])
        elif len(c.shape)==3:
            D=c[ind,:,:]-c[ind-1,:,:]
            if vmin==None :
                vmin=np.min(D)
                vmax=np.max(D)
                v=np.max(np.abs(np.array([vmin,vmax])))
                sp=ax.pcolor(self.grid.elons,self.grid.lats,D,cmap=cmap,vmin=-v,vmax=v)
            else :
                sp=ax.pcolor(self.grid.elons,self.grid.lats,D,cmap=cmap,vmin=vmin,vmax=vmax)
            fig.colorbar(sp,label=colorbar_title)
            plt.axis('scaled')
            if 'topo' in self.topo.__dict__:
                plt.contour(self.grid.elons,self.grid.lats,self.topo.topo_pres.squeeze(),levels=[0])
        plt.title(title)
        plt.grid()
        plt.axis('equal')
        return ax,D_tot,D,fig

    def plot_ice(self,time_plot,vmin=None,vmax=None):
        self.check_self_loaded()
        ind=np.where(self.grid.time_step==time_plot)
        if vmin==None :
            plt.pcolor(self.grid.elons,self.grid.lats,self.ice.ice[ind,:,:].squeeze(),cmap='Blues',vmin=vmin,vmax=vmax)
        else :
            plt.pcolor(self.grid.elons,self.grid.lats,self.ice.ice[ind,:,:].squeeze(),cmap='Blues')
        plt.colorbar()
        plt.axis('equal')
        plt.contour(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),levels=[0])
        plt.show()

    def plot_topo(self,time_plot,vmin=None,vmax=None):
        self.check_self_loaded()
        ind=np.where(self.grid.time_step==time_plot)
        if vmin!=None :
            plt.pcolor(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),cmap='jet',vmin=vmin,vmax=vmax)
        else :
            plt.pcolor(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),cmap='jet')
        plt.colorbar()
        plt.axis('equal')
        plt.contour(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),levels=[0])
        plt.show()

    def plot_at_point(self,C_lm,theta,phi,ax=None):
        y=np.zeros(self.time_step_number)+0j
        if ax==None:
            ax=plt.subplot()
        for i in range(self.time_step_number):
            y[i]=np.real(calc_at_point(C_lm[i,:],self,theta,phi))
        ax.plot(self.grid.time_step[1:],y)
        plt.show()
        return ax,y

    
    def plot_at_point_derivation(self,C_lm,theta,phi,ax=None):
        y=np.zeros(self.time_step_number)+0j
        if ax==None:
            ax=plt.subplot()
        for i in range(1,self.time_step_number):
            y[i]=calc_at_point((C_lm[i,:]-C_lm[i-1,:])/(self.time_step[i]-self.time_step[i-1])/1000,self,theta,phi)
        ax.plot(self.grid.time_step[1:],y)
        plt.show()
        return ax,y

    def plot_3D_at_time(self,data,time_plot,vmin=None,vmax=None):
        self.check_self_loaded()
        if not('P_lm_res' in self.__dict__):
            self.init_resolution_plot()

        
        ind=np.where(self.grid.time_step==time_plot)

        if len(data.shape)>2:
            data_at_time=sphericalobject(data[ind,:,:].squeeze()).grdtocoeff(self)
        elif len(data.shape)==2 :
            data_at_time=sphericalobject(data[ind,:].squeeze(),'coeff')
        else :
            data_at_time=sphericalobject(data,'coeff')
        data_grid=data_at_time.coefftogrd_hd(self)
        if vmin == None :
            min_grid=data_grid.min()
            max_grid=data_grid.max()
            vmin=-np.abs(np.array([min_grid,max_grid])).max()
            vmax=-vmin
        lon,lat=np.meshgrid(self.elons_res, self.lats_res)
        xs, ys, zs = mapping_map_to_sphere(lon,lat)

        topo=self.topo.topo[ind,:,:]
        topo=interp_on(topo,self.grid.elons,self.grid.lats,self.elons_res,self.lats_res)
        cs = plt.contour(self.elons_res, self.lats_res, np.real(topo.squeeze()),levels=[0])

        list_line=[]

        for item in cs.collections:
            x=np.array([])
            y=np.array([])
            for i in item.get_paths():
                v = i.vertices
                x = np.concatenate((x,v[:, 0]))
                y = np.concatenate((y,v[:, 1]))
            xsl,ysl,zsl=mapping_map_to_sphere(np.array(x) , np.array(y))

            list_line.append(go.Scatter3d(x=xsl,y=ysl,z=zsl,mode='lines'))
        




        topo_sphere=dict(type='surface',
            x=xs,
            y=ys,
            z=zs,
            colorscale='rdbu',
            surfacecolor=data_grid,
            cmin=vmin,
            cmax=vmax)
        

        noaxis=dict(showbackground=False,
            showgrid=False,
            showline=False,
            showticklabels=False,
            ticks='',
            title='',
            zeroline=False)

        titlecolor = 'white'
        bgcolor = 'black'
        

        layout = go.Layout(
        autosize=False, width=1200, height=800,
        title = '3D spherical topography map',
        titlefont = dict(family='Courier New', color=titlecolor),
        showlegend = False,
        scene = dict(
            xaxis = noaxis,
            yaxis = noaxis,
            zaxis = noaxis,
            aspectmode='manual',
            aspectratio=go.layout.scene.Aspectratio(
            x=1, y=1, z=1)),
        paper_bgcolor = bgcolor,
        plot_bgcolor = bgcolor)

        list_line.append(topo_sphere)
        plot_data=list_line[-1::-1]
        fig = go.Figure(data=plot_data, layout=layout)
        plot(fig, validate = False, filename='SphericalTopography.html',
        auto_open=True)
