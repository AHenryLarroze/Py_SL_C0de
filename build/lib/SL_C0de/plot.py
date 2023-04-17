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

class world_plot(object):
    def __init__(self):
        ploting='on'

    def init_resolution_plot(self,N=64) :
        self.maxdeg_hd=N
        if N>self.maxdeg :
            self.P_lm_res=legendre(N,self.nb_workers)
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

    def check_model_loaded(self):
        if not('grid' in self.__dict__) : 
            logging.error('error: ', "no model created in this object please use SL_model and SL_model.init_grided_parameters()")
            sys.exit(1)
        
    def plot_RSL(self,time_plot,vmin=-1000,vmax=1000) :
        self.check_model_loaded()
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
    
    def plot_ice(self,time_plot,vmin=None,vmax=None):
        self.check_model_loaded()
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
        self.check_model_loaded()
        ind=np.where(self.grid.time_step==time_plot)
        if vmin==None :
            plt.pcolor(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),cmap='jet',vmin=vmin,vmax=vmax)
        else :
            plt.pcolor(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),cmap='jet')
        plt.colorbar()
        plt.axis('equal')
        plt.contour(self.grid.elons,self.grid.lats,self.topo.topo[ind,:,:].squeeze(),levels=[0])
        plt.show()

    def plot_RSL_3D(self,time_plot):
        self.check_model_loaded()
        if not('RSL' in self.SL.__dict__):
            self.calc_RSL()
        if not('P_lm_res' in self.__dict__):
            self.init_resolution_plot()
        # cmap=plt.get_cmap('seismic')
        # Ctopo=cmap(np.linspace(0,1,200))[:,:3]
        cmin=-6
        cmax=6
        ind=np.where(self.grid.time_step==time_plot)

        RSL_at_time=sphericalobject(self.SL.RSL[ind,:,:].squeeze()).grdtocoeff(self)
        RSL_grid=RSL_at_time.coefftogrd_hd(self)
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
            surfacecolor=RSL_grid,
            cmin=cmin,
            cmax=cmax)
        

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
