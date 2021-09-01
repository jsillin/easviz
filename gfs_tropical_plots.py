# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 09:04:50 2021

@author: jacks
"""

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as col
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import os.path
from os import path
import sys
import geopandas as gpd
import metpy
import hrrr_dash_cbars as cbars
import matplotlib.patches as mpatches
import numpy as np
from scipy.ndimage import gaussian_filter
import scipy.ndimage
import metpy.calc as mpcalc
import supplementary_tools as spt

def addvortcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for CAPE

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.195
    bottom = axes_bbox.y0 - 0.011
    width = 0.185
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs,
                        orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Absolute Vorticity (10^5 s^-1)', size=8)  # MODIFY THIS for other fields!!
    
def addrefcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0
    bottom = axes_bbox.y0- 0.011
    width = 0.185
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=range(5,80,5),
                        orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Composite Reflectivity (dBZ)', size=8)  # MODIFY THIS for other fields!!
    
def addpvucolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.05*abs(axes_bbox.x1-axes_bbox.x0)
    bottom = axes_bbox.y0+0.002
    width = 0.9*abs(axes_bbox.x1-axes_bbox.x0)
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=range(-3,7,1),
                        orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Potential Vorticity Units', size=8)  # MODIFY THIS for other fields!!
    
def addsstcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0+ 0.05*abs(axes_bbox.x1-axes_bbox.x0)
    bottom = axes_bbox.y0+0.002
    width = 0.9*abs(axes_bbox.x1-axes_bbox.x0)
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=range(10,35,2),
                        orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    #cbar.set_label('Degrees Celsius', size=8)  # MODIFY THIS for other fields!!
    
def addrhcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0+ 0.05*abs(axes_bbox.x1-axes_bbox.x0)
    bottom = axes_bbox.y0+0.002
    width = 0.9*abs(axes_bbox.x1-axes_bbox.x0)
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=np.arange(1,105,5),
                        orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    #cbar.set_label('Degrees Celsius', size=8)  # MODIFY THIS for other fields!!
   
   
rundate, runtime = spt.get_init_time('GFS')
   
sst_mdate = spt.get_init_time('SST')[0]
sst_myr = sst_mdate[0:6]

for i in range(0,81):
    fhr = 3*i
    print(fhr)
    if fhr <10:
        pfhr = '00'+str(fhr)
    elif fhr <100:
        pfhr = '0'+str(fhr)
    else:
        pfhr = str(fhr)
    
    fname = 'gfs.t'+runtime+'z.pgrb2.0p25.f'+pfhr
    #fname = 'gfs.t00z.pgrb2.0p25.f240'
       
    ds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    mslpds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'meanSea'})
    tmwds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'heightAboveGround','level':10})
    sstds = xr.open_dataset('https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/V2.1/AVHRR/'+sst_myr+'/oisst-avhrr-v02r01.'+sst_mdate+'_preliminary.nc')
    
    u10 = tmwds['u10']*1.94384
    v10 = tmwds['v10']*1.94384
    mslp = mslpds['prmsl']/100
    gph = ds['gh']
    u = ds['u']
    v = ds['v']
    h5 = ds['gh'].sel(isobaricInhPa=500)
    u2 = ds['u'].sel(isobaricInhPa=200)
    v2 = ds['v'].sel(isobaricInhPa=200)
    u7 = ds['u'].sel(isobaricInhPa=700)
    v7 = ds['u'].sel(isobaricInhPa=700)
    u4 = ds['u'].sel(isobaricInhPa=400)
    v4 = ds['v'].sel(isobaricInhPa=400)
    rh = ds['r']
    sst = sstds['sst'].squeeze()
    td_steer_u = u.sel(isobaricInhPa=[850,800,750,700]).mean(dim='isobaricInhPa')
    td_steer_v = v.sel(isobaricInhPa=[850,800,750,700]).mean(dim='isobaricInhPa')
    hu_steer_u = u.sel(isobaricInhPa=[850,800,750,700,650,600,550,500,450,400,350]).mean(dim='isobaricInhPa')
    hu_steer_v = v.sel(isobaricInhPa=[850,800,750,700,650,600,550,500,450,400,350]).mean(dim='isobaricInhPa')
    mlrh = rh.sel(isobaricInhPa=[700,650,600,550,500,450,400]).mean(dim='isobaricInhPa')
    mlshr_u = u4-u7
    mlshr_v = v4-v7
    
    wind_slice = slice(8,-8,8)
    smwind_slice = slice(15,-15,15)
    x, y = u2.longitude,u2.latitude
    
    
    natl_dom = (-105,-5,0,50)
    watl_dom = (-98,-55,10,37)
    
    if sys.argv[1]=="pv":
        pvds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'potentialVorticity'})
        pvds = pvds.isel(potentialVorticity=1)
        
        pvlev = 300
        pv = pvt =  mpcalc.potential_vorticity_barotropic(gph.sel(isobaricInhPa=pvlev),u.sel(isobaricInhPa=pvlev),v.sel(isobaricInhPa=pvlev))
        pvs = gaussian_filter(pvt*10**8,sigma=2,order=0)
    
        norm_pv, cmap_pv = ctables.registry.get_with_steps('WVCIMSS', -80., .5)
        pvcmap = ListedColormap(cmap_pv(range(0,117)))
        pvnorm = Normalize(-2.5,4)    
    
        fig = plt.figure(figsize=(15,15))
        prj = ccrs.PlateCarree()
        ax = fig.add_subplot(111,projection=prj)
        ax.set_extent(natl_dom)
        
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'),edgecolor='gray')
        ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
        ax.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
        
        pvoc = ax.contourf(pvt.longitude,pvt.latitude,pvs,levels=np.linspace(-3,4,60),extend='both',cmap=pvcmap,norm=pvnorm)
        addpvucolorbar(ax,fig,pvoc,np.linspace(-2,6,16))
        ax.contour(mslp.longitude,mslp.latitude,mslp,levels=range(900,1012,4),colors='white',linewidths=2)
        ax.set_title('\n Valid: '+pvds.valid_time.dt.strftime('%a %b %d %HZ').item(),fontsize=11,loc='right')
        ax.set_title('\n GFS Init: '+pvds.time.dt.strftime('%Y-%m-%d %HZ').item(),fontsize=11,loc='left')
        ax.set_title(str(pvlev)+'mb Potential Vorticity and MSLP <1012mb')
        plt.savefig('natl_pv_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
        ax.set_extent(watl_dom)
        plt.savefig('watl_pv_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
    
    elif sys.argv[1]=='ov':
        rfds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'heightAboveGround','level':1000})
        prj = ccrs.PlateCarree()        
        ref = rfds['refd']
        vo8 = ds['absv'].sel(isobaricInhPa=850)*(10**5)
        
        norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
        newcmap = ListedColormap(cmap_ref(range(40, 194)))
        newnorm = Normalize(0,77)
    
        fig1 = plt.figure(figsize=(15,15))
        ax1 = fig1.add_subplot(projection=prj)
        
        ax1.coastlines(resolution='50m',edgecolor='gray')
        ax1.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
        ax1.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
        ax1.set_extent(natl_dom)
        
        vorc = ax1.contourf(vo8.longitude,vo8.latitude,vo8,levels=range(5,50,1),extend='both',cmap='Wistia',transform=prj)
        ax1.contour(h5.longitude,h5.latitude,h5,levels=range(4800,6000,60),colors='black',transform=prj)
        ax1.barbs(x[wind_slice],y[wind_slice],u2[wind_slice,wind_slice],v2[wind_slice,wind_slice], length=5,color='blue')
        ax1.contour(mslp.longitude,mslp.latitude,mslp,levels=range(900,1080,4),colors='dimgray',linewidths=2)
        ax1.barbs(x[wind_slice],y[wind_slice],u10[wind_slice,wind_slice],v10[wind_slice,wind_slice], length=5,color='purple')
        refc = ax1.contourf(ref.longitude,ref.latitude,ref,norm=newnorm,cmap=newcmap,levels=range(5,75,1),alpha=0.4)
        ax1.contour(sst.lon,sst.lat,sst,levels=[26,28,30],colors=['orange','red','magenta'],linewidths=2.5)
        plt.savefig('natl_ov_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
        ax1.set_extent(watl_dom)
        plt.savefig('watl_ov_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
    
    elif sys.argv[1]=='rh':
        prj = ccrs.PlateCarree()
        fig2 = plt.figure(figsize=(15,15))
        ax2 = fig2.add_subplot(projection=prj) 
        ax2.set_extent(natl_dom)
        ax2.coastlines(resolution='50m',edgecolor='gray')
        ax2.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
        ax2.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
        relc = ax2.contourf(mlrh.longitude,mlrh.latitude,mlrh,levels=range(0,102,2),cmap='BrBG')
        ax2.barbs(x[smwind_slice],y[smwind_slice],mlshr_u[smwind_slice,smwind_slice],mlshr_v[smwind_slice,smwind_slice],color='red',length=5)
        plt.savefig('natl_rh_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
        ax2.set_extent(watl_dom)
        plt.savefig('watl_rh_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
    
    elif sys.argv[1]=='sst':
        fig4 = plt.figure(figsize=(15,15))
        prj = ccrs.PlateCarree()
        ax4 = fig4.add_subplot(projection=prj) 
        ax4.set_extent(natl_dom)
        ax4.coastlines(resolution='50m',edgecolor='gray')
        ax4.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
        ax4.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
        sstc = ax4.contourf(sst.lon,sst.lat,sst,levels=range(10,32,1),cmap='coolwarm',extend='both')
        ax4.contour(sst.lon,sst.lat,sst,levels=[26.5],colors=['yellow'],linewidths=2.5)
        ax4.contour(mslp.longitude,mslp.latitude,mslp,levels=range(900,1012,4),colors='white',linewidths=2)
        plt.savefig('natl_sst_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
        ax4.set_extent(watl_dom)
        plt.savefig('watl_sst_'+fhr+'.png',bbox_inches='tight',pad_inches=0.1)
