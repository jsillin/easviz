# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 04:36:27 2021

@author: jacks
"""
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as col
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
from scipy.ndimage import gaussian_filter
import scipy.ndimage
import numpy as np
import metpy.calc as mpcalc
import supplementary_tools as spt 

rundate, runtime = spt.get_init_time('GFS')

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
    cbar = plt.colorbar(im, cax=cax, ticks=range(-2,7,1),
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

def pad_fhr(fhr):
    fhr = 3*i
    if fhr <10:
        pfhr = '00'+str(fhr)
    elif fhr <100:
        pfhr = '0'+str(fhr)
    else:
        pfhr = str(fhr)
    return pfhr

for i in range(0,80):
    fhr = i*3
    pfhr = pad_fhr(i)
    fname = '/home/jhs389/plotting/gfs.t'+runtime+'z.pgrb2.0p25.f'+pfhr
    
    sst_mdate = spt.get_init_time('SST')[0]
    sst_myr = sst_mdate[0:6]
    
    ds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    rfds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'heightAboveGround','level':1000})
    mslpds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'meanSea'})
    tmwds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'heightAboveGround','level':10})
    pvds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'potentialVorticity'})
    sstds = xr.open_dataset('https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/V2.1/AVHRR/'+sst_myr+'/oisst-avhrr-v02r01.'+sst_mdate+'_preliminary.nc')
    
    ref = rfds['refd']
    u10 = tmwds['u10']*1.94384
    v10 = tmwds['v10']*1.94384
    mslp = mslpds['prmsl']/100
    pvds = pvds.isel(potentialVorticity=1)
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
    vo8 = ds['absv'].sel(isobaricInhPa=850)*(10**5)
    mlrh = rh.sel(isobaricInhPa=[700,650,600,550,500,450,400]).mean(dim='isobaricInhPa')
    mlshr_u = u4-u7
    mlshr_v = v4-v7
    
    wind_slice = slice(8,-8,8)
    smwind_slice = slice(15,-15,15)
    x, y = u2.longitude,u2.latitude
    
    pv = mpcalc.potential_vorticity_barotropic(gph.sel(isobaricInhPa=200),u.sel(isobaricInhPa=200),v.sel(isobaricInhPa=200))
    pvt = pv
    
    norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
    newcmap = ListedColormap(cmap_ref(range(40, 194)))
    newnorm = Normalize(0,77)
    
    pvs = gaussian_filter(pvt*10**8,sigma=1,order=0)
    
    sst = sstds['sst'].squeeze()
    td_steer_u = u.sel(isobaricInhPa=[850,800,750,700]).mean(dim='isobaricInhPa')
    td_steer_v = v.sel(isobaricInhPa=[850,800,750,700]).mean(dim='isobaricInhPa')
    hu_steer_u = u.sel(isobaricInhPa=[850,800,750,700,650,600,550,500,450,400,350]).mean(dim='isobaricInhPa')
    hu_steer_v = v.sel(isobaricInhPa=[850,800,750,700,650,600,550,500,450,400,350]).mean(dim='isobaricInhPa')
    
    fig = plt.figure(figsize=(44,15))
    zH5_crs=ccrs.PlateCarree()
    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
    gs.update(hspace=0.1,wspace=0.01)
    ax1 = fig.add_subplot(gs[:, 1], projection = zH5_crs)
    ax2 = fig.add_subplot(gs[0, 0], projection = zH5_crs)
    ax3 = fig.add_subplot(gs[1, 0], projection = zH5_crs)
    ax4 = fig.add_subplot(gs[0, 2], projection = zH5_crs)
    ax5 = fig.add_subplot(gs[1, 2], projection = zH5_crs)
    
    ax1.coastlines(resolution='50m',edgecolor='gray')
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
    ax1.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray')
    
    ax2.coastlines(resolution='50m',edgecolor='gray')
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
    #ax2.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='gray')
    
    ax3.coastlines(resolution='50m',edgecolor='gray')
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
    #ax3.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='gray')
    
    ax4.coastlines(resolution='50m',edgecolor='gray')
    ax4.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
    #ax4.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='gray')
    
    ax5.coastlines(resolution='50m',edgecolor='gray')
    ax5.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray')
    #ax5.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='gray')
    
    ax1.set_extent((-105,-5,0,50))
    ax2.set_extent((-100,-5,0,50))
    ax3.set_extent((-100,-5,0,50))
    ax4.set_extent((-100,-5,0,50))
    ax5.set_extent((-100,-5,0,50))
    
    vorc = ax1.contourf(vo8.longitude,vo8.latitude,vo8,levels=range(5,50,1),extend='both',cmap='Wistia',transform=zH5_crs)
    ax1.contour(h5.longitude,h5.latitude,h5,levels=range(4800,6000,60),colors='black',transform=zH5_crs)
    ax1.barbs(x[wind_slice],y[wind_slice],u2[wind_slice,wind_slice],v2[wind_slice,wind_slice], length=5,color='blue')
    ax1.contour(mslp.longitude,mslp.latitude,mslp,levels=range(900,1080,4),colors='dimgray',linewidths=2)
    ax1.barbs(x[wind_slice],y[wind_slice],u10[wind_slice,wind_slice],v10[wind_slice,wind_slice], length=5,color='purple')
    refc = ax1.contourf(ref.longitude,ref.latitude,ref,norm=newnorm,cmap=newcmap,levels=range(5,75,1),alpha=0.4)
    ax1.contour(sst.lon,sst.lat,sst,levels=[26,28,30],colors=['orange','red','magenta'],linewidths=2.5)
    
    relc = ax2.contourf(mlrh.longitude,mlrh.latitude,mlrh,levels=range(0,102,2),cmap='BrBG')
    ax2.quiver(x[smwind_slice],y[smwind_slice],mlshr_u[smwind_slice,smwind_slice],mlshr_v[smwind_slice,smwind_slice],color='red',scale=900)
    		   
    pvoc = ax3.contourf(pvt.longitude,pvt.latitude,pvs,levels=np.linspace(-2,6,16),extend='both',cmap='tab20b')
    ax3.contour(mslp.longitude,mslp.latitude,mslp,levels=range(900,1012,4),colors='white',linewidths=2)
    5
    sstc = ax4.contourf(sst.lon,sst.lat,sst,levels=range(10,32,1),cmap='coolwarm',extend='both')
    ax4.contour(sst.lon,sst.lat,sst,levels=[26.5],colors=['yellow'],linewidths=2.5)
    ax4.contour(mslp.longitude,mslp.latitude,mslp,levels=range(900,1012,4),colors='white',linewidths=2)
    
    addvortcolorbar(ax1,fig,vorc,range(0,50,5))
    addrefcolorbar(ax1,fig,refc,range(5,75,1))
    addrhcolorbar(ax2,fig,relc,range(0,102,2))
    addpvucolorbar(ax3,fig,pvoc,np.linspace(-2,6,16))
    addsstcolorbar(ax4,fig,sstc,range(10,32,1))
    
    '''cbar2 = fig.colorbar(relc,orientation='horizontal',pad=0.01,ax=ax2,aspect=50,extendrect=False, ticks=np.arange(1,105,5),shrink=0.85)
    cbar2.set_label('%',fontsize=14)
    
    cbar3 = fig.colorbar(pvoc,orientation='horizontal',pad=0.01,ax=ax3,aspect=50,extendrect=False, ticks=[-2,-1,0,1,2,3,4,5,6,7],shrink=0.85)
    cbar3.set_label('PVU',fontsize=14)
    
    cbar4 = fig.colorbar(sstc,orientation='horizontal',pad=0.01,ax=ax4,aspect=50,extendrect=False, ticks=range(10,35,5),shrink=0.85)
    cbar4.set_label('C',fontsize=14)
    
    #cbar5 = fig.colorbar(difc,orientation='vertical',pad=0.01,ax=ax5,aspect=50,extendrect=False, ticks=np.arange(-50,50,5),shrink=0.9)
    #cbar5.set_label('dBZ',fontsize=14)
    '''
    ax5.quiver(x[smwind_slice],y[smwind_slice],td_steer_u[smwind_slice,smwind_slice],td_steer_v[smwind_slice,smwind_slice],color='red',scale=500)
    ax5.quiver(x[smwind_slice],y[smwind_slice],hu_steer_u[smwind_slice,smwind_slice],hu_steer_v[smwind_slice,smwind_slice],color='magenta',scale=500)
    
    ax1.set_title('GFS Tropical Cyclone Forecasting Dashboard')
    ax1.set_title('\n Valid: '+pvds.valid_time.dt.strftime('%a %b %d %HZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GFS Init: '+pvds.time.dt.strftime('%Y-%m-%d %HZ').item(),fontsize=11,loc='left')
    
    ax2.set_title('700-400mb RH and Shear')
    ax3.set_title('200mb Potential Vorticity and MSLP <1012mb')
    ax4.set_title('Sea Surface Temperatures and MSLP <1012mb')
    ax5.set_title('Steering Winds')
    
    orange = mpatches.Patch(color='orange', label='SST = 26C')
    red = mpatches.Patch(color='red', label='SST = 28C')
    pink = mpatches.Patch(color='magenta', label='SST = 30C')
    black = mpatches.Patch(color='black', label='500mb Heights (m)')
    gray = mpatches.Patch(color='dimgray', label='MSLP (mb)')
    blue = mpatches.Patch(color='blue', label='200mb Wind (kts)')
    white = mpatches.Patch(color='purple', label='10m Wind (kts)')
    
    leg = ax1.legend(handles=[orange,red,pink,blue,white,black,gray],loc=4,title='Overlay Legend',framealpha=1)
    leg.set_zorder(100)
    
    red1 = mpatches.Patch(color='red', label='Weaker Storm')
    pink1 = mpatches.Patch(color='magenta', label='Stronger Storm')
    leg = ax5.legend(handles=[red1,pink1],loc=4,title='Steering',framealpha=1)
    leg.set_zorder(100)
    
    yellow = mpatches.Patch(color='yellow', label='26C Isotherm')
    leg = ax5.legend(handles=[red1,pink1],loc=4,title='Steering',framealpha=1)
    leg.set_zorder(100)
    
    plt.savefig('natl_hurr_dash_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax1.set_extent((-98,-55,10,35))
    ax2.set_extent((-98,-55,10,37))
    ax3.set_extent((-98,-55,10,37))
    ax4.set_extent((-98,-55,10,37))
    ax5.set_extent((-98,-55,10,37))
    plt.savefig('watl_hurr_dash_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
    
