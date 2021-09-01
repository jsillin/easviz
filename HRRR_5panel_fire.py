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
#import wget
import os.path
from os import path
import sys
import geopandas as gpd
import metpy
import supplementary_tools as spt

mdate, runtime = spt.get_init_time('HRRR')
print(runtime)

def hainesf(fname):
    ds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    ods = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'surface'})
    
    elev = ods['orog']*3.28084
    
    t9 = ds['t'].sel(isobaricInhPa=950)
    t8 = ds['t'].sel(isobaricInhPa=850)
    t7 = ds['t'].sel(isobaricInhPa=700)
    t5 = ds['t'].sel(isobaricInhPa=500)

    stability_low_raw = t9-t8
    stability_mid_raw = t8-t7
    stability_hi_raw = t7-t5

    stability_low1 = stability_low_raw.where(stability_low_raw>3,-1)
    stability_low2 = stability_low1.where(stability_low1<8,-3)
    stability_low3 = stability_low2.where(stability_low2<0,2)
    stability_low = abs(stability_low3)

    stability_mid1 = stability_mid_raw.where(stability_mid_raw>5,-1)
    stability_mid2 = stability_mid1.where(stability_mid1<11,-3)
    stability_mid3 = stability_mid2.where(stability_mid2<0,2)
    stability_mid = abs(stability_mid3)

    stability_hi1 = stability_hi_raw.where(stability_hi_raw>3,-1)
    stability_hi2 = stability_hi1.where(stability_hi1<8,-3)
    stability_hi3 = stability_hi2.where(stability_hi2<0,2)
    stability_high = abs(stability_hi3)

    td8 = ds['dpt'].sel(isobaricInhPa=850)
    td7 = ds['dpt'].sel(isobaricInhPa=700)

    moisture_low_raw = t8-td8
    moisture_mid_raw = t8-td8
    moisture_hi_raw = t7-td7

    moisture_low1 = moisture_low_raw.where(moisture_low_raw>5,-1)
    moisture_low2 = moisture_low1.where(moisture_low1<10,-3)
    moisture_low3 = moisture_low2.where(moisture_low2<0,2)
    moisture_low = abs(moisture_low3)

    moisture_mid1 = moisture_mid_raw.where(moisture_mid_raw>5,-1)
    moisture_mid2 = moisture_mid1.where(moisture_mid1<13,-3)
    moisture_mid3 = moisture_mid2.where(moisture_mid2<0,2)
    moisture_mid = abs(moisture_mid3)

    moisture_hi1 = moisture_hi_raw.where(moisture_hi_raw>14,-1)
    moisture_hi2 = moisture_hi1.where(moisture_hi1<21,-3)
    moisture_hi3 = moisture_hi2.where(moisture_hi2<0,2)
    moisture_high = abs(moisture_hi3)

    haines_low = moisture_low+stability_low
    haines_mid = moisture_mid+stability_mid
    haines_high = moisture_high+stability_high

    low_elev_haines = haines_low.where(elev<1000,0)
    high_elev_haines = haines_high.where(elev>3000,0)
    haines_nomid = high_elev_haines+low_elev_haines
    haines = haines_nomid.where(haines_nomid!=0,haines_mid)
    
    return haines


for i in range(1,49):
    if i<10:
        phr='0'+str(i)
    else:
        phr=str(i)

    natf = '/home/jhs389/plotting/hrrr.t'+runtime+'z.wrfnatf'+phr+'.grib2'
    prsf = '/home/jhs389/plotting/hrrr.t'+runtime+'z.wrfprsf'+phr+'.grib2'

    ns_ds = xr.open_dataset(natf,engine='cfgrib',backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'heightAboveGround','level':8}))
    wd_ds = xr.open_dataset(natf,engine='cfgrib',backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'heightAboveGround','level':10}))
    t2_ds = xr.open_dataset(natf,engine='cfgrib',backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'heightAboveGround','level':2}))
    mlrds = xr.open_dataset(natf,engine='cfgrib',backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'heightAboveGround','level':4000}))
    llrds = xr.open_dataset(natf,engine='cfgrib',backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'heightAboveGround','level':1000}))
    cpds = xr.open_dataset(prsf,engine='cfgrib',filter_by_keys={'typeOfLevel': 'pressureFromGroundLayer', 'shortName':'cape'})
    rf_ds = xr.open_dataset(natf,engine='cfgrib',filter_by_keys={'stepType': 'instant', 'typeOfLevel': 'heightAboveGround'})
    ta_ds = xr.open_dataset(natf,engine='cfgrib',filter_by_keys={'typeOfLevel': 'atmosphereSingleLayer','stepType': 'instant'})
    lgds = xr.open_dataset(prsf,engine='cfgrib',filter_by_keys={'stepType':'instant','typeOfLevel':'atmosphere'})

    rh = t2_ds['r2']
    mlref = mlrds['refd']
    llref = rf_ds['refd'].isel(heightAboveGround=0)
    diff = mlref-llref
    lgt = lgds['ltng']
    ref = rf_ds['refd'].isel(heightAboveGround=0)
    ns_mden = ns_ds['unknown']
    ts_mden = ta_ds['unknown']
    u10 = wd_ds['u10']*1.94384
    v10 = wd_ds['v10']*1.94384
    mucape = cpds.isel(pressureFromGroundLayer=2).to_array().squeeze()

    per = gpd.read_file('https://opendata.arcgis.com/datasets/2191f997056547bd9dc530ab9866ab61_0.geojson')

    haines = hainesf(prsf)
    
    dtfs = str(t2_ds.valid_time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    wind_slice = slice(18,-18,18)
    x, y = u10.longitude,u10.latitude

    norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
    newcmap = ListedColormap(cmap_ref(range(40, 194)))
    newnorm = Normalize(0,77)

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(44,15))
    zH5_crs=ccrs.PlateCarree()
    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
    gs.update(hspace=0.1,wspace=0.01)
    ax1 = fig.add_subplot(gs[:, 1], projection = zH5_crs)
    ax2 = fig.add_subplot(gs[0, 0], projection = zH5_crs)
    ax3 = fig.add_subplot(gs[1, 0], projection = zH5_crs)
    ax4 = fig.add_subplot(gs[0, 2], projection = zH5_crs)
    ax5 = fig.add_subplot(gs[1, 2], projection = zH5_crs)

    ax1.coastlines(resolution='10m',edgecolor='lemonchiffon')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'),edgecolor='lemonchiffon')
    ax1.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='lemonchiffon')

    ax2.coastlines(resolution='10m',edgecolor='lemonchiffon')
    ax2.add_feature(cfeature.BORDERS.with_scale('10m'),edgecolor='lemonchiffon')
    ax2.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='lemonchiffon')

    ax3.coastlines(resolution='10m',edgecolor='lemonchiffon')
    ax3.add_feature(cfeature.BORDERS.with_scale('10m'),edgecolor='lemonchiffon')
    ax3.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='lemonchiffon')

    ax4.coastlines(resolution='10m',edgecolor='lemonchiffon')
    ax4.add_feature(cfeature.BORDERS.with_scale('10m'),edgecolor='lemonchiffon')
    ax4.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='lemonchiffon')

    ax5.coastlines(resolution='10m',edgecolor='lemonchiffon')
    ax5.add_feature(cfeature.BORDERS.with_scale('10m'),edgecolor='lemonchiffon')
    ax5.add_feature(cfeature.STATES.with_scale('10m'),edgecolor='lemonchiffon')

    ax1.contourf(ts_mden.longitude,ts_mden.latitude,ts_mden*1000,levels=range(0,750,10),extend='max',cmap='bone',transform=ccrs.PlateCarree())
    refc = ax1.contourf(ref.longitude,ref.latitude,ref,norm=newnorm,cmap=newcmap,levels=range(5,75,1),alpha=0.7)
    per.boundary.plot(ax=ax1,edgecolor='red')

    ax1.barbs(x[wind_slice,wind_slice],y[wind_slice,wind_slice],u10[wind_slice,wind_slice],v10[wind_slice,wind_slice], length=6,color='yellow')

    ax2.contourf(ns_mden.longitude,ns_mden.latitude,ns_mden,levels=range(0,250,10),extend='max',cmap='gnuplot2',transform=ccrs.PlateCarree())
    per.boundary.plot(ax=ax2,edgecolor='red')

    ax3.contourf(rh.longitude,rh.latitude,rh,levels=range(0,100,2),cmap='BrBG')
    per.boundary.plot(ax=ax3,edgecolor='red')

    ax4.contourf(haines.longitude,haines.latitude,haines,levels=[2,3,4,5,6],colors=['darkgreen','forestgreen','lightgreen','gold','darkorange','magenta'],transform=ccrs.PlateCarree())
    per.boundary.plot(ax=ax4,edgecolor='red')

    ax5.contourf(mucape.longitude,mucape.latitude,mucape,levels=range(0,3000,100),cmap='YlOrBr',transform=ccrs.PlateCarree())
    ax5.contourf(diff.longitude,diff.latitude,diff,levels=range(-50,50,2),cmap='PuOr_r',transform=ccrs.PlateCarree())
    ax5.contourf(lgt.longitude,lgt.latitude,lgt,levels=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5],cmap='spring',transform=ccrs.PlateCarree())
    per.boundary.plot(ax=ax5,edgecolor='red')

    ax1.set_extent((-127,-105,40,49.5))
    ax2.set_extent((-125,-107,41.3,49.5))
    ax3.set_extent((-125,-107,41.3,49.5))
    ax4.set_extent((-125,-107,41.3,49.5))
    ax5.set_extent((-125,-107,41.3,49.5))

    ax1.set_title('HRRR Fire Weather Forecasting Dashboard')
    ax1.set_title('\n Valid: '+t2_ds.valid_time.dt.strftime('%a %b %d %HZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n HRRR Init: '+t2_ds.time.dt.strftime('%Y-%m-%d %HZ').item(),fontsize=11,loc='left')

    ax2.set_title('Near-Surface Smoke')
    ax3.set_title('2m Relative Humidity')
    ax4.set_title('Haines Index')
    ax5.set_title('Thunderstorms')
    plt.savefig('NW_5panel_fire_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
