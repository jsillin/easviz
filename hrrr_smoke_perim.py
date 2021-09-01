import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeatres
import pandas as pd
#import wget
import os.path
from os import path
import sys
import geopandas as gpd
import supplementary_tools as spt
#from google.cloud import storage

mdate, runtime = spt.get_init_time('HRRR')
print(runtime)
#storage_client = storage.Client()
#noaa_hrrr_bucket = storage_client.bucket("high-resolution-rapid-refresh")
#blob = bucket.blob(

per = gpd.read_file('https://opendata.arcgis.com/datasets/2191f997056547bd9dc530ab9866ab61_0.geojson')

#mdate = spt.get_init_hour('HRRR')

def plotsmoke(ds,fhr,slev):
    mden = ds['unknown']
    dtfs = pd.to_datetime(mden.step.values+mden.time.values).strftime('%Y-%m-%d_%H%MZ') 
    fig = plt.figure(figsize=(12,12))
    prj = projection=ccrs.PlateCarree()
    ax = fig.add_subplot(111,projection=prj)
    ax.add_feature(cfeatres.BORDERS,edgecolor='lemonchiffon')
    ax.add_feature(cfeatres.COASTLINE,edgecolor='lemonchiffon')
    ax.add_feature(cfeatres.STATES,edgecolor='lemonchiffon')

    if slev=='ts':
        mdc = ax.contourf(mden.longitude,mden.latitude,mden*1000,levels=range(0,750,10),extend='max',cmap='bone',transform=ccrs.PlateCarree())
        cbr = fig.colorbar(mdc, orientation = 'vertical', ax = ax, extendrect=False, pad=0.01, ticks = range(0,750,50), shrink = 0.55)
    elif slev=='ns':
        mdc = ax.contourf(mden.longitude,mden.latitude,mden,levels=range(0,300,10),extend='max',cmap='bone',transform=ccrs.PlateCarree())
        cbr = fig.colorbar(mdc,orientation='vertical',ax=ax,extendrect=False,pad=0.01,ticks=range(0,300,10),shrink=0.55)
    
    cbr.set_label('mg/m^2')

    if slev =='ts':
        ax.set_title('Total Atmosphere Smoke and Active Fire Perimeters')
    elif slev=='ns':
        ax.set_title('Near-Surface Smoke and Active Fire Perimeters')
    ax.set_title(dtfs,loc='right')
    ax.set_title('HRRR',loc='left')    
    per.plot(ax=ax,edgecolor='red')
    
    ax.set_extent((-130,-100,30,49))
    plt.savefig('WC_'+slev+'_perim_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-125,-110,40,49.5))
    plt.savefig('NW_'+slev+'_perim_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-100,-65,30,49))
    plt.savefig('EC_'+slev+'_perim_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-130,-65,24,50))
    plt.savefig('US_'+slev+'_perim_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)

for i in range(1,49):
    if i <10:
        phr = '0'+str(i)
    else:
        phr = str(i)    
    nsds = xr.open_dataset('/home/jhs389/plotting/hrrr.t'+runtime+'z.wrfnatf'+phr+'.grib2',backend_kwargs=dict(filter_by_keys={'typeOfLevel':'heightAboveGround','level':8}))
    tsds = xr.open_dataset('/home/jhs389/plotting/hrrr.t'+runtime+'z.wrfprsf'+phr+'.grib2',filter_by_keys={'typeOfLevel':'atmosphereSingleLayer','stepType':'instant'})
   
    plotsmoke(nsds,i,'ns')
    plotsmoke(tsds,i,'ts')
