# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 05:01:50 2021

@author: jacks
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
from matplotlib.path import Path
from cartopy.mpl.patch import geos_to_path
from datetime import datetime, timedelta
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.colors as col
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import supplementary_tools as spt 

def load_basemap(map_path):
    img=plt.imread(map_path)
    img_proj = ccrs.PlateCarree()
    img_extent = (-180, 180, -90, 90)
    return img, img_proj,img_extent
	
map_path='/home/jhs389/plotting/RenderData.tif'
bm_img,bm_proj,bm_extent=load_basemap(map_path)

rundate, runtime = spt.get_init_time('GFS')

for i in range(0,120):
	fhr = i
	if fhr <10:
		pfhr = '00'+str(fhr)
	elif fhr <100:
		pfhr = '0'+str(fhr)
	else:
		pfhr = str(fhr)
		
	fname = '/home/jhs389/plotting/gfs.t'+runtime+'z.pgrb2.0p25.f'+pfhr
	
	ds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
	t8 = ds['t'].sel(isobaricInhPa=850)-273.15	
	
	# Define the colormaps
	norm_ref, cmap_ref = ctables.registry.get_with_steps('rainbow', -80., .5)
	newcmap = ListedColormap(cmap_ref(range(0, 175)))
	newnorm = Normalize(-44,48)

	dtfs = pd.to_datetime(t8.time.values+t8.step.values).strftime('%Y-%m-%d_%H%MZ') 
	fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw={'projection': ccrs.NearsidePerspective(central_longitude=-76.5,central_latitude=42.5),'transform':ccrs.PlateCarree()})
	ax.add_feature(cfeature.BORDERS,edgecolor='k')
	ax.add_feature(cfeature.COASTLINE,edgecolor='k')
	ax.add_feature(cfeature.STATES,edgecolor='k')

	t8c = ax.contourf(t8.longitude,t8.latitude,t8,norm=newnorm,cmap=newcmap,antialiased=True,levels=range(-40,40,1),alpha=0.7,transform=ccrs.PlateCarree())
	ax.contour(t8.longitude,t8.latitude,t8,levels=[0],colors=['blue'],transform=ccrs.PlateCarree())

	cbar = fig.colorbar(t8c,orientation='vertical',pad=0.01,ax=ax,aspect=50,extendrect=False,ticks=range(-40,40,5),shrink=0.65)
	ax.set_extent((-157,0,0,90))

	ax.set_title('850hPa Temperature (C)',fontsize=16)
	ax.set_title('Valid: '+dtfs,loc='right',fontsize=11)
	ax.set_title('GFS Forecast via NOAA \nPlot by Jack Sillin',loc='left',fontsize=11)
	plt.savefig('nwh_gfs_t8temp_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
