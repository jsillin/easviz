# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 04:55:40 2021

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
	
	ds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'atmosphereSingleLayer'})
	pwat = ds['pwat']*0.0393701
	
	# Define the colormaps
	colors = [(0,1,1,c) for c in np.linspace(0.1,.6,80)]
	cmappwat = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=8)	
	dtfs = pd.to_datetime(pwat.time.values+pwat.step.values).strftime('%Y-%m-%d_%H%MZ') 
	fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw={'projection': ccrs.Robinson(),'transform':ccrs.PlateCarree()})
	ax.imshow(bm_img,extent=bm_extent,transform=bm_proj,origin='upper')
	ax.contourf(pwat.longitude,pwat.latitude,pwat,levels=np.linspace(0.1,2.5,50),antialiased=True,cmap=cmappwat,extend='max',transform=ccrs.PlateCarree())
	ax.set_title('Total Column Water',fontsize=16)
	ax.set_title('Valid: '+dtfs,loc='right',fontsize=11)
	ax.set_title('GFS Forecast via NOAA \nPlot by Jack Sillin',loc='left',fontsize=11)
	plt.savefig('global_gfs_pwat_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
