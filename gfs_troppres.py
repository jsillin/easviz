# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 05:03:20 2021

@author: jacks
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.colors import Normalize

import pandas as pd
import xarray as xr
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import supplementary_tools as spt
import matplotlib.animation as manim

rundate, runtime = spt.get_init_time('GFS')
animfig = plt.figure(figsize=(15,15))
ims = []

for i in range(0,120):
	fhr = i
	if fhr <10:
		pfhr = '00'+str(fhr)
	elif fhr <100:
		pfhr = '0'+str(fhr)
	else:
		pfhr = str(fhr)
		
	fname = 'gfs.t'+runtime+'z.pgrb2.0p25.f'+pfhr
	
	pvds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'potentialVorticity'})
	dtfs = pd.to_datetime(pvds.time.values+pvds.step.values).strftime('%Y-%m-%d_%H%MZ') 
	trop = pvds.isel(potentialVorticity=1)
	trop_pres = trop['pres']/100
	
	norm_ref, cmap_ref = ctables.registry.get_with_steps('WVCIMSS', -80., .5)
	newcmap = ListedColormap(cmap_ref(range(0, 145)))
	newnorm = Normalize(100,800)

	fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw={'projection': ccrs.NearsidePerspective(central_longitude=-76.5,central_latitude=42.5),'transform':ccrs.PlateCarree()})
	ax.add_feature(cfeature.BORDERS,edgecolor='k')
	ax.add_feature(cfeature.COASTLINE,edgecolor='k')
	ax.add_feature(cfeature.STATES,edgecolor='k')

	pvc = ax.contourf(trop_pres.longitude,trop_pres.latitude,trop_pres,transform=ccrs.PlateCarree(),levels=range(100,800,20),antialiased=True,extend='both',cmap=newcmap,norm=newnorm)
	cbar = fig.colorbar(pvc,orientation='vertical',pad=0.01,ax=ax,aspect=50,extendrect=False,ticks=range(100,800,50),shrink=0.65)
	ax.set_extent((-157,0,0,90))

	ax.set_title('2PVU Pressure',fontsize=16)
	ax.set_title('Valid: '+dtfs,loc='right',fontsize=11)
	ax.set_title('GFS Forecast via NOAA \nPlot by Jack Sillin',loc='left',fontsize=11)

	ims.append(ax)

	plt.savefig('nwh_troppres_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)

anim = manim.ArtistAnimation(animfig, ims, 20, 200,True)

anim.save('troppres_anim_test.mp4', fps=12, codec='h264', dpi=120)
