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

#rundate, runtime = spt.get_init_time('GFS')
runtime='06'
animfig = plt.figure()
ims = []

norm_ref, cmap_ref = ctables.registry.get_with_steps('WVCIMSS', -80., .5)
newcmap = ListedColormap(cmap_ref(range(0, 145)))
newnorm = Normalize(100,800)

def anim_pv():
    plt.ioff()
    f = plt.figure(figsize = (8, 4), dpi = 200)
    f.clf()

    # Create plot that can handle geographical data
    ax = f.add_subplot(111, projection = ccrs.PlateCarree())  # Add axes to figure

    # Set plot title to something meaningful
    fmt = '2PVU Pressure (hPa): '
    fmt += ' %Y-%m-%d %H:%M UTC'

    fname = '/home/jhs389/plotting/gfs.t'+runtime+'z.pgrb2.0p25.f000'
   	
    pvds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'potentialVorticity'})
    dtfs = pd.to_datetime(pvds.time.values+pvds.step.values).strftime('%Y-%m-%d_%H%MZ') 
    trop = pvds.isel(potentialVorticity=1)
    trop_pres = trop['pres']/100
	
    trop_pres.plot.contourf(ax=ax,transform=ccrs.PlateCarree(),levels=range(100,800,20),antialiased=True,extend='both',cmap=newcmap,norm=newnorm)
    
    frames = range(1,120,1)

    def anim(i):
        plt.ioff()
        ax.cla()
        ax.coastlines(resolution='50m')
        #ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))   
        
        fhr = i
        if fhr <10:
            pfhr = '00'+str(fhr)
       	elif fhr <100:
            pfhr = '0'+str(fhr)
       	else:
       		pfhr = str(fhr)
        		
        fname = '/home/jhs389/plotting/gfs.t'+runtime+'z.pgrb2.0p25.f'+pfhr
       	
       	pvds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'potentialVorticity'})
        dtfs = pd.to_datetime(pvds.time.values+pvds.step.values).strftime('%Y-%m-%d_%H%MZ') 
        trop = pvds.isel(potentialVorticity=1)
        trop_pres = trop['pres']/100
        	
        ax.set_title('2PVU Pressure',fontsize=12)
        ax.set_title('Valid: '+dtfs,loc='right',fontsize=9)
        ax.set_title('GFS Forecast via NOAA \nPlot by Jack Sillin',loc='left',fontsize=9)
        ax.set_extent((-157,0,0,90))
        
        ax.contourf(trop_pres.longitude,trop_pres.latitude,trop_pres,transform=ccrs.PlateCarree(),levels=range(100,800,20),antialiased=True,extend='both',cmap=newcmap,norm=newnorm,add_colorbar=True)
        #cbar = fig.colorbar(pvc,orientation='vertical',pad=0.01,ax=ax,aspect=50,extendrect=False,ticks=range(100,800,50),shrink=0.65)
        
        plt.ion()
        plt.draw()
        
        	#plt.savefig('nwh_troppres_'+str(fhr)+'.png',bbox_inches='tight',pad_inches=0.1)
        
        
    anim = manim.FuncAnimation(f, anim,frames, repeat=False)
        
    anim.save('dt_anim_test1.mp4', fps=12, codec='h264', dpi=120)
    plt.ion()
    
anim_pv()