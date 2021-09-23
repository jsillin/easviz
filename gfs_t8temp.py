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
from cartopy.util import add_cyclic_point 

def load_basemap(map_path):
    img=plt.imread(map_path)
    img_proj = ccrs.PlateCarree()
    img_extent = (-180, 180, -90, 90)
    return img, img_proj,img_extent
	
map_path='/home/jhs389/plotting/RenderData.tif'
bm_img,bm_proj,bm_extent=load_basemap(map_path)

import matplotlib.animation as manim

rundate, runtime = spt.get_init_time('GFS')
animfig = plt.figure(figsize=(15,15))
ims = []

# Define the colormaps
norm_ref, cmap_ref = ctables.registry.get_with_steps('rainbow', -80., .5)
newcmap = ListedColormap(cmap_ref(range(0, 175)))
newnorm = Normalize(-44,48)


def anim_t8():
    plt.ioff()
    f = plt.figure(figsize = (8, 4), dpi = 200, tight_layout=True)
    f.clf()

    # Create plot that can handle geographical data
    ax = f.add_subplot(111, projection = ccrs.PlateCarree())  # Add axes to figure

    # Set plot title to something meaningful
    fmt = '850hPa Temperature (K): '
    fmt += ' %Y-%m-%d %H:%M UTC'

    fname = '/home/jhs389/plotting/gfs.t'+runtime+'z.pgrb2.0p25.f000'
	
    
    frames = range(0,120,1)
    def anim(i):
        plt.ioff()
        ax.cla()
        ax.coastlines(resolution='50m')
        fhr = i
        if fhr <10:
            pfhr = '00'+str(fhr)
        elif fhr <100:
            pfhr = '0'+str(fhr)            
        else:
            pfhr =  str(fhr)
		
        fname = '/home/jhs389/plotting/gfs.t'+runtime+'z.pgrb2.0p25.f'+pfhr
        ds = xr.open_dataset(fname,engine='cfgrib',filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
        t8 = ds['t'].sel(isobaricInhPa=850)-273.15	
        t8 = add_cyclic_point(t8,coord=t8.lon)
        dtfs = pd.to_datetime(t8.time.values+t8.step.values).strftime('%Y-%m-%d_%H%MZ') 
        ax.set_extent((-157,0,0,90))

        ax.contourf(t8.longitude,t8.latitude,t8,norm=newnorm,cmap=newcmap,antialiased=True,levels=range(-40,40,1),alpha=0.7,transform=ccrs.PlateCarree(),add_colorbar=True)
        ax.contour(t8.longitude,t8.latitude,t8,levels=[0],colors=['blue'],transform=ccrs.PlateCarree())

        ax.set_title('850hPa Temperature (C)',fontsize=11)
        ax.set_title('Valid: '+dtfs,loc='right',fontsize=8)
        ax.set_title('GFS Forecast via NOAA \nPlot by Jack Sillin',loc='left',fontsize=8)
        
        plt.ion()
        plt.draw()

    anim = manim.FuncAnimation(f, anim,frames, repeat=False)
        
    anim.save('t8_anim_test2.mp4', fps=12, codec='h264', dpi=150)
    plt.ion()

anim_t8()