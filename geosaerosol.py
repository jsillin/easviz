import xarray as xr
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

ds = xr.open_dataset('https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/inst1_2d_hwl_Nx/inst1_2d_hwl_Nx.20210808_12')

orgm = ds['ocexttau']
salt = ds['ssexttau']
dust = ds['duexttau']
bcar = ds['bcexttau']

def load_basemap(map_path):
    img=plt.imread(map_path)
    img_proj = ccrs.PlateCarree()
    img_extent = (-180, 180, -90, 90)
    return img, img_proj,img_extent
	
# Define the colormaps
colors = [(1,0.5,0,c) for c in np.linspace(0.1,.9,80)]
cmapdust = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=8)
colors = [(1,0,0,c) for c in np.linspace(0.1,.9,80)]
cmapsmoke = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=8)
colors = [(0.4,1,1,c) for c in np.linspace(0,.9,80)]
cmapsalt = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=8)
colors = [(1,1,0.4,c) for c in np.linspace(0,.9,80)]
cmaporgm = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=8)

for i in range(1,123):
    map_path='RenderData.tif'
    bm_img,bm_proj,bm_extent=load_basemap(map_path)
    fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw={'projection': ccrs.Robinson(),'transform':ccrs.PlateCarree()})
    ax.imshow(bm_img,extent=bm_extent,transform=bm_proj,origin='upper')
    dtfs = pd.to_datetime(bcar.isel(time=i).time.values).strftime('%Y-%m-%d_%H%MZ') 
    print(dtfs)
    #ax.add_feature(cfeature.BORDERS,edgecolor='lemonchiffon')
    #ax.add_feature(cfeature.COASTLINE,edgecolor='lemonchiffon')
    #ax.add_feature(cfeature.STATES,edgecolor='lemonchiffon')
    ax.contourf(orgm.lon,orgm.lat,(orgm.isel(time=i)+bcar.isel(time=i)),levels=np.linspace(0.02,1.5,125),extend='max',antialiased=True,cmap=cmaporgm,transform=ccrs.PlateCarree())
    ax.contourf(salt.lon,salt.lat,salt.isel(time=i),levels=np.linspace(0.02,.5,125),extend='max',antialiased=True,cmap=cmapsalt,transform=ccrs.PlateCarree())
    ax.contourf(dust.lon,dust.lat,dust.isel(time=i),levels=np.linspace(0.05,1.25,125),extend='max',antialiased=True,cmap=cmapdust,transform=ccrs.PlateCarree())
    #ax.contourf(bcar.longitude,bcar.latitude,bcar.isel(step=i),levels=np.linspace(0.02,1,75),extend='max',antialiased=True,cmap=cmapsmoke,transform=ccrs.PlateCarree())

    #ax.set_extent((230,350,-3,50))
    ax.set_title('Total Column Aerosol Optical Depth',fontsize=16)
    ax.set_title('Valid: '+dtfs,loc='right',fontsize=11)
    ax.set_title('CAMS Forecast via Copernicus \nPlot by Jack Sillin',loc='left',fontsize=11)
    fig.set_facecolor('white')

    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.005
    bottom = axes_bbox.y0 + 0.001

    #ax.text(left,bottom,"Plot by Jack Sillin, data via Copernicus",fontsize=11,color='white')
    ax.set_extent((-182,182,-90,90))
    blue = mpatches.Patch(color='cyan',label='Sea Salt')
    red = mpatches.Patch(color='yellow',label='Smoke')
    orange = mpatches.Patch(color='orange',label='Dust')
    ax.legend(handles=[red,orange,blue],title='Aerosol Legend',loc=4,framealpha=1)
    #plt.savefig('geos_aerosol_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    #ax.set_extent((230,360,-3,60))
    plt.savefig('global_geos_aerosol_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((230,360,-3,60))
    plt.savefig('natl_geos_aerosol_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
