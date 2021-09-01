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
import pandas as pd
import numpy as np
import xarray as xr


import cdsapi

c = cdsapi.Client()

c.retrieve(
    'cams-global-atmospheric-composition-forecasts',
    {
        'variable': [
            'black_carbon_aerosol_optical_depth_550nm', 'dust_aerosol_optical_depth_550nm', 'organic_matter_aerosol_optical_depth_550nm',
            'sea_salt_aerosol_optical_depth_550nm',
        ],
        'date': '2021-08-05/2021-08-05',
        'time': '12:00',
        'leadtime_hour': [
            '0', '1', '10',
            '100', '101', '102',
            '103', '104', '105',
            '106', '107', '108',
            '109', '11', '110',
            '111', '112', '113',
            '114', '115', '116',
            '117', '118', '119',
            '12', '120', '13',
            '14', '15', '16',
            '17', '18', '19',
            '2', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '3',
            '30', '31', '32',
            '33', '34', '35',
            '36', '37', '38',
            '39', '4', '40',
            '41', '42', '43',
            '44', '45', '46',
            '47', '48', '49',
            '5', '50', '51',
            '52', '53', '54',
            '55', '56', '57',
            '58', '59', '6',
            '60', '61', '62',
            '63', '64', '65',
            '66', '67', '68',
            '69', '7', '70',
            '71', '72', '73',
            '74', '75', '76',
            '77', '78', '79',
            '8', '80', '81',
            '82', '83', '84',
            '85', '86', '87',
            '88', '89', '9',
            '90', '91', '92',
            '93', '94', '95',
            '96', '97', '98',
            '99',
        ],
        'type': 'forecast',
        'format': 'grib',
    },
    'cams_aerosol.grib')
	
ds = xr.open_dataset('cams_aerosol.grib',engine='cfgrib')

orgm = ds['omaod550']
salt = ds['ssaod550']
dust = ds['duaod550']
bcar = ds['bcaod550']

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


map_path='RenderData.tif'
bm_img,bm_proj,bm_extent=load_basemap(map_path)
print(np.shape(bm_img))

#alphas = Normalize(0, .3, clip=True)(bcar.isel(step=1))
#alphas = np.clip(alphas, .4, 1)  # alpha value clipped at the bottom at .4

for i in range(len(bcar.step)):
    fig, ax = plt.subplots(1, 1, figsize=(15, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.imshow(bm_img,extent=bm_extent,transform=bm_proj,origin='upper')
    dtfs = pd.to_datetime(bcar.isel(step=i).step.values+bcar.time.values).strftime('%Y-%m-%d_%H%MZ') 

    #ax.add_feature(cfeature.BORDERS,edgecolor='lemonchiffon')
    #ax.add_feature(cfeature.COASTLINE,edgecolor='lemonchiffon')
    #ax.add_feature(cfeature.STATES,edgecolor='lemonchiffon')
    ax.contourf(orgm.longitude,orgm.latitude,(orgm.isel(step=i)+bcar.isel(step=i)),levels=np.linspace(0.02,1.5,125),extend='max',antialiased=True,cmap=cmaporgm,transform=ccrs.PlateCarree())
    ax.contourf(salt.longitude,salt.latitude,salt.isel(step=i),levels=np.linspace(0.02,.75,125),extend='max',antialiased=True,cmap=cmapsalt,transform=ccrs.PlateCarree())
    ax.contourf(dust.longitude,dust.latitude,dust.isel(step=i),levels=np.linspace(0.02,1.5,125),extend='max',antialiased=True,cmap=cmapdust,transform=ccrs.PlateCarree())
    #ax.contourf(bcar.longitude,bcar.latitude,bcar.isel(step=i),levels=np.linspace(0.02,1,75),extend='max',antialiased=True,cmap=cmapsmoke,transform=ccrs.PlateCarree())

    ax.set_extent((230,350,-3,50))
    ax.set_title('Total Column Aerosol Optical Depth',fontsize=14)
    ax.set_title('Valid: '+dtfs,loc='right',fontsize=11)
    ax.set_title('CAMS Forecast',loc='left',fontsize=11)
    ax.text(-129,-4,"Plot by Jack Sillin, data via Copernicus",fontsize=11,color='white')
    fig.set_facecolor('white')

    blue = mpatches.Patch(color='cyan',label='Sea Salt')
    red = mpatches.Patch(color='yellow',label='Smoke')
    orange = mpatches.Patch(color='orange',label='Dust')
    ax.legend(handles=[red,orange,blue],title='Aerosol Legend',loc=4,framealpha=1)
    plt.savefig('US_aerosol_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((230,360,-3,60))
    plt.savefig('natl_cams_aerosol_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
