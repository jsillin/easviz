import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import geopandas as gpd

per = gpd.read_file('https://opendata.arcgis.com/datasets/2191f997056547bd9dc530ab9866ab61_0.geojson')

runtime = '18'
for i in range(0,49):
    if i<10:
        pfhour = '0'+str(i)
    else:
        pfhour = str(i)
    ds = xr.open_dataset('hrrr.t'+runtime+'z.wrfprsf'+pfhour+'.grib2',engine='cfgrib',filter_by_keys={'typeOfLevel':'atmosphereSingleLayer','stepType':'instant'})
    mden = ds['unknown']
    print(mden)
    dtfs = pd.to_datetime(mden.time.values).strftime('%Y-%m-%d_%H%MZ')
    dtfs2 = pd.to_datetime((mden.step.values+mden.time.values)).strftime('%Y-%m-%d_%H%MZ')
    print(dtfs2)
    fig = plt.figure(figsize=(12,12))
    prj = ccrs.PlateCarree()
    ax = fig.add_subplot(111,projection=prj)
    ax.add_feature(cfeature.BORDERS,edgecolor='lemonchiffon')
    ax.add_feature(cfeature.COASTLINE,edgecolor='lemonchiffon')
    ax.add_feature(cfeature.STATES,edgecolor='lemonchiffon')
    
    mdc = ax.contourf(mden.longitude,mden.latitude,mden*1000,levels=range(0,750,10),extend='max',cmap='bone',transform=prj)
    cbr = fig.colorbar(mdc,orientation='vertical',ax=ax,extendrect=False,pad=0.01,ticks=range(0,750,50),shrink=0.60)
    cbr.set_label('mg/m^2')
    ax.set_title('Total Atmosphere Smoke',fontsize=14)
    ax.set_title(dtfs2,loc='right',fontsize=12)
    ax.set_title('HRRR Forecast init '+dtfs,fontsize=12)

    outfname = 'HRRR_total_smoke_'+str(i)
    ax.set_extent((-130,-100,30,49))
    plt.savefig('HRRR_total_smoke_USWest_f'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-125,-110,40,49.5))
    plt.savefig('HRRR_total_smoke_USNW_f'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
