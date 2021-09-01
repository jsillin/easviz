import wget
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

for i in range(1,73):
    if i <10:
        phr = '0'+str(i)
    else:
        phr = str(i)
    url='https://ftpprd.ncep.noaa.gov/data/nccf/com/aqm/prod/aqm.20210803/aqm.t06z.pm25.f'+phr+'.148.grib2'
    print(url)
    wget.download(url)
    ds = xr.open_dataset('aqm.t06z.pm25.f'+phr+'.148.grib2',engine='cfgrib')
    pm2 = ds['pmtf']


    fig = plt.figure(figsize=(15,12))
    prj = ccrs.PlateCarree()
    ax = fig.add_subplot(111,projection=prj)
    ax.add_feature(cfeature.BORDERS,edgecolor='lightgray')
    ax.add_feature(cfeature.COASTLINE,edgecolor='lightgray')
    ax.add_feature(cfeature.STATES,edgecolor='lightgray')

    dtfs = pd.to_datetime(pm2.step.values+pm2.time.values).strftime('%Y-%m-%d_%H%MZ')
    pmc = ax.contourf(pm2.longitude,pm2.latitude,pm2,levels=[0,12,35.5,55.5,150.5,250.5,350.5],colors=['green','yellow','orange','red','purple','maroon'])
    cbr = fig.colorbar(pmc,orientation='horizontal',ax=ax,extendrect=False,pad=0.01,ticks=[0,12,35.5,55.5,150.5,250.5,350.5],shrink=0.85)
    cbr.set_label('10^-6 g/m^3')
    ax.set_title('PM2.5 Concentration Shaded by AQI Category')
    ax.set_title(dtfs,loc='right')
    ax.set_title('NOAA AQM',loc='left')
    ax.set_extent((-130,-100,30,49))
    plt.savefig('WC_aqi_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-125,-110,40,49.5))
    plt.savefig('NW_aqi_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-100,-65,30,49))
    plt.savefig('EC_aqi_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((-130,-65,24,50))
    plt.savefig('US_smoke_perim_'+str(i)+'.png',bbox_inches='tight',pad_inches=0.1)

