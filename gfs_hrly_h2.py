import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import xarray as xr
import metpy
from datetime import datetime
import datetime as dt
from metpy.units import units
import scipy.ndimage as ndimage
from metpy.plots import USCOUNTIES
import cartopy
from scipy.ndimage.filters import generic_filter as gf


def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

startTime=datetime.now()

m_date='20200903'
m_hour='12'

year = startTime.year

if startTime.month <10:
    month = '0'+str(startTime.month)
else:
    month = str(startTime.month)

if startTime.day <10:
    day = '0'+str(startTime.day)
else:
    day = str(startTime.day)

if startTime.hour <10:
    hour = '0'+str(startTime.hour)
else:
    hour = str(startTime.hour)

mdate = str(year)+str(month)+str(day)

def get_init_hr(hour):
    if int(hour) <6:
        init_hour = '00'
    elif int(hour) <11:
        init_hour = '06'
    elif int(hour) <17:
        init_hour = '12'
    elif int(hour) <22:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs'+mdate+'/gfs_0p25_1hr_'+get_init_hr(hour)+'z'
init_hour = get_init_hr(hour)
'''
for i in range(119):
    fhr = i+1
'''
# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFS')
#Parse data using MetPy
ds = xr.open_dataset(url)
init_hr = dt.datetime(int(year),int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(15,70,0.25)
lons = np.arange(220,310,0.25)

for i in range(1,120):
    fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)
    #Rename variables to useful things
    data = data.rename({
        'absvprs':'avort',
        'hgtprs':'gph',
        'rhprs':'rh',
        'tmpprs':'temp',
        'ugrdprs':'u',
        'vgrdprs': 'v',
    })

    vertical, = data['temp'].metpy.coordinates('vertical')
    time = data['temp'].metpy.time
    zH5_crs = data['temp'].metpy.cartopy_crs

    t2 = data['temp'].sel(lev=250.0,lat=lats,lon=lons)-273.15
    u2 = data['u'].sel(lev=250.0,lat=lats,lon=lons).squeeze()*1.94384449
    v2 = data['v'].sel(lev=250.0,lat=lats,lon=lons).squeeze()*1.94384449
    av2 = data['avort'].sel(lev=250.0,lat=lats,lon=lons).squeeze()*1e5
    rh2 = data['rh'].sel(lev=250.0,lat=lats,lon=lons).squeeze()
    h2 = data['gph'].sel(lev=250.0,lat=lats,lon=lons).squeeze()
    x, y = t2.metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    wind_slice = slice(5,-5,5)

    wspd = ((u2**2)+(v2**2))**.5

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))

    #fig.suptitle("NAM Forecast valid at " + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=36)

    ########## PLOTTING #######################################################
    h2c = ax1.contour(x,y,h2,colors='dimgray', levels = range(9000,11000,60),linewidths=1.5)
    #t2c = ax1.contour(x,y,t2,colors='r', levels = range(-60,0,5),linestyles='dashed',linewidths=1)
    w2c = ax1.contourf(x,y,wspd,cmap='cool',levels=range(50,200,5),alpha=0.8)
    a5cb = fig.colorbar(w2c, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(50,200,20))
    a5cb.set_label('250mb Wind Speed (kts)', fontsize = 12)
    ax1.barbs(x[wind_slice],y[wind_slice],u2[wind_slice,wind_slice],v2[wind_slice,wind_slice], length=7)

    #h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)
    ax1.set_title('250mb Heights (m) and Winds (kts)',fontsize=16)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax1.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GFS/gfs_hrly_h2jet_v1_'+str(i)+'.png')

    ########## PLOT 2 #######################################################
    wind_slice_s = slice (10,-10,10)
    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111,projection=zH5_crs)
    ax2.coastlines(resolution='50m')
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax2.add_feature(cfeature.STATES.with_scale('50m'))
    w2c2 = ax2.contourf(x,y,wspd,cmap='cool',levels=range(50,200,5),alpha=0.8)
    h2c2 = ax2.contour(x,y,h2,colors='dimgray', levels = range(9000,11000,60),linewidths=1.5)
    w2cb2 = fig2.colorbar(w2c2, orientation = 'horizontal', aspect = 80, ax = ax2, pad = 0.01,
                        extendrect=False, ticks = range(50,200,20))
    w2cb2.set_label('250mb Wind Speed (kts)', fontsize = 12)
    ax2.barbs(x[wind_slice_s],y[wind_slice_s],u2[wind_slice_s,wind_slice_s],v2[wind_slice_s,wind_slice_s], length=7)

    #h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)
    ax2.set_title('250mb Heights (m) and Winds (kts)',fontsize=16)
    ax2.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax2.set_title('\n GFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax2.set_extent((225, 300, 20, 65))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GFS/gfs_hrly_h2jetCONUS_v1_'+str(i)+'.png')

    fcst_hr = str(0)
    print('Hour '+str(i)+' completed!')
    plt.close()
    timeelapsed = datetime.now()-startTime
    print(timeelapsed)







'''
url= 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs20200903/gfs_0p25_1hr_12z'
ds = xr.open_dataset(url)
t2m_ds = ds['tmp2m']
init_hr = t2m_ds['time'][0].values
#fc_hr = t2m.ds['time'][i].values
lats = np.arange(20,50,0.25)
lons = np.arange(240,300,0.25)
t2m = t2m_ds.sel(time = init_hr, lat = lats, lon = lons)
print(t2m)

fig = plt.figure(figsize = (12,12))
fig.clf()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_extent((240,300, 20, 50), crs = ccrs.PlateCarree())
t2m_c = ax.contourf(t2m, cmap='RdPu')
plt.savefig('testingnomads6.png')
'''
