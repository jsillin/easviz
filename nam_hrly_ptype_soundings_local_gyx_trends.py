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
from metpy.plots import USCOUNTIES
from metpy.plots import SkewT
import metpy.calc as mpcalc
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import supplementary_tools as spt
# make unique directory to store output
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

def wet_bulb(temp,dewpoint):
    tdd = temp-dewpoint
    wet_bulb = temp-((1/3)*tdd)
    return wet_bulb

def fram(ice,wet_bulb,velocity):
    ilr_p = ice
    ilr_t = (-0.0071*(wet_bulb**3))-(0.039*(wet_bulb**2))-(0.3904*wet_bulb)+0.5545
    ilr_v = (0.0014*(velocity**2))+(0.0027*velocity)+0.7574

    cond_1 = np.ma.masked_where(wet_bulb>-0.35,ice)
    cond_2 = np.ma.masked_where((wet_bulb<-0.35) & (velocity>12.),ice)
    cond_3 = np.ma.masked_where((wet_bulb<-0.35) & (velocity<=12.),ice)

    cond_1 = cond_1.filled(0)
    cond_2 = cond_2.filled(0)
    cond_3 = cond_3.filled(0)

    ilr_1 = (0.7*ilr_p)+(0.29*ilr_t)+(0.01*ilr_v)
    ilr_2 = (0.73*ilr_p)+(0.01*ilr_t)+(0.26*ilr_v)
    ilr_3 = (0.79*ilr_p)+(0.2*ilr_t)+(0.01*ilr_v)

    accretion_1 = cond_1*ilr_1
    accretion_2 = cond_2*ilr_2
    accretion_3 = cond_3*ilr_3

    total_accretion=accretion_1+accretion_2+accretion_3
    return total_accretion

#grabbing data from NOMADS
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
    if int(hour) <5:
        init_hour = '00'
    elif int(hour) <10:
        init_hour = '06'
    elif int(hour) <16:
        init_hour = '12'
    elif int(hour) <21:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)
init_hour = get_init_hr(hour)
url='http://nomads.ncep.noaa.gov:80/dods/nam/nam'+mdate+'/nam_'+get_init_hr(hour)+'z'
purl='http://nomads.ncep.noaa.gov:80/dods/nam/nam'+spt.get_prev_init_time('NAM')[0]+'/nam_'+spt.get_prev_init_time('NAM')[1]+'z'
# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/NAM')


#Parse data using MetPy
ds = xr.open_dataset(url)
pds = xr.open_dataset(purl)

init_hr = dt.datetime(year,int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]
pinit_time = pds['time'][0]

lats = np.arange(12.21990800000,61.20556254545,0.110827275)
lons = np.arange((360-152.87862300000),(360.-49.47263081081),0.11338376)


#Now loop through the rest to come up with a ptype map each hour and
#ideally a four-panel accumulation map
for i in range(0,26):
    fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    pdata = pds.metpy.parse_cf()
    pdata = pdata.isel(time=i+2)

    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tcdcclm':'tcc',
        'tmpprs': 'temperature',
        'ugrd10m': 'u',
        'vgrd10m': 'v',
        'hgtprs': 'height',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf',
        'rhprs':'rh'
    })

    pdata = pdata.rename({
        'tmpprs': 'temperature',
        'tmp2m':'sfc_temp',
        'refcclm':'radar',
        'rhprs':'rh'
    })

    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

    #This extends each ptype one gridpoint outwards to prevent a gap between
    #different ptypes
    radius = 1
    kernel = np.zeros((2*radius+1,2*radius+1))
    y1,x1 = np.ogrid[-radius:radius+1,-radius:radius+1]
    mask=x1**2+y1**2 <=radius**2
    kernel[mask]=1

    snowc= gf(catsnow,np.max,footprint=kernel)
    icec = gf(catice,np.max,footprint=kernel)
    sleetc = gf(catsleet,np.max,footprint=kernel)
    rainc = gf(catrain,np.max,footprint=kernel)

    #Coordinate stuff
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    zH5_crs = data['temperature'].metpy.cartopy_crs

    #data['temperature'].metpy.convert_units('degC')
    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    pt2m = pdata['sfc_temp'].squeeze()
    pt2m = ((pt2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    wb2m = wet_bulb(t2m,td2m)
    wb2mc = (wb2m-32.)*(5./9.)

    reflectivity = data['radar'].squeeze()
    '''
    rain = np.ma.masked_where(rainc==0,new_precip)
    sleet = np.ma.masked_where(catsleet==0,new_precip)
    ice = np.ma.masked_where(catice==0,new_precip)
    snow = np.ma.masked_where(snowc==0,new_precip)
    '''
    rain = np.ma.masked_where(rainc==0,reflectivity)
    sleet = np.ma.masked_where(sleetc==0,reflectivity)
    ice = np.ma.masked_where(icec==0,reflectivity)
    snow = np.ma.masked_where(snowc==0,reflectivity)

    #Smooth rain
    #rain = ndimage.gaussian_filter(rain,sigma=1,order=0)

    mslp = data['mslp']/100.
    mslpc = mslp.squeeze()
    mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)

    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))


    ########## PLOTTING #######################################################
    #tmp_2m = ax1.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    ptmp_2m32 = ax1.contour(x,y,pt2m,colors='steelblue', alpha = 0.8, levels = [32])

    #cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
    #                    extendrect=False, ticks = range(-20,100,5))
    #cbr.set_label('2m Temperature (F)', fontsize = 14)

    #h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1080,4),linewidths=1,alpha=0.7)
    #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    #q_levs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25]
    q_levs = [0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.5]
    qarlevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]
    qaslevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5]
    qazlevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3]
    qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
    qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
    qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
    qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']
    qra_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00','#ffcc00','#ff9500','#ff4800','#ff2900','#ff1200','#ff0000','#cc0000','#990000','#990033','#b3003b','#ff3333','#ff6666','#ffffff']
    qrs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5','#4f01f6','#7a00f5','#9e00f5','#b833ff','#d280ff','#cc00f1','#ad00cc','#820099','#4700b3']
    qzr_cols = ['#ff0066','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333','#ff6666','#ff9999','#ffcccc','#ffffff']
    qip_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc','#9933ff','#bf80ff','#e6ccff','#ffffff']
    q_levs_r = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3]

    try:
        ra = ax1.contourf(x,y,rain,colors=qr_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no rain')
    try:
        sn = ax1.contourf(x,y,snow,colors=qs_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no snow')
    try:
        ip = ax1.contourf(x,y,sleet,colors=qi_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no sleet')
    try:
        zr = ax1.contourf(x,y,ice, colors=qz_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no ice')
    #ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)

    ax1.set_title('Precipitation Type and Selected Soundings',fontsize=16)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n NAM Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax1.set_extent((287.25, 291.25, 42.5, 45.75))#, crs = zH5_crs)    # Set a title and show the plot

    soundlat = 35
    soundlon = 360-99

    sound_pres = data.lev
    ptop=300
    startlon=69.2
    londelt=0.76
    sound_lons = []
    r=5
    for i in range(0,r):
        lon = -startlon-(londelt*i)
        sound_lons.append(lon)
    sound_lats = [42.7,43.2,43.7,44.25,44.75,45.3]
    for i in range(1,r):
        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.2,.15,.1))

        soundlat = sound_lats[0]
        soundlon = -(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15
        psound_temps = pdata['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        psound_rh = pdata['rh'].interp(lat=soundlat,lon=soundlon)

        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
        psound_dp = mpcalc.dewpoint_from_relative_humidity(psound_temps.data*units.degC,psound_rh.data*units.percent)
        skew.plot(sound_pres,psound_dp,'midnightblue',linewidth=2)
        skew.plot(sound_pres,sound_dp,'cornflowerblue',linewidth=2)
        skew.plot(sound_pres,psound_temps,'indianred',linewidth=2)
        skew.plot(sound_pres,sound_temps,'r',linewidth=2)

        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')
        #skew.ax.text(0.1,0.1,str(soundlat)+'  '+str(soundlon))

    for i in range(0,r):
        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.3,.15,.1))

        soundlon = -(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15
        psound_temps = pdata['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        psound_rh = pdata['rh'].interp(lat=soundlat,lon=soundlon)

        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
        psound_dp = mpcalc.dewpoint_from_relative_humidity(psound_temps.data*units.degC,psound_rh.data*units.percent)
        skew.plot(sound_pres,psound_dp,'midnightblue',linewidth=2)
        skew.plot(sound_pres,sound_dp,'cornflowerblue',linewidth=2)
        skew.plot(sound_pres,psound_temps,'indianred',linewidth=2)
        skew.plot(sound_pres,sound_temps,'r',linewidth=2)

        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')
    for i in range(0,r):
        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.4,.15,.1))

        soundlat = sound_lats[2]
        soundlon = -(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15
        psound_temps = pdata['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        psound_rh = pdata['rh'].interp(lat=soundlat,lon=soundlon)

        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
        psound_dp = mpcalc.dewpoint_from_relative_humidity(psound_temps.data*units.degC,psound_rh.data*units.percent)
        skew.plot(sound_pres,psound_dp,'midnightblue',linewidth=2)
        skew.plot(sound_pres,sound_dp,'cornflowerblue',linewidth=2)
        skew.plot(sound_pres,psound_temps,'indianred',linewidth=2)
        skew.plot(sound_pres,sound_temps,'r',linewidth=2)

        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')
    for i in range(0,r):
        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.5,.15,.1))
        soundlat = sound_lats[3]
        soundlon = -(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15
        psound_temps = pdata['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        psound_rh = pdata['rh'].interp(lat=soundlat,lon=soundlon)

        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
        psound_dp = mpcalc.dewpoint_from_relative_humidity(psound_temps.data*units.degC,psound_rh.data*units.percent)
        skew.plot(sound_pres,psound_dp,'midnightblue',linewidth=2)
        skew.plot(sound_pres,sound_dp,'cornflowerblue',linewidth=2)
        skew.plot(sound_pres,psound_temps,'indianred',linewidth=2)
        skew.plot(sound_pres,sound_temps,'r',linewidth=2)

        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')
    for i in range(0,r):
        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.6,.15,.1))
        soundlat = sound_lats[4]
        soundlon = -(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15
        psound_temps = pdata['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        psound_rh = pdata['rh'].interp(lat=soundlat,lon=soundlon)

        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
        psound_dp = mpcalc.dewpoint_from_relative_humidity(psound_temps.data*units.degC,psound_rh.data*units.percent)
        skew.plot(sound_pres,psound_dp,'midnightblue',linewidth=2)
        skew.plot(sound_pres,sound_dp,'cornflowerblue',linewidth=2)
        skew.plot(sound_pres,psound_temps,'indianred',linewidth=2)
        skew.plot(sound_pres,sound_temps,'r',linewidth=2)

        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')
    for i in range(0,r):
        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.7,.15,.1))
        soundlat = sound_lats[5]
        soundlon = -(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15
        psound_temps = pdata['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        psound_rh = pdata['rh'].interp(lat=soundlat,lon=soundlon)

        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
        psound_dp = mpcalc.dewpoint_from_relative_humidity(psound_temps.data*units.degC,psound_rh.data*units.percent)
        skew.plot(sound_pres,psound_dp,'midnightblue',linewidth=2)
        skew.plot(sound_pres,sound_dp,'cornflowerblue',linewidth=2)
        skew.plot(sound_pres,psound_temps,'indianred',linewidth=2)
        skew.plot(sound_pres,sound_temps,'r',linewidth=2)

        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')    #rows = ax1.gridlines(ylocs=sound_lats,linewidth=2, linestyle='--', edgecolor='dimgrey',draw_labels=True)
    #cols = ax1.gridlines(xlocs=sound_lons,linewidth=2, linestyle='--', edgecolor='dimgrey',draw_labels=True)

    dashed_red_line = lines.Line2D([], [], linestyle='solid', color='r', label='Latest Forecast Temperature')
    dashed_ired_line = lines.Line2D([], [],linestyle='solid',color='indianred',label='Previous Forecast Temperature')
    dashed_purple_line = lines.Line2D([],[],linestyle='dashed',color='purple',label='0C Isotherm')
    dashed_green_line = lines.Line2D([], [], linestyle='solid', color='cornflowerblue', label='Latest Forecast Dew Point')
    dashed_igreen_line = lines.Line2D([], [], linestyle='solid', color='midnightblue', label='Latest Forecast Dew Point')

    grey_line = lines.Line2D([], [], color='darkgray', label='MSLP (hPa)')
    blue_line = lines.Line2D([], [], color='b',label='2m 0C Isotherm')
    leg = ax1.legend(handles=[dashed_red_line,dashed_ired_line,dashed_green_line,dashed_igreen_line,dashed_purple_line],title='Sounding Legend',loc=4,framealpha=1)
    leg.set_zorder(500)
    plt.savefig(output_dir+'/NAM/nam_hrly_pytpe_gyx_sound_trend_v3_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.close()
    plt.clf()
    '''
    ######################################################################

    fig2 = plt.figure(figsize=(25,15))
    ax2 = fig2.add_subplot(221,projection=zH5_crs)
    ax3 = fig2.add_subplot(222,projection=zH5_crs)
    ax4 = fig2.add_subplot(223,projection=zH5_crs)
    ax5 = fig2.add_subplot(224,projection=zH5_crs)

    ax2.coastlines(resolution='10m')
    ax2.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax2.add_feature(cfeature.STATES.with_scale('10m'))

    ax3.coastlines(resolution='10m')
    ax3.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax3.add_feature(cfeature.STATES.with_scale('10m'))

    ax4.coastlines(resolution='10m')
    ax4.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax4.add_feature(cfeature.STATES.with_scale('10m'))

    ax5.coastlines(resolution='10m')
    ax5.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax5.add_feature(cfeature.STATES.with_scale('10m'))

    fig2.suptitle('GFS Accumulated Precipitation By Precipitation Type \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14)

    try:
        ra = ax5.contourf(x,y,acc_rain,colors=qra_cols,levels=qarlevs,alpha=0.7)
        rac = fig2.colorbar(ra,orientation = 'vertical', aspect = 20, ax = ax5, pad = 0.01,
                            extendrect=False, ticks = [0.5,1,1.5,2,3,4,5,7,9],shrink=0.7)
    except:
        print('no rain')
    try:
        sn = ax2.contourf(x,y,acc_snow,colors=qrs_cols,levels=qaslevs,alpha=0.7)
        snc = fig2.colorbar(sn,orientation = 'vertical', aspect = 20, ax = ax2, pad = 0.01,
                            extendrect=False, ticks = range(0,5,1),shrink=0.7)
    except:
        print('no snow')
    try:
        ip = ax3.contourf(x,y,acc_sleet,colors=qip_cols,levels=qazlevs,alpha=0.7)
        ipc = fig2.colorbar(ip,orientation = 'vertical', aspect = 20, ax = ax3, pad = 0.01,
                            extendrect=False, ticks = [0.01,0.1,0.5,1,1.5,2,3],shrink=0.7)
    except:
        print('no sleet')
    try:
        zr = ax4.contourf(x,y,acc_ice, colors=qzr_cols,levels=qazlevs,alpha=0.7)
        zrc = fig2.colorbar(zr,orientation = 'vertical', aspect = 20, ax = ax4, pad = 0.01,
                            extendrect=False, ticks = [0.01,0.1,0.5,1,1.5,2,3],shrink=0.7)
    except:
        print('no ice')

    ax2.set_title('Accumulated Liquid Equivalent Snow')
    ax3.set_title('Accumulated Liquid Equivalent Sleet')
    ax4.set_title('Accumulated Liquid Equivalent Freezing Rain')
    ax5.set_title('Accumulated Rain')
    ax2.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    fig2.tight_layout()
    plt.savefig(output_dir+'/GFS/ec_accum_ptype_'+dtfs+'.png')
    ax2.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax3.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax4.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax5.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    plt.savefig(output_dir+'/GFS/ne_accum_ptype_'+dtfs+'.png')
    plt.close()
    plt.clf()

    fig3 = plt.figure(figsize=(15,15))
    ax6 = fig3.add_subplot(111,projection=zH5_crs)
    try:
        zr_acc = ax6.contourf(x,y,acc_fram,colors=qzr_cols,levels=qazlevs,alpha=0.8)
        zrc = fig3.colorbar(zr_acc,orientation = 'horizontal', aspect = 80, ax = ax6, pad = 0.01,
                                extendrect=False, ticks = [0.01,0.1,0.5,1,1.5,2,3],shrink=0.9)
        zrc.set_label('Inches of Horizontal Accretion')
    except:
        print('no ice')

    ax6.coastlines(resolution='10m')
    ax6.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax6.add_feature(cfeature.STATES.with_scale('10m'))
    ax6.set_title('GFS Ice Accretion Forecast (FRAM)',fontsize=14)
    ax6.set_title(' \n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax6.set_title(' \n Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax6.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GFS/ec_ice_accretion_fram_'+dtfs+'.png')
    ax6.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax6.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFS/ne_ice_accretion_fram_'+dtfs+'.png')
    ax6.set_extent((260,271,36,42))#, crs = zH5_crs)    # Set a title and show the plot
    #plt.savefig(output_dir+'/GFS/kc_outer_ice_accretion_fram_'+dtfs+'.png')
    #ax6.set_extent((263,268,37.5,40.5))#, crs = zH5_crs)    # Set a title and show the plot
    #plt.savefig(output_dir+'/GFS/kc_ice_accretion_fram_'+dtfs+'.png')
    plt.close()
    plt.clf()

    print('Hour '+str(i)+' completed!')
    timeelapsed = datetime.now()-startTime
    print(timeelapsed)
    '''