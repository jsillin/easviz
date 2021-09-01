import xarray as xr
from astral.sun import sun
import datetime
from astral import LocationInfo
from timezonefinder import TimezoneFinder
import pandas as pd
from geopy.geocoders import Nominatim
from geopy.extra.rate_limiter import RateLimiter
from geopy.extra.rate_limiter import RateLimiter

def get_wxdata(lat,lon,date):
    ds = xr.open_dataset('http://nomads.ncep.noaa.gov:80/dods/blend/blend'+date+'/blend_1hr_12z')
    tmax = ((ds['tmax2m']-273.15)*1.8)+32
    tmin = ((ds['tmp2m']-273.15)*1.8)+32
    pop = ds['apcp254gtsfc']
    wind = ds['wind10m']*2.23694
    pot = ds['tstmsfc']
    
    varlist = [tmax,tmin,pop,wind,pot]
    varnames = ['tmax','tmin','pop','wind','pot']
    pointvalues = []

    for i in range(len(varlist)):
        var = varlist[i]
        varname = varnames[i]
        if varname == 'tmin':
            pointvar = var.interp(lat=lat,lon=lon)
            pointvalues.append(round(float(pointvar.min(dim='time').values),0))        
        else:
            pointvar = var.interp(lat=lat,lon=lon)
            pointvalues.append(round(float(pointvar.max(dim='time').values),0))

    return pointvalues

def get_sunset(lat,lon,date):
    tf = TimezoneFinder()
    timezone = tf.timezone_at(lng=lon,lat=lat)
    loc = LocationInfo('loc','region',timezone,lat,lon)
    s = sun(loc.observer,date=datetime.date(int(date[0:4]),int(date[4:6]),int(date[6:8])),tzinfo=loc.timezone)
    sunset = s['sunset'].strftime('%I:%M')
    return sunset
	
def get_venuedata(vdate):
    df = pd.read_csv('alanis_full.csv')
    try:
        venue = df.loc[df['Date']==vdate]
        lat = float(venue['Lat'])
        lon = float(venue['Lon'])
    except:
        lat='No Concert Today'
        lon = 'No Concert Today'
    return lat,lon

current_et = datetime.datetime.utcnow()-datetime.timedelta(hours=5)	
date = current_et.strftime('%Y%m%d')
vdate = current_et.strftime('%-m/%-d')

lat,lon=get_venuedata(vdate)
wx_data = get_wxdata(lat,lon,date)
wx_data.insert(2,get_sunset(30.135799,-97.646202,'20210809'))
data_df = pd.DataFrame(wx_data)
data_df.to_csv('tourforecast.csv')

