export PATH=${PATH}:$HOME/gsutil

YEAR=$(date +"%Y%m%d")
echo {$YEAR}
HRRRNATURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t18z.wrfnatf*.grib2"
HRRRPRSURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t18z.wrfprsf*.grib2"
GFSURL="gs://global-forecast-system/gfs.${YEAR}/18/atmos/gfs.t18z.pgrb2.0p25.f*"

gsutil cp -R ${GFSURL} /home/jhs389/plotting/

echo "DATA RETRIEVED"
