export PATH=${PATH}:$HOME/gsutil

YEAR=$(date +"%Y%m%d")
HRRRNATURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t00z.wrfnatf*.grib2"
HRRRPRSURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t00z.wrfprsf*.grib2"
GFSURL="gs://global-forecast-system/gfs.${YEAR}/00/atmos/gfs.t00z.pgrb2.0p25.f*"

gsutil cp -R ${GFSURL} /home/jhs389/plotting/

echo ${URL}
