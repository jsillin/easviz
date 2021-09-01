export PATH=${PATH}:$HOME/gsutil

YEAR=$(date +"%Y%m%d")
HRRRNATURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t06z.wrfnatf*.grib2"
HRRRPRSURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t06z.wrfprsf*.grib2"
GFSURL="gs://global-forecast-system/gfs.${YEAR}/06/atmos/gfs.t06z.pgrb2.0p25.f*"

gsutil cp -R ${HRRRNATURL} /home/jhs389/plotting/
gsutil cp -R ${HRRRPRSURL} /home/jhs389/plotting/
gsutil cp -R ${GFSURL} /home/jhs389/plotting/

echo ${URL}
