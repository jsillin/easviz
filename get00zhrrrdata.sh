export PATH=${PATH}:$HOME/gsutil

YEAR=$(date +"%Y%m%d")
HRRRNATURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t00z.wrfnatf*.grib2"
HRRRPRSURL="gs://high-resolution-rapid-refresh/hrrr.${YEAR}/conus/hrrr.t00z.wrfprsf*.grib2"

gsutil cp -R ${HRRRNATURL} /home/jhs389/plotting/
gsutil cp -R ${HRRRPRSURL} /home/jhs389/plotting/

echo ${URL}
