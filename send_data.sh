export PATH=${PATH}:$HOME/gsutil

gsutil cp *ts*.png gs://nwp-plots/psmoke/
gsutil cp *ns*.png gs://nwp-plots/psmoke/
gsutil cp *geos*.png gs://nwp-plots/camsaerosol/
gsutil cp *panel*.png gs://nwp-plots/hrrrdash/
gsutil cp *hurr*.png gs://nwp-plots/gfsdash/
gsutil cp *t8*.png gs://nwp-plots/gfsmisc/
gsutil cp *pv*.png gs://nwp-plots/gfsmisc/
