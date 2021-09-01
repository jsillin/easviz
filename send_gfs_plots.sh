export PATH=${PATH}:$HOME/gsutil

gsutil cp /home/jhs389/plotting/*hurr*.png gs://nwp-plots/gfsdash/ &
gsutil cp /home/jhs389/plotting/*t8*.png gs://nwp-plots/gfsmisc/ &
gsutil cp /home/jhs389/plotting/*pv*.png gs://nwp-plots/gfsdash/ &
gsutil cp /home/jhs389/plotting/*rh*.png gs://nwp-plots/gfsdash/ &
gsutil cp /home/jhs389/plotting/*ov*.png gs://nwp-plots/gfsdash/ &
gsutil cp /home/jhs389/plotting/*sst*.png gs://nwp-plots/gfsdash/ &
