export PATH=${PATH}:$HOME/gsutil

gsutil cp /home/jhs389/*ts*.png gs://nwp-plots/psmoke/
gsutil cp /home/jhs389/*ns*.png gs://nwp-plots/psmoke/
gsutil cp /home/jhs389/*panel*.png gs://nwp-plots/hrrrdash/

