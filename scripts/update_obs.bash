#!/bin/bash
#
# F. Massonnet

# Update observational verification data

set -o nounset
set -o errexit
set -x 

if [[ $HOSTNAME = pelican ]]
then
  source ~/module_load.txt
fi

# Dates defining the name of the SIPN South experiment and the period to analyze
# Very important to be consistent here!!!
target="TOP2022"
dateStart=20220501
dateEnd=20220831

# Step 1 : Retrieving the raw data

./retrieve_NSIDC-0081.bash $dateStart $dateEnd
./retrieve_OSI-401-b.bash  $dateStart $dateEnd

# Step 2: Formatting to TECLIM compliant format

python3 ./format_NSIDC-0081.py $dateStart $dateEnd
python3 ./format_OSI-401-b.py  $dateStart $dateEnd

mkdir -p ../data/${target}/netcdf

cp $TECLIM_CLIMATE_DATA/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/siconc_SIday_NSIDC-0081_r1i1p1_${dateStart}-${dateEnd}_sh.nc ../data/${target}/netcdf/NSIDC-0081_000_concentration.nc

cp $TECLIM_CLIMATE_DATA/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_${dateStart}-${dateEnd}_sh.nc ../data/${target}/netcdf/OSI-401-b_000_concentration.nc


exit
python3 obs2CSV.py

