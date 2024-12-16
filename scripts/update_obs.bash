#!/bin/bash
#
# F. Massonnet

# Update observational verification data

set -o nounset
set -o errexit  # This is commented out because otherwise, the script stop when not finding obs after the latest date (when we are in the middle of the forecasting season)
set -x 

if [[ $HOSTNAME = pelican ]]
then
  source ~/module_load.txt
elif [[ $HOSTNAME = mac-SE24-401.local ]]
then
  python3=/usr/local/bin/python3
else
  python3=python3
fi

# Dates defining the name of the SIPN South experiment and the period to analyze
# Very important to be consistent here!!!
target="2024-2025"
dateStart=20241201
dateEnd=20241213

# Step 1 : Retrieving the raw data

./retrieve_NSIDC-0081.bash $dateStart $dateEnd
./retrieve_OSI-401-b.bash  $dateStart $dateEnd

# Step 2: Formatting to TECLIM compliant format

$python3 ./format_NSIDC-0081.py $dateStart $dateEnd
$python3 ./format_OSI-401-b.py  $dateStart $dateEnd

mkdir -p ../data/${target}/netcdf

cp $TECLIM_CLIMATE_DATA/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/siconc_SIday_NSIDC-0081_r1i1p1_${dateStart}-${dateEnd}_sh.nc ../data/${target}/netcdf/NSIDC-0081_000_${dateStart}-${dateEnd}_concentration.nc

cp $TECLIM_CLIMATE_DATA/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_${dateStart}-${dateEnd}_sh.nc ../data/${target}/netcdf/OSI-401-b_000_${dateStart}-${dateEnd}_concentration.nc

$python3 obs2CSV.py $target $dateStart $dateEnd

