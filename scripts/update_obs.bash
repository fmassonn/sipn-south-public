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
  #python3=/usr/local/bin/python3
  #python3=/usr/bin/python3
  python3=/opt/anaconda3/bin/python3
else
  python3=python3
fi

# Dates defining the name of the SIPN South experiment and the period to analyze
# It is very important to be consistent here (i.e., to ensure that the three variables target, dateStart, dateEnd, correspond to the same period).
# These strings will be used to define the text files containing the data, and will be placed in the folder with the same name.

target="2025-2026"
dateStart=20251201
dateEnd=20260228 # This can be a date later than today. The download will just not happen if that data does not exist.


# Step 1 : Retrieving the raw data
./retrieve_NSIDC-0803.bash $dateStart $dateEnd
./retrieve_OSI-401-b.bash  $dateStart $dateEnd


# Step 2: Formatting to ELIC compliant format

$python3 ./format_NSIDC-0803.py $dateStart $dateEnd
$python3 ./format_OSI-401-b.py  $dateStart $dateEnd

mkdir -p ../data/${target}/netcdf

cp $TECLIM_CLIMATE_DATA/obs/ice/siconc/NSIDC/NSIDC-0803/processed/native/siconc_SIday_NSIDC-0803_r1i1p1_${dateStart}-${dateEnd}_sh.nc ../data/${target}/netcdf/NSIDC-0803_000_${dateStart}-${dateEnd}_concentration.nc

cp $TECLIM_CLIMATE_DATA/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_${dateStart}-${dateEnd}_sh.nc ../data/${target}/netcdf/OSI-401-b_000_${dateStart}-${dateEnd}_concentration.nc

$python3 obs2CSV.py $target $dateStart $dateEnd

