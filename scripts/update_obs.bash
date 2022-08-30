#!/bin/bash
#
# F. Massonnet

# Update obs verif data

# !!!!
#
# MAKE SURE to update the first year in the *.R and *.py scripts below
#
#    AND 
#
# TO UPDATE THE YEAR IN obs2CSV.py as well
#
# !!!

set -o nounset
set -o errexit
#set -x 

if [[ $HOSTNAME = pelican ]]
then
  source ~/module_load.txt
fi

./retrieve_NSIDC-0081.bash
./retrieve_OSI-401-b.bash
Rscript format_NSIDC-0081.R
python3 ./format_OSI-401-b.py

# Year of 1 Dec of initialization
yearb=2021
yearbp1=$(( $yearb + 1 ))

#target=$yearb-$yearbp1
target="TOP2022"

isleap() { date -d $1-02-29 &>/dev/null && true  || false ; }


date -d $yearb-02-29 &>/dev/null && leap="T" || leap="F"

echo $leap

if [[ $leap == "T" ]]
then
  if [[ $target == "TOP2022" ]]
  then
    t1=487  # Subsetting May 1 - Aug 31 of the second year
    t2=609
  else
    t1=336
    t2=425
  fi
elif [[ $leap == "F" ]]
then
  if [[ $target == "TOP2022" ]]
  then
    t1=486
    t2=608
  else
    t1=335
    t2=424
  fi
fi

echo $t1
echo $t2

mkdir -p ../data/${yearb}-${yearbp1}/netcdf

ncks -F -O -d time,$t1,$t2 $TECLIM_CLIMATE_DATA/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/siconc_SIday_NSIDC-0081_r1i1p1_${yearb}0101-${yearbp1}1231_sh.nc ../data/${target}/netcdf/NSIDC-0081_000_concentration.nc

ncks -F -O -d time,$t1,$t2 $TECLIM_CLIMATE_DATA/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_${yearb}0101-${yearbp1}1231_sh.nc ../data/${target}/netcdf/OSI-401-b_000_concentration.nc


python3 obs2CSV.py

