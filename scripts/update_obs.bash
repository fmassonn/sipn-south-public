#!/bin/bash
#
# F. Massonnet
# Update obs verif data

set -o nounset
set -o errexit
set -x 

source ~/module_load.txt
#. ./retrieve_NSIDC-0081.bash
#. ./retrieve_OSI-401b.bash
#Rscript format_NSIDC-0081.R
#python ./format_OSI-401b.py

# Year of 1 Dec of initialization
yearb=2020
yearbp1=$(( $yearb + 1 ))
isleap() { date -d $1-02-29 &>/dev/null && true  || false ; }

date -d $yearb-02-29 &>/dev/null && leap="T" || leap="F"

echo $leap

if [[ $leap == "T" ]]
then
  t1=336
  t2=425
elif [[ $leap == "F" ]]
then
  t1=335
  t2=424
fi

echo $t1
echo $t2

ncks -F -O -d time,$t1,$t2 ../data/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/siconc_SIday_NSIDC-0081_r1i1p1_${yearb}0101-${yearbp1}1231_sh.nc ../data/${yearb}-${yearbp1}/netcdf/NSIDC-0081_000_concentration.nc

ncks -F -O -d time,$t1,$t2 ../data/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_${yearb}0101-${yearbp1}1231_sh.nc ../data/${yearb}-${yearbp1}/netcdf/OSI-401-b_000_concentration.nc

python3 obs2CSV.py

