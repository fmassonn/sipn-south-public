#!/bin/bash
#
# F. Massonnet
# Downloads sea ice concentration of NSIDC-0803 (v2) :  NASA AMSR2 polar gridded sea ice concentrations
# (replacement of SSMIS-based products, expected to cease on 31 July 2025)
# https://nsidc.org/data/user-resources/help-center/guide-nsidc-0802-and-nsidc-0803


set -o nounset
set -o errexit
set -x 

if [[ $# -ne 2 ]]
then
	echo "./retrieve_NSIDC-0803.bash YYYYMMDD YYYYMMDD"
	echo "retrieves the raw data for product NSIDC-0803 between the two dates"
	exit
fi

rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/NSIDC/NSIDC-0803/raw/

mkdir -p $outdir

currentDate=$1
hemisphere="S"

while [ "$currentDate" -le "$2" ]
do
  echo $currentDate
 
  thisYear=`date  -j -f "%Y%m%d" $currentDate "+%Y"`
  thisMonth=`date -j -f "%Y%m%d" $currentDate "+%m"`
  thisDay=`date   -j -f "%Y%m%d" $currentDate "+%d"`


  # Check file existence

  urlRoot=https://daacdata.apps.nsidc.org/pub/DATASETS
  productSpec=nsidc0803_daily_a2_seaice_conc_v2
  fileName=NSIDC-0803_SEAICE_AMSR2_${hemisphere}_${thisYear}${thisMonth}${thisDay}_v2.0.nc

  url=$urlRoot/$productSpec/$fileName

  options='--keep-session-cookies --no-check-certificate --auth-no-challenge=on --no-clobber'

  # Get status (file availability)
  status=$(curl -L -s -o /dev/null -w "%{http_code}" \
  -b ~/.urs_cookies -c ~/.urs_cookies \
  -n "$url")

  if [ "$status" = "200" ]; then
    wget $options $url -P $outdir
  else
    echo "File absent: $url"
  fi

  # Increment
  currentDate=`date -j -v+1d -f "%Y%m%d" $currentDate "+%Y%m%d"`
done

# Download grid files
productSpec=nsidc0771_polarstereo_anc_grid_info/

# LatLon
fileName=NSIDC0771_LatLon_P${hemisphere}_${hemisphere}25km_v1.1.nc
url=${urlRoot}/$productSpec/$fileName

wget $options $url -P $outdir

# CellArea
fileName=NSIDC0771_CellArea_P${hemisphere}_${hemisphere}25km_v1.1.nc
url=${urlRoot}/$productSpec/$fileName
wget $options $url -P $outdir

# Mask
productSpec=nsidc0780_seaice_masks_v1/netcdf/
fileName=NSIDC-0780_SeaIceRegions_P${hemisphere}-${hemisphere}25km_v1.0.nc
url=${urlRoot}/$productSpec/$fileName

wget $options $url -P $outdir
