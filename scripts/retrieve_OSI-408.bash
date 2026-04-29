#!/bin/bash
#
# Downloads and updates the 
# Sea ice Concentration product for OSISAF
# AMSR2 Sea Ice Concentration Maps on 10 km Polar Stereographic Grid
# OSI-408
# https://osi-saf.eumetsat.int/products/osi-408
#
# www.osisaf.met.no
#
# Contact: François Massonnet - francois.massonnet@uclouvain.be
#

set -o nounset
set -o errexit
#set -x

if [[ $# -ne 2 ]]
then
        echo "./retrieve_OSI-408.bash YYYYMMDD YYYYMMDD"
        echo "retrieves the raw OSI-408 data between the two dates"
        exit
fi

rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/OSI-SAF/OSI-408/raw/

mkdir -p $outdir

currentDate=$1

while [ "$currentDate" -le "$2" ]
do
  echo $currentDate


  thisYear=`date  -j -f "%Y%m%d" $currentDate "+%Y"`
  thisMonth=`date -j -f "%Y%m%d" $currentDate "+%m"`
  thisDay=`date   -j -f "%Y%m%d" $currentDate "+%d"`

  rootaddress="ftp://osisaf.met.no/archive/ice/conc_amsr/"
  url=$rootaddress/$thisYear/$thisMonth/ice_conc_sh_polstere-100_amsr2_${thisYear}${thisMonth}${thisDay}1200.nc

  if curl -s --fail  $url > /dev/null ; then # have to use curl for ftp servers, otherwise it downloads the file
    echo "File found" 
    wget -N -c $url -P $outdir 
  else
    echo "File not found , skipping: $url"
  fi

  # Increment
  currentDate=`date -j -v+1d -f "%Y%m%d" $currentDate "+%Y%m%d"`
done


# Retrieve grid cell area (Thomas Lavergne's  given link on 28 April 2026)

urllink="ftp://osisaf.met.no/docs/tools/lmask_sh_stere_100.nc.gz"

wget -N -c $urllink -N -P $outdir

gunzip -f $outdir/lmask_sh_stere_100.nc.gz > $outdir/lmask_sh_stere_100.nc
