#!/bin/bash
#
# Downloads and updates the 
# Sea ice Concentration product for OSISAF
# SSMIS Sea Ice Concentration Maps on 10 km Polar Stereographic Grid
# OSI-401-b
# http://osisaf.met.no/p/ice/index.html#conc-ssmis
# http://osisaf.met.no/docs/osisaf_cdop2_ss2_pum_ice-conc_v1p4.pdf
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
        echo "./retrieve_OSI-401-b.bash YYYYMMDD YYYYMMDD"
        echo "retrieves the raw OSI-401-b data between the two dates"
        exit
fi

rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/OSI-SAF/OSI-401-b/raw/

mkdir -p $outdir

currentDate=$1

while [ "$currentDate" -le "$2" ]
do
  echo $currentDate


  thisYear=`date  -j -f "%Y%m%d" $currentDate "+%Y"`
  thisMonth=`date -j -f "%Y%m%d" $currentDate "+%m"`
  thisDay=`date   -j -f "%Y%m%d" $currentDate "+%d"`

  rootaddress="ftp://osisaf.met.no/archive/ice/conc/"
  url=$rootaddress/$thisYear/$thisMonth/ice_conc_sh_polstere-100_multi_${thisYear}${thisMonth}${thisDay}1200.nc

  wget -N -c $url -P $outdir 

  # Increment
  currentDate=`date -j -v+1d -f "%Y%m%d" $currentDate "+%Y%m%d"`
done
