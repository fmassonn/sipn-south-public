#!/bin/bash
#
# F. Massonnet
# Downloads sea ice concentration of NASA Team near real time (NSIDC 0081)

set -o nounset
set -o errexit
#set -x 

if [[ $# -ne 2 ]]
then
	echo "./retrieve_NSIDC-0081.bash YYYYMMDD YYYYMMDD"
	echo "retrieves the raw NSIDC data between the two dates"
	exit
fi

rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/NSIDC/NSIDC-0081/raw/

mkdir -p $outdir

currentDate=$1

threDate=`date -j -f "%Y%m%d" 20160401 + "%Y%m%d"` # Threshold date corresponding to sensor change


while [ "$currentDate" != $2 ]
do
  echo $currentDate

  if [[ $currentDate -ge $threDate ]]
  then
    sensor="f18"
  else
    sensor="f17"
  fi
 
  thisYear=`date  -j -f "%Y%m%d" $currentDate "+%Y"`
  thisMonth=`date -j -f "%Y%m%d" $currentDate "+%m"`
  thisDay=`date   -j -f "%Y%m%d" $currentDate "+%d"`

  url=https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0081.002/${thisYear}.${thisMonth}.${thisDay}/NSIDC0081_SEAICE_PS_S25km_${currentDate}_v2.0.nc

  wget -ncd --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off $url -P $outdir

  # Increment
  currentDate=`date -j -v+1d -f "%Y%m%d" $currentDate "+%Y%m%d"`
done

# Download grid files

wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25lats_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25lons_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25area_v3.dat -P ${outdir}


wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25lats_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25lons_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25area_v3.dat -P ${outdir}
