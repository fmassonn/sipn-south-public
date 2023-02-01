#!/bin/bash
#
# F. Massonnet
# Downloads sea ice concentration of NASA Team near real time (NSIDC 0081)
set -o nounset
#set -o errexit
set -x 
rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/NSIDC/NSIDC-0081/raw/

mkdir -p $outdir

yearb=2022
yeare=2023

for hemi in south 
do
  for month in `seq 1 12`
  do
    month=$(printf "%02d" $month)

    for year in `seq $yearb $yeare`
    do

      for day in `seq 1 31`
      do

        thisDate=`date -j  -f "%F" ${year}-${month}-${day} +"%Y%m%d"`
        threDate=`date -j  -f "%F" 2016-04-01 +"%Y%m%d"`


        if [[ $thisDate -ge $threDate ]]
        then
          sensor="f18"
        else
          sensor="f17"
        fi
        
        dday=$(printf "%02d" $day)
  
        #wget -ncd --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0081.002/${year}.${month}.${dday}/nt_${year}${month}${dday}_${sensor}_nrt_s.bin -P $outdir
	url=https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0081.002/${year}.${month}.${dday}/NSIDC0081_SEAICE_PS_S25km_${year}${month}${dday}_v2.0.nc

	wget -ncd --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off $url -P $outdir
      done
    done
  done
done
# Download grid files

wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25lats_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25lons_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25area_v3.dat -P ${outdir}


wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25lats_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25lons_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25area_v3.dat -P ${outdir}
