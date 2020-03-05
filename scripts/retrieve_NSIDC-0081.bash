#!/bin/bash
#
# F. Massonnet
# Downloads sea ice concentration of NASA Team near real time (NSIDC 0081)
set -o nounset
set -o errexit
set -x 
rootdir=../data/ #$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/NSIDC/NSIDC-0081/raw/

mkdir -p $outdir

for hemi in south 
do
  for month in 01 02 03 04 05 06 07 08 09 10 11 12
  do
    for year in `seq 2019 2020`
    do
      echo ""
      wget -nc ftp://sidads.colorado.edu/pub/DATASETS/nsidc0081_nrt_nasateam_seaice/${hemi}/nt_${year}${month}??_f??_nrt_?.bin -P ${outdir}
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
