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
set -x

yearb=2022
monb=1

yeare=2022
mone=12

ftype="multi" # multi (= operational, OSI-401-b) 

rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/OSI-SAF/OSI-401-b/raw

mkdir -p $outdir

#------------------------

for year in `seq $yearb $yeare`
do
  firstMonth=1
  lastMonth=12
  if [[ $year == $yearb ]]
  then
    firstMonth=$monb
  fi
  if [[ $year == $yeare ]]
  then
    lastMonth=$mone
  fi

  for month in `seq $firstMonth $lastMonth`
  do
    month=$(printf "%02d" $month)
    rootaddress="ftp://osisaf.met.no/archive/ice/conc/"
    wget -N -c $rootaddress/${year}/${month}/ice_conc_sh_polstere-100_${ftype}_${year}${month}??1200.nc -P $outdir
  done # month
done # year
