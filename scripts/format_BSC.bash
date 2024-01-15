#!/bin/bash
#
# Converts EC-Earth decadal forecasts initialized on Nov 1 to SIPN South compliant output

# The input data must be in the current directory

set -o nounset
set -o errexit
set -x
yearb=2023
yeare=$(( $yearb + 1 ))
forecastingSeason=${yearb}-${yeare}

indir="../data/${forecastingSeason}/raw/BSC/"

# Month of initialization
monb=11

# Month of end
mone=10

nmemb=50

modName=BSC

# Where to take the geometrical parameters from
refGridFile=/Users/massonnetf/CLIMDATA/grid/mesh_mask_nemo.N3.6_ORCA1L75.nc


for imemb in `seq 1 $nmemb`
do
  iXXX=$(printf "%03d" $imemb)
  inFile=$indir/siconc_SImon_EC-Earth3_dcppB-forecast_s${yearb}-r${imemb}i4p1f1_gn_${yearb}${monb}-${yeare}${mone}.nc
  outFile=${modName}_${iXXX}_concentration_monthly.nc

  # Extract months December to Feb and dump in file; also ignore the variables and dimensions of no interest 
  # 
  ncks -F -O -d time,2,4 -v siconc,longitude,latitude $inFile $outFile
  

  # Extract geometrical data
  ncks -F -O -v tmaskutil,e1t,e2t $refGridFile tmpFile.nc
  # Remove time dimension
  ncwa -O -a t tmpFile.nc tmpFile.nc
  
  # Convert mask to %
  ncap2 -A -s "sftof=tmaskutil * 100" tmpFile.nc tmpFile.nc
  ncap2 -A -s "areacello=e1t*e2t"     tmpFile.nc tmpFile.nc

  ncks -F -A -v areacello,sftof tmpFile.nc $outFile
done
