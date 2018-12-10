#!/bin/bash
set -o nounset
set -o errexit
set -x
# Downloads the SIPN South 2018 forecast data at the appropriate location
#
# F. Massonnet, June 2018
#               December 2019: update to download 2018-2019 data


url=https://nextcloud.cism.ucl.ac.be/s/gTL53xhjp4iQMM8


for mystring in 2018-2019
  do
  for file in `cat list_files_${mystring}.txt`
  do
    wget -N -c  ${url}/download?path=%2F${mystring}/netcdf/$file -O ../data/${mystring}/netcdf/$file
  done
done

