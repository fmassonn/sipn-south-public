#!/bin/bash
set -o nounset
set -o errexit

# Downloads the SIPN South 2018 forecast data at the appropriate location
#
# F. Massonnet, June 2018
#               December 2019: update to download 2018-2019 data


url=https://nextcloud.cism.ucl.ac.be/s/gTL53xhjp4iQMM8

for year in `seq 2017 2018`
  mystring="$year-$((year + 1))"
  for file in `cat list_files_${mystring}.txt`
  do
    wget -N -c ${url}/download?path=%2F${mystring}/netcdf/$file -O ../data/netcdf/${my_string}/$file
  done
done

