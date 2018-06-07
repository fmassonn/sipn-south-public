#!/bin/bash
set -o nounset
set -o errexit

# Downloads the SIPN South 2018 forecast data at the appropriate location
#
# F. Massonnet, June 2018

url=https://nextcloud.cism.ucl.ac.be/s/gTL53xhjp4iQMM8

for file in `cat list_files.txt`
do
  wget -N -c ${url}/download?path=%2F2018/netcdf/$file -O $file
done

