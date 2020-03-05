#!/bin/bash
#
# F. Massonnet
# Update obs verif data

set -o nounset
set -o errexit
set -x 

#. ./retrieve_NSIDC-0081.bash

#. ./retrieve_OSI-401b.bash

module load R
Rscript format_NSIDC-0081.R

#Storm
#module load Python/3.6.1-intel-2018
python ./format_OSI-401b.py


module load NCO
ncks -F -O -d time,335,424 ../data/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/siconc_SIday_NSIDC-0081_r1i1p1_20190101-20201231_sh.nc ../data/2019-2020/netcdf/NSIDC-0081_000_concentration.nc

ncks -F -O -d time,335,424 ../data/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_20190101-20201231_sh.nc ../data/2019-2020/netcdf/OSI-401-b_000_concentration.nc

./obs2CSV.py

