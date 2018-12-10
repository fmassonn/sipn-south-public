#!/bin/bash
# assembles figures

set -o nounset
set -o errexit

cd ../figs

for day in `seq 1 28`
do
  dd=$(printf "%02d" $day)
  echo $dd
  # The probabilities
  convert +append nrl_prob-15_concentration_d${dd}.png \
                  nasa-gmao_prob-15_concentration_d${dd}.png \
                  MetOffice_prob-15_concentration_d${dd}.png \
                  ${dd}_X.png

  convert +append ucl_prob-15_concentration_d${dd}.png \
                  emc_prob-15_concentration_d${dd}.png \
                  ${dd}_Y.png

  convert -append ${dd}_X.png ${dd}_Y.png ${dd}.png

  rm -f ${dd}_X.png ${dd}_Y.png
done

convert -delay 50 ??.png probability.gif

# The monthly mean
convert +append nrl_ens-mean_concentration_mon-mean.png \
                nasa-gmao_ens-mean_concentration_mon-mean.png \
                MetOffice_ens-mean_concentration_mon-mean.png \
                X.png

convert +append ucl_ens-mean_concentration_mon-mean.png \
                emc_ens-mean_concentration_mon-mean.png \
                Y.png

convert -append X.png Y.png monthly_mean.png
rm -f X.png Y.png
