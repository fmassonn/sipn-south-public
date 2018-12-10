#!/bin/bash
# assembles figures

set -o nounset
set -o errexit

cd ../figs

for day in `seq 63 90`
do
  dd=$(printf "%02d" $day)
  echo $dd
  # The probabilities
  convert +append nrl_prob-15_concentration_d${dd}.png \
                  nasa-gmao_prob-15_concentration_d${dd}.png \
                  MetOffice_prob-15_concentration_d${dd}.png \
                  ${dd}_X.png

  convert +append Lamont_prob-15_concentration_d${dd}.png \
                  Nico-Sun_prob-15_concentration_d${dd}.png \
                  ${dd}_Y.png

  convert -append ${dd}_X.png ${dd}_Y.png ${dd}.png

  rm -f ${dd}_X.png ${dd}_Y.png
done

#convert -delay 50 ??.png probability.gif

# The monthly mean
convert +append nrl_ens-mean_concentration_Feb2019-mean.png \
                nasa-gmao_ens-mean_concentration_Feb2019-mean.png \
                MetOffice_ens-mean_concentration_Feb2019-mean.png \
                X.png

convert +append Lamont_ens-mean_concentration_Feb2019-mean.png \
                Nico-Sun_ens-mean_concentration_Feb2019-mean.png \
                Y.png

convert -append X.png Y.png monthly_mean.png
rm -f X.png Y.png
