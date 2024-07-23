#!/bin/bash
# assembles figures

set -o nounset
set -o errexit
set -x

cd ../figs

# The monthly mean
convert +append AWI-CPS_ens-mean_concentration_Feb2024-mean.png \
                gfdl_ens-mean_concentration_Feb2024-mean.png \
                IOCAS-SIPNet_ens-mean_concentration_Feb2024-mean.png \
                X.png

convert +append Lamont_ens-mean_concentration_Feb2024-mean.png \
                MetOffice_ens-mean_concentration_Feb2024-mean.png \
                NicoSun_ens-mean_concentration_Feb2024-mean.png \
                Y.png

convert +append SYSU-SML-MLM_ens-mean_concentration_Feb2024-mean.png \
                SYSU-SML-KNN_ens-mean_concentration_Feb2024-mean.png \
                cmcc_ens-mean_concentration_Feb2024-mean.png \
                Z.png
convert +append emc_ens-mean_concentration_Feb2024-mean.png \
                T.png

convert -append X.png Y.png Z.png T.png monthly_mean.png
rm -f X.png Y.png Z.png 

exit

# Movie

for day in `seq 63 90`
do
  dd=$(printf "%02d" $day)
  echo $dd
  # The probabilities
  convert +append MetOffice_prob-15_concentration_d${dd}.png \
                  ucl_prob-15_concentration_d${dd}.png \
                  CNRM_prob-15_concentration_d${dd}.png \
                  ${dd}_X.png

  convert +append NicoSun_prob-15_concentration_d${dd}.png  \
                  Lamont_prob-15_concentration_d${dd}.png   \
                  barreira_prob-15_concentration_d${dd}.png \
                  ${dd}_Y.png
 
  convert +append cmcc_prob-15_concentration_d${dd}.png \
                  gfdl_prob-15_concentration_d${dd}.png \
                  SYSU_prob-15_concentration_d${dd}.png \
                  ${dd}_Z.png

  convert -append ${dd}_X.png ${dd}_Y.png ${dd}_Z.png ${dd}.png

  rm -f ${dd}_X.png ${dd}_Y.png ${dd}_Z.png
done

convert -delay 50 ??.png probability.gif
