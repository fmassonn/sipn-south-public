#!/bin/bash
# assembles figures

set -o nounset
set -o errexit

cd ../figs

          ["MetOffice",  range(1, 42 + 1), [0.00, 0.64, 0.85],   "$^d$"   ],  \
          ["ucl",        range(1, 10 + 1), [0.00, 0.84, 0.95],   "$^d$"   ],  \
          ["CNRM",       range(1, 51 + 1), [0.00, 0.94, 0.95],   "$^d$"   ],  \
          ["NicoSun",    range(1, 3 + 1),  [0.89, 0.42, 0.04],   "$^s$"   ],  \
          ["barreira",   [1]             , [0.89, 0.62, 0.24],   "$^s$"   ],  \
          ["Lamont",     [1],              [0.89, 0.82, 0.44],   "$^{s,i}$"], \
          ["cmcc"  ,     range(1, 50 + 1), [0.89, 0.32, 0.44],   "$^d$"   ],  \

# The monthly mean
convert +append CNRM_ens-mean_concentration_Feb2021-mean.png \
                MetOffice_ens-mean_concentration_Feb2021-mean.png \
                Lamont_ens-mean_concentration_Feb2021-mean.png
                X.png

convert +append NicoSun_ens-mean_concentration_Feb2021-mean.png \
                ucl_ens-mean_concentration_Feb2021-mean.png \
                barreira_ens-mean_concentration_Feb2021-mean.png
                Y.png

convert         cmcc_ens-mean_concentration_Feb2021-mean.png \
                Z.png

convert -append X.png Y.png Z.png monthly_mean.png
rm -f X.png Y.png Z.png 

# Movie

for day in `seq 63 90`
do
  dd=$(printf "%02d" $day)
  echo $dd
  # The probabilities
  convert +append nasa-gmao_prob-15_concentration_d${dd}.png \
                  MetOffice_prob-15_concentration_d${dd}.png \
                  ${dd}_X.png

  convert +append NicoSun_prob-15_concentration_d${dd}.png \
                  ucl_prob-15_concentration_d${dd}.png     \
                  ${dd}_Y.png
 
  convert -append ${dd}_X.png ${dd}_Y.png ${dd}.png

  rm -f ${dd}_X.png ${dd}_Y.png 
done

convert -delay 50 ??.png probability.gif
