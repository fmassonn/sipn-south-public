#!/bin/bash
#
# F. Massonnet
# Downloads sea ice concentration of NASA Team near real time (NSIDC 0081)
set -o nounset
#set -o errexit
set -x 
rootdir=$TECLIM_CLIMATE_DATA
outdir=${rootdir}/obs/ice/siconc/NSIDC/NSIDC-0081/raw/

mkdir -p $outdir

yearb=2022
yeare=2022

listNotAvailable=( "20210220" "20210221" "20220901" "20220902" "20220903" "20220904" "20220905" "20220906" "20220907" "20220908" "20220909" "20220910" "20220911" "20220912" "20220913" "20220914" "20220915" "20220916" "20220917" "20220918" "20220919" "20220920" "20220921" "20220922" "20220923" "20220924" "20220925" "20220926" "20220927" "20220928" "20220929" "20220930" ) 

for hemi in south 
do
  for month in `seq 11 12`
  do
    month=$(printf "%02d" $month)

    for year in `seq 2015 2022`
    do
      if [[ $month ==  "01" ]] || [[ $month == "03" ]] || [[ $month == "05" ]] || [[ $month == "07" ]] || [[ $month == "08" ]] || [[ $month == "10" ]] || [[ $month == "12" ]]
      then
        nDay=31
      elif [[ $month ==  "04" ]] || [[ $month == "06" ]] || [[ $month == "09" ]] || [[ $month == "11" ]]
      then
        nDay=30
      else
        nDay=28
      fi

      date -d $yearb-02-29 &>/dev/null && leap="T" || leap="F"
      if [[ $leap == "T" ]]
      then
        nDay=29
      fi

      for day in `seq 1 $nDay`
      do

        echo $month
        thisDate=`date -j  -f "%F" ${year}-${month}-${day} +"%Y%m%d"`
        threDate=`date -j  -f "%F" 2016-04-01 +"%Y%m%d"`

        skipDate=False

        for item in ${listNotAvailable[@]}
        do
          if [[ $thisDate == $item ]]
          then
            skipDate=True
          fi 
        done
     
        echo $listNotAvailable

        if [[ $thisDate -ge $threDate ]]
        then
          sensor="f18"
        else
          sensor="f17"
        fi
        
        dday=$(printf "%02d" $day)
        if [[ $skipDate == False ]]
        then
          wget -ncd --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0081.001/${year}.${month}.${dday}/nt_${year}${month}${dday}_${sensor}_nrt_s.bin -P $outdir
        fi
      done

    done
  done
done

# Download grid files

wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25lats_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25lons_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/psn25area_v3.dat -P ${outdir}


wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25lats_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25lons_v3.dat -P ${outdir}
wget -nc ftp://sidads.colorado.edu/pub/DATASETS/brightness-temperatures/polar-stereo/tools/geo-coord/grid/pss25area_v3.dat -P ${outdir}
