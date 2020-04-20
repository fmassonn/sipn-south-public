

;plot monthly statistic of predicted sea ice area by longitude, for SIPN predictions


begin

  stat   = "mean"  ;statistic to plot - "mean" or "min"
  
  isAnom = True  ;plot climatological anomaly (True) or abs value (False)
  isObs  = True  ;plot Obs ?
  colfile  = "~/scripts/ANTARCTIC/SIPN/SIPN_colorfile.csv" ;filepath of standard colormap

  callyear = 2020
  year     = 2020
  month    = "Feb"
  nday     = 28 ;   days_in_month(year,mInd)   ;number of days in month
  mInd     = 2
  date0    = (/2019,12,1/)  ;first date of prediction file, formatted year,month,day

  nlon     = 36

  sNam  = (/"King Hakon VII","E. Antarctic","Ross","Amundsen","Bellingshausen","Weddell" /)
  nsec  = dimsizes(sNam)
  xmin  = (/360.,90.,160.,210.,260.,300./)


  oType = "ps"
  diro  = "~/plots/ANTARCTIC/SIPN/"
  filo  = "SIPN_"+month+year+"_"+stat+"_SIA_bylon-corrected"

;**************************
;get climatology
;*************************



  fi   = addfile("~/data/NSIDC/NSIDC_nasa_team_SH_icearea_daily.nc","r")
  time = fi->time
  mon  = (time%10000)/100
  tInd = ind(mon.eq.mInd)

  SIA  = fi->SIA(tInd,:)

;remap
  xres = 10
  lon  = ispan(5,355,xres)
  nlon = dimsizes(lon)
  ilon = fi->lon

  tmp = SIA(:,:nlon-1)

  lonmin = lon-(xres/2)
  lonmax = lon+(xres/2)

  do xx = 0, nlon-1
    ii        = ind(ilon.gt.lonmin(xx).and.ilon.le.lonmax(xx))
    tmp(:,xx) = dim_sum(SIA(:,ii))
  end do

;get daily-climatology
  day      = time(tInd)%100
  climtmp  = new((/nday,nlon/), typeof(tmp))
  
  do i = 0, nday-1
    ii          := ind(day.eq.(i+1))
    climtmp(i,:) = dim_avg_n(tmp(ii,:), 0)
  end do
 
  climtmp = climtmp * 1.e-6


;calculate statistic
  if (stat.eq."mean") then
    clim = dim_avg_n(climtmp,0)
  end if
  if (stat.eq."min") then
    clim = dim_min_n(climtmp,0)
  end if
  

  delete([/SIA,climtmp,time,mon,day,ilon,lonmin,lonmax,ii,tmp/])


;******************************************
;function to read SIPN csv formatted files
;******************************************

 undef("read_csv")
 function read_csv(fnam[1]:string)
;fnam : scalar string representing input filepath

 local dat, ndat, line, ncol,out,i

 begin

   dat    = asciiread(fnam,-1,"string")     ;read data as a string
   dsplit = str_split_csv(dat,",",0) ;split entries using comma
   out    = stringtofloat(dsplit)              ; convert to float
   return(out)
 end 



;******************
;get file lists
;find contributors and their ensembles
;**************************************

  diri         = "/Users/whobbs/data/SIPN/"+callyear+"/"
  fili         = systemfunc("ls "+diri+"*_???_regional-area.txt")
  nfil         = dimsizes(fili)


  ;find contributors
  filsplt      = str_split_csv(fili, "_",0)           ;split up filepaths
  filsplt(:,0) = str_sub_str( filsplt(:,0), diri, "") ;strip directory off filepaths
  src          = get_unique_values(filsplt(:,0))      ;get list of contributors
  nsrc         = dimsizes(src)

  ;get indices of members for each group
  memList      = NewList("lifo")
  Nmem         = new(nsrc, integer)

  do i = 0, nsrc-1
    ListAppend(memList, ind(filsplt(:,0).eq.src(i)) )
    Nmem(i) = dimsizes(memList[i])
  end do

  delete(filsplt)

;****************************
;read forecast data
;****************************
  

;find indices of analysis period, based on perdiction file start date

  ttmp     = yyyymmdd_time(date0(0),year,"integer")
  ind0     = ind(ttmp.eq.date0(0)*10000+date0(1)*100+date0(2))
  stdat    = ind(ttmp.eq.year*10000+mInd*100+1) -ind0
  endat    = stdat + nday -1



;read data
  data = new((/nfil,nlon,nday/), float)

  do i = 0, nfil-1
    tmp         = read_csv(fili(i))
    data(i,:,:) = tmp(:,stdat:endat)
  end do




  delete([/fili,tmp,ttmp,ind0/])



;******************
;get monthly statistic
;******************8

  
  if (stat.eq."mean") then
    monMean = dim_avg(data)
  end if
  if (stat.eq."min") then
    monMean = dim_min(data)
  end if


  if (isAnom) then
    monMean = monMean - conform(monMean, clim, 1)
  end if

  delete(data)


;****************************
;read Obs
;****************************

  if (isObs) then

    fili         = systemfunc("ls "+diri+"OBS/*_???_regional-area.txt")
    nObs         = dimsizes(fili)

    filsplt      = str_split_csv(fili, "_",0)           ;split up filepaths
    filsplt(:,0) = str_sub_str( filsplt(:,0), diri+"OBS/", "") ;strip directory off filepaths
    osrc         = get_unique_values(filsplt(:,0)) 

    Obs = new((/nObs,nlon,nday/), float)

    do i = 0, nObs-1
      tmp         = read_csv(fili(i))
      Obs(i,:,:)  = tmp(:,stdat:endat)
    end do

    delete(fili)


    if (stat.eq."mean") then
      oMean = dim_avg(Obs)
    end if
    if (stat.eq."min") then
      oMean = dim_min(Obs)
    end if


    if (isAnom) then
      oMean = oMean - conform(oMean, clim, 1)
    end if

    delete(Obs)

  end if


;******************************
;get ensemble-means and spread
;******************************
  
  if (isObs) then
    EnsMn          = new((/nsrc+nObs,nlon/), typeof(monMean))
    EnsMn(nsrc:,:) = oMean
  else
    EnsMn = new((/nsrc,nlon/), typeof(monMean))
  end if

  StDev = new((/nsrc,nlon/), typeof(monMean))
  Range = new((/nsrc,2,nlon/), typeof(monMean))

  
  do i = 0, nsrc-1

    if (Nmem(i).gt.1) then

      EnsMn(i,:)   = dim_avg_n( monMean(memList[i],:), 0)
      StDev(i,:)   = dim_stddev_n( monMean(memList[i],:), 0)

      Range(i,0,:) = dim_min_n( monMean(memList[i],:), 0)
      Range(i,1,:) = dim_max_n( monMean(memList[i],:), 0)

    else
      EnsMn(i,:)   = monMean(memList[i],:)
    end if

  end do


  ;regressions against Obs
  if (isObs) then
    reg = regCoef(EnsMn(nsrc,:),EnsMn) 
  end if

do i = 0, nsrc-1
print(src(i)+", "+reg(i))
end do

;**************************
;plot
;************************


  wks = gsn_open_wks(oType, diro+filo)

  
  res                        = True
  res@vpHeightF              = .6
  res@vpWidthF               = .9
  res@gsnMaximize            = True
  res@gsnFrame               = False
  res@gsnDraw                = False

  if (isAnom) then
    res@tiYAxisString          = "SIA anomaly (10~S~6~N~ km~S~2~N~)"
    res@trYMinF                = 0.1 * floor(10.*min((/min(monMean),min(oMean)/)))
  else
    res@tiYAxisString          = "SIA (10~S~6~N~ km~S~2~N~)"
    res@trYMinF                = 0.
  end if
  res@trYMaxF                = 0.1 * ceil(10.* max((/max(monMean),max(oMean)/)))
  res@gsnYRefLine            = 0.
  res@gsnYRefLineThicknessF  = 2.

  res@tiXAxisString          = "longitude (~S~o~N~E)"
  res@trXMinF                = 0
  res@trXMaxF                = 360
  res@gsnXRefLine            = xmin   ;sector lines
  res@gsnXRefLineThicknessF  = 2.
  res@gsnXRefLineDashPattern = 3.

  res@xyMonoDashPattern      = True




;colors 
  tmp    := asciiread(colfile,-1,"string")     ;read data as a string
  tmp   := str_split_csv(tmp(1:),",",0) ;split entries using comma

  cInd                      = new(nsrc, integer)
  do i = 0, nsrc-1
    cInd(i) = ind(tmp(:,0).eq.src(i))
  end do



  cmap                      = stringtofloat(tmp(cInd,1:)) ;nsrc, 3
  delete(tmp)

;plot ranges
  res@xyLineThicknessF      = 0.01              
      
  res@xyLineColors          = cmap
  res@gsnXYFillOpacities    = .3


  i=0
 ;initial plot needs to have a range, i.e. more than 1 member
  ; this ensure first plot has that
  do while (Nmem(i).eq.1) 
    i = i+1
  end do
  res@gsnXYFillColors       = cmap(i,:)  
  res@gsnLeftString         = "";month+" "+year
  plot                      = gsn_csm_xy(wks, lon, Range(i,:,:), res)
  res@gsnLeftString         = ""

  do i = i+1, nsrc-1
    if (Nmem(i).eq.1) then
      continue
    end if
    res@gsnXYFillColors       = cmap(i,:)
    oplot = gsn_csm_xy(wks, lon, Range(i,:,:), res)
    overlay(plot, oplot)
  end do
  delete([/res@gsnXYFillOpacities, res@gsnXYFillColors /])


;plot ensemble means

  
  if (isObs) then
    res@xyExplicitLegendLabels = array_append_record(src,osrc,0)
    res@xyLineThicknesses      = array_append_record(conform_dims(nsrc,3.5,-1),conform_dims(nObs,5.,-1),0)
    lcols                      = new((/nsrc+nObs,3/),float) 
    lcols(:nsrc-1,:)           = cmap
    lcols(nsrc:,:)             = 0.
    res@xyLineColors          := lcols
    res@xyDashPatterns         = array_append_record(conform_dims(nsrc+1,0,-1),(/14,15,16/),0)
    res@xyMonoDashPattern      = False

  else
    res@xyExplicitLegendLabels = src
    res@xyLineThicknessF       = 3.5
    res@xyDashPattern          = 0
  end if


 ;legend  
  res@pmLegendDisplayMode      = "Always"
  res@lgPerimOn                = True

  res@pmLegendWidthF           = 0.22
  res@pmLegendHeightF          = 0.24
  res@lgLabelFontHeightF       = 0.018
  res@pmLegendSide             = "Top"
  res@pmLegendParallelPosF     = 0.19  ;x-position
  if (isAnom) then
    res@pmLegendOrthogonalPosF = -1.03 ;y-position
  else
    res@pmLegendOrthogonalPosF = -0.5 ;y-position
  end if


  oplot =  gsn_csm_xy(wks, lon, EnsMn, res)
  overlay(plot, oplot)



;plot climatology

  if (.not.isAnom) then
    res@xyLineThicknessF      = 5.
    res@xyLineColors(0)       = (/0,0,0,1/)

    res@pmLegendDisplayMode    = "Never"
    oplot =  gsn_csm_xy(wks, lon, clim, res)
    overlay(plot, oplot)
  end if

;add sector markers
  xmax  = array_append_record(xmin(1:),xmin(0),0)
  xpos = dim_avg_n((/xmax,xmin/),0)
  xpos = where(xmin.gt.xmax,xpos-180., xpos) 
  
  ypos  = conform(xpos, res@trYMaxF - .05, -1)


  rest               = True
  rest@txFontHeightF = 0.012
  
  text               = gsn_add_text(wks, plot, sNam, xpos, ypos, rest)
  


  draw(plot)
  frame(wks)


end
