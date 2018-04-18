
;plot monthly statistic of predicted sea ice area by longitude, for SIPN predictions


begin

  stat   = "min"  ;statistic to plot - "mean" or "min"
  
  isAnom = True  ;plot climatological anomaly (True) or abs value (False)
  isObs  = True  ;plot Obs ?
  

  year  = 2018
  month = "Feb"
  nday  = 28
  nlon  = 36

  sNam  = (/"King Hakon VII","E. Antarctic","Ross","Amundsen","Bellingshausen","Weddell" /)
  nsec  = dimsizes(sNam)
  xmin  = (/360.,90.,160.,210.,260.,300./)


  oType = "ps"
  diro  = "~/plots/ANTARCTIC/SIPN/"
  filo  = "SIPN_"+month+year+"_"+stat+"_SIA_bylon"

;**************************
;get climatology
;*************************



  fi   = addfile("~/data/NSIDC/NSIDC_nasa_team_SH_icearea_daily.nc","r")
  time = fi->time
  mon  = (time%10000)/100
  tInd = ind(mon.eq.2)

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
  

  delete([/SIA,climtmp,time,mon,day,ilon,lonmin,lonmax,ii/])


;******************
;get file lists
;find contributors and their ensembles
;**************************************

  diri         = "/Users/whobbs/data/SIPN/"+year+"/"
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

  data = new((/nfil,nlon,nday/), float)

  do i = 0, nfil-1
    data(i,:,:) = readAsciiTable(fili(i),nday,"float",0)
  end do

  delete(fili)





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
      Obs(i,:,:) = readAsciiTable(fili(i),nday,"float",0)
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
    res@trYMinF                = 0.1 * floor(10.*min(monMean))
  else
    res@tiYAxisString          = "SIA (10~S~6~N~ km~S~2~N~)"
    res@trYMinF                = 0.
  end if
  res@trYMaxF                = 0.1 * ceil(10.*max(monMean))
  res@gsnYRefLine            = 0.
  res@gsnYRefLineThicknessF  = 2.

  res@tiXAxisString          = "longitude (~S~o~N~E)"
  res@trXMinF                = 0
  res@trXMaxF                = 360
  res@gsnXRefLine            = xmin   ;sector lines
  res@gsnXRefLineThicknessF  = 2.
  res@gsnXRefLineDashPattern = 3.

  res@xyMonoDashPattern      = True




;colors - this needs to be better
  cmap  = read_colormap_file("srip_reanalysis")
  ncol  = 12
  stcol = 2



;plot ranges
  res@xyLineThicknessF      = 0.01                     
  res@xyLineColors          = cmap(:nsrc*2-1:2,:)
  res@gsnXYFillOpacities    = .3



  i=1
  res@gsnXYFillColors       = res@xyLineColors(i)  
  plot                      = gsn_csm_xy(wks, lon, Range(i,:,:), res)

  do i = 2, nsrc-1
    res@gsnXYFillColors       = res@xyLineColors(i) 
    oplot = gsn_csm_xy(wks, lon, Range(i,:,:), res)
    overlay(plot, oplot)
  end do
  delete([/res@gsnXYFillOpacities, res@gsnXYFillColors /])


;plot ensemble means

  
  if (isObs) then
    res@xyExplicitLegendLabels = array_append_record(src,osrc,0)
    res@xyLineThicknesses      = array_append_record(conform_dims(nsrc,3.5,-1),conform_dims(nObs,5.,-1),0)
    lcols                      = new((/nsrc+nObs,4/),float) 
    lcols(:nsrc-1,:)           = res@xyLineColors
    lcols(nsrc:,:)             = conform_dims((/nObs,4/),(/0.,0.,0.,1/),1)
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
  res@pmLegendSide             = "Top"
  res@pmLegendParallelPosF     = 0.13  ;x-position
  if (isAnom) then
    res@pmLegendOrthogonalPosF = -1. ;y-position
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
