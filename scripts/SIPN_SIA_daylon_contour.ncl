


;plot predicted sea ice area anomaly by day and longitude, for SIPN predictions


begin


  isObs  = True  ;read and plot observations?

  year  = 2018
  month = "Feb"
  nday  = 28
  nlon  = 36


  oType = "ps"
  diro  = "~/plots/ANTARCTIC/SIPN/"
  filo  = "SIPN_"+month+year+"_SIA_daylon"


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
  lon@long_name = "longitude"
  lon@units     = "degrees east"
  nlon = dimsizes(lon)
  ilon = fi->lon

  tmp = SIA(:,:nlon-1)

  lonmin = lon-(xres/2)
  lonmax = lon+(xres/2)

  do xx = 0, nlon-1
    ii        = ind(ilon.gt.lonmin(xx).and.ilon.le.lonmax(xx))
    tmp(:,xx) = dim_sum(SIA(:,ii))
  end do
  delete(ii)

;get daily-climatology
  day   = time(tInd)%100
  clim  = new((/nday,nlon/), typeof(tmp))
  
  do i = 0, nday-1
    ii        = ind(day.eq.(i+1))
    clim(i,:) = dim_avg_n(tmp(ii,:), 0)
  end do
  delete(ii)
  clim = clim * 1.e-6

  delete([/SIA,time,mon,day,ilon,tInd,tmp/])


;*************************************
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
;read data
;****************************

  tmp = new((/nfil,nlon,nday/), float)

  do i = 0, nfil-1
    tmp(i,:,:) = readAsciiTable(fili(i),nday,"float",0)
  end do


  tmp!0 = "src"
  tmp!1 = "lon"
  tmp!2 = "day"


  SIA   = tmp(src|:,day|:,lon|:)
  SIA   = SIA - conform(SIA,clim,(/1,2/))


  delete([/fili,tmp/])



;******************************
;get ensemble-means and spread
;******************************

  EnsMn     = SIA(:nsrc-1,:,:)
  EnsMn&lon = lon  
  EnsMn&day = ispan(1,nday,1)

  sig       = EnsMn
  
  do i = 0, nsrc-1

    if (Nmem(i).gt.1) then
      EnsMn(i,:,:) = (/ dim_avg_n( SIA(memList[i],:,:), 0) /)

      ;find where all members match with ensemble mean
      tst          = conform(SIA(memList[i],:,:),EnsMn(i,:,:),(/1,2/)) * SIA(memList[i],:,:) ;+ve where member matches sign of EnsMn
      sig(i,:,:)   = where(dim_min_n(tst,0).gt.0.,1.,0.)      
      delete(tst)
    else
      EnsMn(i,:,:) = (/ SIA(memList[i],:,:) /)
    end if

  end do
  
  
  delete([/SIA,memList/])
  


;***********************
;read observations
;************************

  if (isObs) then

    fili         = systemfunc("ls "+diri+"OBS/*_???_regional-area.txt")
    nObs         = dimsizes(fili)

    filsplt      = str_split_csv(fili, "_",0)           ;split up filepaths
    filsplt(:,0) = str_sub_str( filsplt(:,0), diri+"OBS/", "") ;strip directory off filepaths
    osrc         = get_unique_values(filsplt(:,0)) 

    tmp = new((/nObs,nlon,nday/), float)

    do i = 0, nObs-1
      tmp(i,:,:) = readAsciiTable(fili(i),nday,"float",0)
    end do


    tmp!0 = "src"
    tmp!1 = "lon"
    tmp!2 = "day"

    Obs   = tmp(src|:,day|:,lon|:)
    Obs   = Obs - conform(Obs,clim,(/1,2/))


    delete([/tmp,fili,filsplt/])

  end if
  
 


;**************************
;plot
;************************


  wks = gsn_open_wks(oType, diro+filo)

  
  res                        = True
  res@vpHeightF              = .6
  res@vpWidthF               = .9
  res@gsnFrame               = False
  res@gsnDraw                = False

  res@cnLinesOn              = False
  res@cnFillOn               = True
  res@cnFillPalette          = "BlueDarkRed18"
  symMinMaxPlt(EnsMn, 18, False, res)

  res@lbLabelBarOn           = False
  res@cnInfoLabelOn          = False
  res@cnLineLabelsOn         = False

;significance resources
  ores                       = True
  ores@gsnFrame              = False
  ores@gsnDraw               = False
  ores@gsnLeftString         = ""
  ores@gsnRightString        = ""
  ores@cnFillOn              = True
  ores@cnFillPattern         = 17
  ores@cnFillColors          = (/"transparent","black"/)
  ores@cnLevelSelectionMode  = "ExplicitLevels"
  ores@cnLevels              = .5
  ores@cnMonoFillPattern     = True
  ores@lbLabelBarOn          = False
  ores@cnInfoLabelOn         = False
  ores@cnLineLabelsOn        = False

;panel resources
  resP                       = True
  resP@gsnPanelLabelBar      = True
  resP@lbOrientation         = "Vertical"
  resP@gsnMaximize           = True

  

;plot
  if (isObs) then
    nplot  = nsrc+nObs
    RowCol = (/nplot/nObs,nObs/) 
  else
    nplot = nsrc
    RowCol = (/nplot/2,2/)
  end if  


  plot   = new(nplot, graphic)

  do i = 0, nsrc-1

    res@gsnLeftString  = src(i)
    res@gsnRightString = Nmem(i)+" members"

    plot(i) = gsn_csm_contour(wks,EnsMn(i,:,:), res)
    
    if (Nmem(i).gt.1) then
      oplot = gsn_csm_contour(wks,sig(i,:,:), ores)
      overlay(plot(i), oplot)
    end if

  end do

  if (isObs) then
    res@gsnRightString = ""

    do i = 0, nObs-1
      res@gsnLeftString = osrc(i)
      plot(i+nsrc)      = gsn_csm_contour(wks,Obs(i,:,:), res)
    end do
  end if



;add base map of Antarctic continent to bottom of each column
  mpres             = True
  mpres@gsnDraw     = False
  mpres@gsnFrame    = False
  mpres@mpLandFillColor = "black"
  mpres@mpMaxLatF   = -55.
  mpres@mpMinLatF   = -90.
  mpres@mpCenterLonF = 180.

  ;annotation resources
  ares                     = True
  ares@gsnAttachPlotsXAxis = True


  map  = new(RowCol(1), graphic)
  anno = new(RowCol(1), graphic)

  do i = 0, dimsizes(map)-1
    map(i)  = gsn_csm_map_ce(wks, mpres)
    anno(i) = gsn_attach_plots(plot(nplot-i-1), map(i), False,ares)
  end do



;panel
  gsn_panel(wks,plot,RowCol , resP)


end