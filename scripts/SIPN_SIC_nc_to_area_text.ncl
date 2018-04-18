

;convert a SIC netcdf file to csv files in SIPN-south 'standard' format
;used mostly for actual observations

begin

  src  = "OSI-401-b_000"
;src = "NSIDC-0081_000"
  year = 2018


;***********************
;read SIC and grid area
;*************************


  fi = addfile(year+"/"+src+"_concentration.nc","r")
  
  SIC = fi->siconc
  A   = fi->areacello

  ;convert concentration from % to fraction
  if (max(SIC).gt.10.) then
    SIC = SIC/100.
  end if


;****************************
;calculate total area
;write to file
;*************************

  SIC = SIC * conform(SIC,A,(/1,2/)) * 1.e-12   ; Area/gid box (10^6 km^2)

  SIA = dim_sum_n(SIC, (/1,2/))

  SIAstr = str_join(sprintf("%3.2f", SIA), ",")
  asciiwrite(year+"/"+src+"_total-area.txt", SIAstr) 




;*********************************
;calculate SIA by longitude bin 
;write to file
;********************************


  nlon = 36
  olon = ispan(5,355,10)


  nday  = dimsizes(SIC&time)
  SIA  := new((/nlon,nday/), typeof(SIC))


  lon = fi->longitude

  ;-180-to-180 longitudes
  if (min(lon).lt.0.) then
    olon = where(olon.gt.180, olon-360, olon)
  end if


  rank = dimsizes(dimsizes(lon))

  if (rank.eq.1) then ;rectilinear grid
    
    do xx = 0, nlon-1
      SIA(xx,:) = dim_sum_n(SIC(:,:,{olon(xx)-5:olon(xx)+5}), (/1,2/))
    end do

  else                 ;curvilinear grid

    lon1d = ndtooned(lon)
    SIC2d = reshape(SIC, (/nday,dimsizes(lon1d)/))

    do xx = 0, nlon-1 
      rInd     :=ind(lon1d.gt.(olon(xx)-5).and.lon1d.le.(olon(xx)+5))
      SIA(xx,:) = dim_sum(SIC2d(:,rInd))
    end do

  end if

   

;write to file
  SIAstr := new(nlon, string)
  
  do xx = 0, nlon-1
    SIAstr(xx) = str_join(sprintf("%3.2f", SIA(xx,:)), ",")
  end do
  asciiwrite(year+"/"+src+"_regional-area.txt", SIAstr) 



end 