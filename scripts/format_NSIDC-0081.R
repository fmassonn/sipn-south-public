#!/usr/bin/env Rscript
rm(list=ls())
library(ncdf4)
# Francois Massonnet, April 2016
# Updated for TECLIM, May 2017
#
# Reads binaries from NASA nasateam near real time sea ice concentration and converts them in NCDF
#
#

yearb <- 2015          # The first year to process
yeare <- 2019          # The last  year to process

TECLIMDIR <- Sys.getenv(x = "TECLIM_CLIMATE_DATA")

rootdir <- paste(TECLIMDIR, "/obs/ice/siconc/NSIDC/NSIDC-0081/raw/", sep = "/") # The base directory

outdir <- paste(TECLIMDIR, "/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/", sep = "/") # Where to save the files

# Create if not exists
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


scalef_lon_lat <- 1e-5 # The lon and lat are magnified by 100,000 
scalef_are     <- 1e3  # The area files are 1000 times larger than they are supposed to be
                       # and are meant to be in km2. Thus they are in 1e-3 km2 = 1e3 m2.
                       # To retrieve m2 we have to multiply the data by 1000.

scalef_conc    <- 1 / 250 # The data are 250 times larger than they are supposed to be. And they
                       # are meant to be expressed in % [0-100]. Thus scale them by 1000.
nbytes_lon_lat <- 4    # The lon and lat are stored as 4-byte int
nbytes_are     <- 4    # The area      are stored as 4-byte int
nbytes_sic     <- 1    # The data      is stored as 4-byte int
nbytes_header  <- 300  # The number of bytes allocated to header
read_mode <- "little"  # not endian convention

nyear <- yeare - yearb + 1

# Count how many days
nt <- 0
for (year in seq(yearb, yeare)) {
  nt <- nt + 365
  isleap <- ((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)
  if (isleap) {
    nt <- nt + 1
  }
}

scaling <- 100.0 # convert fractional to % (to be CMIP compliant)

for (hemi in c("north","south")){
  print(hemi)
  # 1. Read the grid
  # ----------------
  if (hemi == "north"){
    ny <- 304
    nx <- 448
    file_lat <- paste(rootdir,"psn25lats_v3.dat",sep="/")
    file_lon <- paste(rootdir,"psn25lons_v3.dat",sep="/")
    file_are <- paste(rootdir,"psn25area_v3.dat",sep="/")

  } else {
    if (hemi == "south"){
      ny <- 316
      nx <- 332
      file_lat <- paste(rootdir,"pss25lats_v3.dat",sep="/")
      file_lon <- paste(rootdir,"pss25lons_v3.dat",sep="/")
      file_are <- paste(rootdir,"pss25area_v3.dat",sep="/")

    } 
  }

  sic_all = array(NA, dim = c(nx, ny, nt))

  # LAT
  to.read <- file(file_lat,"rb")
  lat     <- readBin(to.read,integer(),n=ny*nx,size=nbytes_lon_lat,endian=read_mode)
  lat     <- lat*scalef_lon_lat
  lat     <- matrix(lat,nrow=nx,byrow=TRUE)
  close(to.read)

  # LON
  to.read <- file(file_lon,"rb")
  lon     <- readBin(to.read,integer(),n=ny*nx,size=nbytes_lon_lat,endian=read_mode)
  lon     <- lon*scalef_lon_lat
  lon     <- matrix(lon,nrow=nx,byrow=TRUE)
  close(to.read)

  # AREA
  to.read <- file(file_are,"rb")
  are     <- readBin(to.read,integer(),n=ny*nx,size=nbytes_are,endian=read_mode)
  are     <- are*scalef_are
  are     <- matrix(are,nrow=nx,byrow=TRUE)
  close(to.read)
  # Read the data now
  # -----------------
  counter <- 0
  for (year in seq(yearb,yeare)){
    for (m in seq(1,12)){
      mm <- sprintf("%02d",m)

      # Sensor name
      if ((year == 2016 & m >= 4) | (year >= 2017)) {
        inst <- "f18"
      } else {
        inst <- "f17"
      }

      if (m == 1 | m == 3 | m == 5 | m == 7 | m == 8 | m == 10 | m == 12) {
        nday <- 31
      } else {
        if (m == 4 | m == 6 | m == 9 | m == 11) {
          nday <- 30
        } else {
          isleap <- ((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0)
          if (isleap) {
            nday <- 29
          } else {
            nday <- 28
          }
        }
      }
      
      for (day in seq(1, nday)) {
        counter <- counter + 1
        dd <- sprintf("%02d", day)
        filebin <- paste(rootdir,"/","nt_",toString(year),mm,dd,"_",inst,"_nrt_",substr(hemi,1,1),".bin",sep="")

        if (! file.exists(filebin)){
          print(paste("File ",filebin," not found",sep=""))
        } else {
          # Process it
          print(paste("Processing ",filebin,sep=""))
          to.read <- file(filebin,"rb")
          # There are first 300 bytes assigned to documentation (see http://nsidc.org/data/nsidc-0051)
          sic1d     <- readBin(to.read,integer(),n=nbytes_header+ny*nx,size=nbytes_sic,signed=FALSE)
          close(to.read)
          sic1d     <- sic1d[seq(nbytes_header+1, length(sic1d))]*scalef_conc

          sic     <- matrix(sic1d,nrow=nx,byrow=TRUE)
          # Put values of 1.004 (NP hole) to 1.0 and 1.2 to -1e24
          sic[sic==1.004] <- 1.0
          sic[sic>1.0]  <- -1e24

          sic_all[, , counter] = sic * scaling # to put in %
       }
     }
   } 
  } # year

  # Create time axis
  delta <- as.Date(paste(toString(yearb), "-01-01", sep = "")) - as.Date("1850-01-01")
  
  # Create mask based on last sic
  mask <- array(0, dim = c(nx, ny))
  # Look for the min over time
  tmp <- apply(sic_all, c(1, 2), max, na.rm = T)
  # if min < 0, mask, otherwise, 100
  mask[tmp >= 0] <- 100.0

  # Save as NetCDF
  fileout <- paste(outdir,"/siconc_SIday_NSIDC-0081_r1i1p1_",toString(yearb),"0101-",toString(yeare),"1231_",substr(hemi,1,1),"h.nc", sep="")
  dimt <- ncdim_def('time','',seq(1,nt), unlim = TRUE)
  dimy <- ncdim_def('y'   ,'',seq(1,ny))
  dimx <- ncdim_def('x'   ,'',seq(1,nx))

  vartime<- ncvar_def("time", paste('days since ','1850-01-01', sep =""), list(dimt), "time")
  varlat <- ncvar_def("latitude","degN",list(dimx,dimy),-1e26,"Latitude")
  varlon <- ncvar_def("longitude","degE",list(dimx,dimy),-1e26,"Longitude")
  varsic <- ncvar_def("siconc","%",list(dimx,dimy,dimt),-1e26,"Sea ice concentration")
  varare <- ncvar_def("areacello","m^2", list(dimx,dimy),-1e26, "cell_area")
  varmsk <- ncvar_def("sftof","%",list(dimx,dimy),-1e26,"sea_area_fraction")
       
  fid <- nc_create(fileout,list(varlat,varlon,varare,varmsk,varsic) )

  ncvar_put(fid,vartime,seq(1,nt) + delta - 1)
  ncatt_put(fid,"time","units", "days since 1850-01-01")
  ncatt_put(fid,"time","calendar", "standard")
  ncatt_put(fid,"time","axis","T")

  ncvar_put(fid,varlat,lat)
  ncvar_put(fid,varlon,lon)
  ncvar_put(fid,varare,are)
  ncvar_put(fid,varmsk,mask)
  ncvar_put(fid,varsic,sic_all)

  nc_close(fid)

  print(paste(fileout, " CREATED", sep = ""))
} # hemi


