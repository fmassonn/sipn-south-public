#!/usr/bin/python3
#
# Francois Massonnet francois.massonnet@uclouvain.be
# 
# Converts near-real time NSIDC sea ice concentration data
# (NSIDC-0081) to TECLIM compliant format
#

def formatData(dateStart, dateEnd):
  """ dateStart and dateEnd being to YYYYMMDD strings defining the period to format data
  """


  import numpy as np
  from netCDF4 import Dataset
  from datetime import datetime, timedelta
  import os
  import struct
  
  # The (machine-dependent) location of CLIMDATA folder
  data_dir = os.environ["TECLIM_CLIMATE_DATA"]
  
  # =========================
  
  scalefLat = 1e-5 # When reading binaries
  scalefLon = 1e-5 # When reading binaries
  scalefAre = 1e3  # When reading binaries

  hemi="sh"    # Hemisphere to proces ("nh" or "sh")
  
  # Where to read the data from
  rootdir = data_dir + "/obs/ice/siconc/NSIDC/NSIDC-0081/raw/"
  # Where to write the data to
  outdir =  data_dir + "/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/"
  
  
  # Check existence of paths
  if not os.path.exists(rootdir):
      sys.exit("Input directory does not exist")
  
  if not os.path.exists(outdir):
      os.makedirs(outdir)
  
  # ========
  # Input grid dimensions
  if hemi == "sh":
      ny, nx = 332, 316
  elif hemi == "nh":
      pass
  else:
      sys.exit("Hemishpere unknown")
  
  # Create time axis
  d1 = datetime.strptime(dateStart, "%Y%m%d")
  d2 = datetime.strptime(dateEnd  , "%Y%m%d")
  daterange = [d1 + timedelta(days=x) for x in range((d2-d1).days + 1)]
  nt = len(daterange)
  
  # Create numpy arrays with the output variables: siconc, longitude, latitude, time
  siconc       = np.empty((nt, ny, nx))
  siconc[:]    = np.nan
  longitude    = np.empty((ny, nx))
  longitude[:] = np.nan
  latitude     = np.empty((ny, nx))
  latitude[:]  = np.nan
  time         = np.empty((nt))
  time[:]      = np.nan
  cellarea     = np.empty((ny, nx))
  cellarea[:]  = np.nan
  mask         = np.empty((ny, nx))
  mask[:]      = np.nan
  
  jt = 0 # Time counter to write the data
  read_geom = True # boolean to read the geometric data only once
  
  for day in daterange:
      print(day)
      # Check file existence
      filein = rootdir + "/NSIDC0081_SEAICE_PS_" + hemi[0].upper() + "25km_" + day.strftime("%Y%m%d") + "_v2.0.nc"
      if not os.path.exists(filein):
          print("File " + filein + " not found")
      else:
          f = Dataset(filein, mode = "r")
          try:
            sic = f.variables["F18_ICECON"][:]
          except KeyError:
            try:
              sic = f.variables["F17_ICECON"][:]
            except KeyError:
              try:
                sic = f.variables["F16_ICECON"][:]
              except KeyError:
                stop()
  
          sic[sic > 1.0] = 0.0
          sic[sic < 0.0] = 0.0
  	
          siconc[jt, :, :] = sic[:] * 100 # Convert to %
  
          if read_geom:
              # We need to read those stupid binaries
              fileLat = rootdir + "/" + "pss25lats_v3.dat"
              thisArray = np.reshape(np.array(struct.unpack_from("<" + "i" * ny * nx, open(fileLat, mode="rb").read())), (ny, nx))
              lat = thisArray * scalefLat
  
              fileLon = rootdir + "/" + "pss25lons_v3.dat"
              thisArray = np.reshape(np.array(struct.unpack_from("<" + "i" * ny * nx, open(fileLon, mode="rb").read())), (ny, nx))
              lon = thisArray * scalefLon
  
              fileAre = rootdir + "/" + "pss25area_v3.dat"
              thisArray = np.reshape(np.array(struct.unpack_from("<" + "i" * ny * nx, open(fileAre, mode="rb").read())), (ny, nx))
              are  = thisArray * scalefAre
  
              msk  = np.full((ny, nx), 100.0)	
  
              latitude[:] = lat
              longitude[:] = lon
              cellarea[:] = are
              mask[:]     = msk
  
              read_geom = False
          f.close()
      jt += 1
  
  
  # Save the data
  # -------------
  
  # 5. Save as NetCDF
  # -----------------
  date_ref =  datetime(1850, 1, 1) # zero-time reference
  
  fileout = outdir + "siconc_SIday_NSIDC-0081_r1i1p1_" + dateStart + "-" + dateEnd + "_" + hemi + ".nc"
  # Create file
  f = Dataset(fileout, mode = "w")
  # Create dimensions
  time = f.createDimension('time', None)
  y   = f.createDimension('y', ny)
  x   = f.createDimension('x', nx)
  # Create variables
  times = f.createVariable('time', np.int32, ('time',))
  times[:] = np.arange((d1 - date_ref).days, (d2 - date_ref).days + 1) # + 1 because of Python indexing
  times.units = "days since " + str(date_ref.year) + "-" + str(date_ref.month) + "-" + str(date_ref.day)
  
  s = f.createVariable('siconc', np.float32, ('time', 'y', 'x'))
  s.units = "%"
  s[:] = siconc
  
  lons = f.createVariable('longitude', np.float32, ('y', 'x'))
  lons.units = "degrees East"
  lons[:] = lon
  
  lats = f.createVariable('latitude', np.float32, ('y', 'x'))
  lats.units = "degrees North"
  lats[:] = lat
  
  ca = f.createVariable('areacello', np.float32, ('y', 'x'))
  ca.units = "m2"
  ca[:] = cellarea
  
  sf  = f.createVariable('sftof', np.float32, ('y', 'x'))
  sf.units = "%"
  sf[:] = mask
  
  # Close
  f.close()
  print("")
  print("File produced: " + fileout)



if __name__ == "__main__":
  import sys

  # If not two arguments are passed
  if len(sys.argv) != (2 + 1):
    sys.exit("format_NSIDC-0081.py: No argument given, two expected.")
  else:
    dateStart = sys.argv[1]
    dateEnd   = sys.argv[2]

    formatData(dateStart, dateEnd)

