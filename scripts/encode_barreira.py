#!/usr/bin/python

import csv
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Encore barreira data in the folder raw of that year

myyear = "2023-2024"
nameContrib = "barreira"

R = 6378000.0

dlon, dlat = 1.0, 1.0

lat_out = np.arange(-90.0 + dlat / 2.0, 90.0   , dlat)
lon_out = np.arange(-180.0 + dlon / 2.0, 180.0 , dlon)

lon_out, lat_out = np.meshgrid(lon_out, lat_out)
ny, nx = lon_out.shape

areacello = R * np.cos( (lat_out * 2.0 * np.pi / 360.0 )) * (dlon * 2.0 * np.pi / 360.0)  * R * (dlat * 2.0 * np.pi / 360.0)

root = "../data/" + myyear + "/raw/" + nameContrib + "/"
nmemb = 1

for jmemb in np.arange(1, nmemb + 1):
  jjmemb = str(jmemb).zfill(3)

  #files = ["./memb" + jjmemb + "/SHN_" + str(j) + "_" + str(jmemb) + ".txt" for j in range(1, 90 + 1)]

  files = list()
  for k in [["dic", 31, myyear[:4]], ["Ene", 31, myyear[5:]], ["Feb", 28, myyear[5:]]]:
  #  thisList = ["./memb" + jjmemb + "/SHN_" + k[0] + "_" + str(dayNam) + "_" + str(jmemb) + "_original.txt" for dayNam in np.arange(1, k[1] + 1)]
    thisList = ["SHN_" + k[0] + "_" + k[2][2:] + " " + str(dayNam) + " " + "concentracion.txt" for dayNam in np.arange(1, k[1] + 1)]
    files += thisList

  print(files)
  nt = len(files)


  for jt in range(nt):
      print(jt)

      data = list()
    
      if os.path.exists(root + files[jt]):

          with open (root + files[jt], "r") as f:
              csvreader = csv.reader(f)
              next(csvreader)
              for row in csvreader:
                  data.append(np.array([float(r) for r in row[0].split()]))
        
          data = np.array(data)
        
          lat = data[:, 0]
          lon = data[:, 1] 
          are = data[:, 2]
          sic = data[:, 3]
        
          # Recenter longitude to [-180, 180]
          lon[lon > 180.0] = lon[lon > 180.0] - 360.0
        
          # Initialize array if first time step
          if jt == 0:
              sic_out = np.full((nt, ny, nx), np.nan)
        
          for jy in np.arange(ny):
              for jx in np.arange(nx):
                  # Locate all points that fall in that box
                  lonmin = lon_out[jy, jx] - dlon / 2.0
                  lonmax = lon_out[jy, jx] + dlon / 2.0
                  latmin = lat_out[jy, jx] - dlat / 2.0
                  latmax = lat_out[jy, jx] + dlat / 2.0
          
                  mask = (lat >= latmin) * (lat < latmax) * (lon >= lonmin) * (lon < lonmax)
          
                  if np.sum(mask) == 0:
                      sic_out[jt, jy, jx] = 0.0
                  else:
                      sic_out[jt, jy, jx] = np.mean(sic[mask])
      else:
          print("File " + root + files[jt] + " not found, skipping")
    
  # Create file
  fileout = "../data/" + myyear + "/netcdf/" + nameContrib + "_" + jjmemb + "_" + myyear[:4] + "1201-" + str(myyear[5:]) + "0228" + "_concentration.nc"
  f = Dataset(fileout, mode = "w")

  # Create dimensions
  time = f.createDimension('time', nt)
  y =  f.createDimension('y', ny)
  x =  f.createDimension('x', nx)
  # Create variables
  times = f.createVariable('time', np.int32, ('time',))
  times.units = "days since 2018 " + myyear[:4]  + "-" + "12" + "-" + "01"

  h = f.createVariable("siconc", np.float32, ('time', 'y', 'x'))
  h.units = "%"
  h[:] = sic_out

  lons = f.createVariable("longitude", np.float32, ("y", "x"))
  lons[:] = lon_out

  lats = f.createVariable("latitude", np.float32, ("y", "x"))
  lats[:] = lat_out

  cellarea = f.createVariable("areacello", np.float32, ("y", "x"))
  cellarea.units = "m2"
  cellarea[:] = areacello

  sftof = f.createVariable("sftof", np.float32, ("y", "x"))
  sftof.units = "%"
  sftof[:] = np.full((ny, nx), 100.0)

  # Close
  f.close()
