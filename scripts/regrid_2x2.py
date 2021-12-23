# Author: F. Massonnet
# Date  : December 2017
#         December 2018: update to account for new data arrival
#                        each year
#         February 2019: first version

#         April    2020: update 

# Purpose: Regridding of SIPN South data to common grid

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import csv
import numpy as np
import matplotlib.pyplot as plt
import sys

from   mpl_toolkits.basemap import Basemap
from   netCDF4 import Dataset
from   glob import glob


# label with the year investigated (2017-2018, 2018-2019, ...)

myyear = "2021-2022"

# Load namelist
exec(open("./namelist_spatial_" + myyear + ".py").read())

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# List of forecasts for each submission
list_for = [info[j_sub][1] for j_sub in range(n_sub)]
n_for    = [len(l) for l in list_for] # Nb of forecasts


# Tweak and add obs to list
obslist = ["OSI-401-b", "NSIDC-0081"]
n_sub += len(obslist)
sub_id = obslist + sub_id
list_for = [[0]] + [[0]] + list_for
n_for = [1, 1] + n_for

# Regridding at finer than 2x2 not possible due to native
# resolution of Lamont

# Target grid specification
# -------------------------
# Earth's radius
R = 6378000.0

# Longitude and latitude steps
dlon, dlat = 2.0, 2.0

lat_out = np.arange(-90.0 +  dlat / 2.0, 90.0  , dlat)
lon_out = np.arange(-180.0 + dlon / 2.0, 180.0 , dlon)

lon_out, lat_out = np.meshgrid(lon_out, lat_out)
ny, nx = lon_out.shape

# Grid cell area of output grid
areacello = R * np.cos( (lat_out * 2.0 * np.pi / 360.0 )) * \
                (dlon * 2.0 * np.pi / 360.0)  * R * \
                (dlat * 2.0 * np.pi / 360.0)

# First we create a land-sea mask for that grid, based on NSIDC's mask

print("Creating land-sea mask from NSIDC product")
f = Dataset("../data/" + myyear + "/netcdf/NSIDC-0081_000_concentration.nc", \
            mode = "r")
mask_in = f.variables["sftof"][:]
lat_mask= f.variables["latitude"][:]
lon_mask= f.variables["longitude"][:]
f.close()

lon_mask[lon_mask > 180.0] = lon_mask[lon_mask > 180.0] - 360.0
lon_mask[lon_mask < -180.0]= lon_mask[lon_mask < - 180.0] + 360.0
mask_out = np.full((ny, nx), np.nan)

# Loop over grid cells of the reference grid and define mask
for jy in np.arange(ny):
  for jx in np.arange(nx):
    
    # Locate all points that fall in that box
    lonmin = lon_out[jy, jx] - dlon / 2.0
    lonmax = lon_out[jy, jx] + dlon / 2.0
    latmin = lat_out[jy, jx] - dlat / 2.0
    latmax = lat_out[jy, jx] + dlat / 2.0

    # Identify the region of the input mask 
    region = (lat_mask >= latmin) * (lat_mask < latmax) * \
             (lon_mask >= lonmin) * (lon_mask < lonmax)
    # Set to land if at least one point of the native grid that is contained in the target grid is land
    if np.sum(region) == 0:
      mask_out[jy, jx] = 0.0
    elif np.min(mask_in[region]) == 0:
      mask_out[jy, jx] = 0.0
    else:
      mask_out[jy, jx] = 100.0


# Interpolation of the SIC files

for j_sub in range(n_sub):
  print("Interpolating " + str(sub_id[j_sub]))
  for j_for in list_for[j_sub]:
    # Open native file
    fileroot = "../data/" + myyear + "/netcdf/"
    filename = sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_concentration.nc"
    file = fileroot + filename

    f = Dataset(file, mode = "r")
    lon_in = f.variables["longitude"][:]
    lat_in = f.variables["latitude"][:]
    sic_in = f.variables["siconc"][:]
    f.close()

    # Recenter longitude to [-180, 180]
    lon_in[lon_in > 180.0]  = lon_in[lon_in > 180.0]  - 360.0
    lon_in[lon_in < -180.0] = lon_in[lon_in < -180.0] + 360.0

    # Expand to 2-D mesh grid if not done yet
    if len(lat_in.shape) == 1 or len(lat_in.shape) == 1:
      lon_in, lat_in = np.meshgrid(lon_in, lat_in)

    # If first member, create mask for that submission (needed only once)
    if j_for == list_for[j_sub][0]:
        
      # The idea here is to go through each grid cell of the target grid.
      # If that grid cell is "ocean" in the target grid, then an array
      # is stored for that cell that corresponds to the grid cells of the
      # *input* (submission) grid that are within the current target grid
      # cell
      list_mask = [list() for jxy in range(ny * nx)]
      for jy in np.arange(ny):
        for jx in np.arange(nx):
          if mask_out[jy, jx] == 100.0:
            # linear index
            jxy = jy * nx + jx

            lonmin = lon_out[jy, jx] - dlon / 2.0
            lonmax = lon_out[jy, jx] + dlon / 2.0
            latmin = lat_out[jy, jx] - dlat / 2.0
            latmax = lat_out[jy, jx] + dlat / 2.0

            mask   = (lat_in >= latmin) * (lat_in < latmax) * \
                     (lon_in >= lonmin) * (lon_in < lonmax)

            list_mask[jxy] = mask
            del mask

    print("Regridding " + filename)
  
    nt_in, ny_in, nx_in = sic_in.shape
    nt = nt_in
  
    # Output field
    sic_out = np.full((nt, ny, nx), np.nan)
    for jt in range(nt):
      print(str(jt))
      for jy in np.arange(ny):
        for jx in np.arange(nx):
          if mask_out[jy, jx] == 100.0:
            jxy = jy * nx + jx
            mask = list_mask[jxy]

            if np.sum(mask) == 0:
              sic_out[jt, jy, jx] = 0.0
            else:
              sic_out[jt, jy, jx] = np.mean(sic_in[jt, mask])       

    # Create file
    fileoutname = filename[:-3] + "_2x2.nc"
    fileout = fileroot + "/regrid/" + fileoutname
  
    print("Writing " + fileout)
    f = Dataset(fileout, mode = "w")
    # Create dimensions
    time = f.createDimension('time', nt)
    y =  f.createDimension('y', ny)
    x =  f.createDimension('x', nx)
  
    # Create variables
    times = f.createVariable('time', np.int32, ('time',))
    times.units = "months since 2018" + "-" + "12" + "-" + "01"
  
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
  
    sftof = f.createVariable("sftof", np.int32, ("y", "x"))
    sftof[:] = mask_out
  
    # Close
    f.close()
    print("Created " + fileout)
