#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : May 2019
# Data: SIPN South contributors

# Maps of differences

# Requires all data to be on the same grid

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import random

from   mpl_toolkits.basemap import Basemap, addcyclic
from   netCDF4 import Dataset
from   datetime import datetime


# Script parameters
myyear = "2018-2019"  # label with the year investigated (2017-2018, 2018-2019, ...)

# Load namelist
exec(open("./namelist_spatial_" + myyear + ".py").read())

# Tweak and add OSISAF as last model (to check baseline error)
info = info + [["OSI-401-b", [0], [0.2, 0.2, 0.2], ""]]
# Submission ID
n_sub = len(info)

sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

n_sub = len(info)

col = [info[j_sub][2] for j_sub in range(n_sub)]


if myyear == "2017-2018":
  inidate = "20180201"
  ndays   = 28 #number of days of the forecast period
  period_name = "February 2018"
  t1, t2 = 0, 0 + ndays # time indices defining the period for diagnostics
                        # (Pythonic convention)

elif myyear == "2018-2019":
  inidate = "20181201"
  ndays   = 90
  period_name = "February 2019"
  #period_name = "1 December 2019"
  period_short_name = "Feb2019"
  #period_short_name = "1Dec2019"
  t1, t2 = 63 - 1, 63 - 1 + 28
  #t1, t2 = 1 - 1, 1 - 1 + 1

# Time axis
time = pd.date_range(pd.to_datetime(inidate, format = "%Y%m%d"), periods = ndays).tolist()

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# List of forecasts for each submission
list_for = [info[j_sub][1] for j_sub in range(n_sub)]
n_for    = [len(l) for l in list_for] # Nb of forecasts

# Observations
obs_name = "NSIDC-0081"

# ======================

# Load observational data
# -----------------------

f = Dataset("../data/" + myyear + "/netcdf/regrid/" + obs_name + "_000_concentration_2x2.nc")
sic_obs = f.variables["siconc"][:]
latitude = f.variables["latitude"][:]
longitude = f.variables["longitude"][:]
cellarea  = f.variables["areacello"][:]
mask_obs  = f.variables["sftof"][:]
nt = f.dimensions["time"].size
f.close()

# Create projection
# -----------------
map = Basemap(projection = "spstere", boundinglat = -50.0, \
                lon_0 = 180.0, resolution = 'l')

# Map the lon, lat to the projection
x, y = map(longitude, latitude)

plt.figure(figsize = (6, 6))
# Load model data
for j_sub in range(n_sub):
  print("Doing " + str(sub_id[j_sub]))

  submission = list()

  for j_for in list_for[j_sub]:
    filein = "../data/" + myyear + "/netcdf/regrid/" + sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_concentration_2x2.nc"

    if not os.path.exists(filein):
      sys.exit(filein + " not found")

    # Open file, read geometric parameters if the first one 
    print("  Loading " + filein)
    f = Dataset(filein, mode = "r")
    sic = f.variables["siconc"][:]
    f.close()

    submission.append(sic)


  mean = np.mean(np.array(submission), axis = 0)


  map.drawcoastlines(linewidth = 0.25)
  map.fillcontinents(color = 'grey', lake_color = 'w')
  map.drawmeridians(np.arange(0, 360, 30))
  map.drawparallels(np.arange(-90, 90, 10))


  cs = map.contourf(longitude, latitude, sic[0, :, :] - sic_obs[0, :, :], np.arange(-33.0, 36.0, 6.0), cmap = plt.cm.RdBu_r, \
                    latlon = True, extend = "both")
  map.contour(longitude, latitude, sic_obs[0, :, :], [15.0], colors = "k", linewidths = 2)
  map.contour(longitude, latitude, sic[0, :, :]    , [15.0], colors = "darkgreen", linewidths = 2)
  cbar = map.colorbar(cs, location = 'bottom', pad = "5%", ticks = [-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0])

  cbar.set_label("%")
  plt.title(sub_id[j_sub])

  plt.savefig("../figs/map_diff_" + sub_id[j_sub] + ".png", dpi = 500)
