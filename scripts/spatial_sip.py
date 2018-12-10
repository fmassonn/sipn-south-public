#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : December 2017
# Data: SIPN South contributors

# Imports
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt#; import matplotlib; matplotlib.use("Agg")
import numpy as np
import os
from netCDF4 import Dataset
import sys
from mpl_toolkits.basemap import Basemap, addcyclic


print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")
print("WARNING: NRL between 0 and 1 for incoming forecasts")

# Load color maps

myyear = "2018-2019"  # label with the year investigated (2017-2018, 2018-2019, ...)
plotobs = False       # Add obs as reference or not (False if forecast mode)

# Load namelist
exec(open("./namelist_spatial_" + myyear + ".py").read())

# Plot setup
clevs = np.arange(0, 105, 5)


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
  period_short_name = "Feb2019"
  t1, t2 = 63 - 1, 63 - 1 + 28

# Time axis
time = pd.date_range(pd.to_datetime(inidate, format = "%Y%m%d"), periods = ndays).tolist()

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# Number of forecasts
n_for = [info[j_sub][1] for j_sub in range(n_sub)]

# Observations
obs_labels = ["OSI-401-b", "NSIDC-0081"]
n_obs = len(obs_labels)

# ======================




# Load observational data
# -----------------------

if plotobs:
  obs = list() # Each item of the list will have one x_array, y_array, SIC for each day
  for j_obs in range(n_obs):
    f = Dataset("../data/" + myyear + "/netcdf/" + obs_labels[j_obs] + "_000_concentration.nc")
    sic = f.variables["siconc"][t1:t2, :, :]
    latitude = f.variables["latitude"][:]
    longitude = f.variables["longitude"][:]
    map = Basemap(projection = "spstere", boundinglat = - 50,lon_0 = 180, resolution = 'l')
    x, y = map(longitude, latitude)
    obs.append((x, y, sic))
    del x, y, sic
    f.close()


for j_sub in range(n_sub):
  print("Doing " + str(sub_id[j_sub]))
  map = Basemap(projection = "spstere", boundinglat = - 50,lon_0 = 180, resolution = 'l')
  map.fillcontinents(color = 'grey', lake_color = 'w')

  # Load the data into a big array "data"
  for j_for in np.arange(1, n_for[j_sub] + 1):
    filein = "../data/" + myyear + "/netcdf/" + sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_concentration.nc"

    if not os.path.exists(filein):
      sys.exit(filein + " not found")

    # Open file, read geometric parameters if the first one 
    print("  Loading " + filein)
    f = Dataset(filein, mode = "r")
    if j_for == 1:
      lat = f.variables["latitude"][:]
      lon = f.variables["longitude"][:]
      nt = f.dimensions["time"].size
      if len(lon.shape) == 1:
        lon, lat = np.meshgrid(lon, lat)

      x, y = map(lon, lat)
      ny, nx = lat.shape
      data = np.empty((n_for[j_sub], nt, ny, nx))

    data[j_for - 1, :, :, :] = f.variables["siconc"][:]
    sftof = f.variables["sftof"][:]
    f.close()


  # Regularize
  data[data > 100.0] = 100.0
  data[data < 0.0  ] = 0.0
  # Do some plots

  # Monthly mean for each member
  for j_for in np.arange(1, n_for[j_sub] + 1):
    fig = plt.figure("fig", figsize = (5, 5))
    sic_monmean = np.mean(data[j_for - 1, t1:t2, :, :], axis = 0)
    cl = map.contour(x, y, sic_monmean, [15.0], latlon = False, colors = '#ffcccc', linewidths = 1, linestyles = "-")
    cs = map.contourf(x, y, sic_monmean, clevs, cmap = plt.cm.PuBu_r, latlon = False, extend = "neither")     
    map.fillcontinents(color = 'grey', lake_color = 'w'); map.drawcoastlines(linewidth = 1.0)
    map.drawmeridians(np.arange(0, 360, 30), color = [0.7, 0.7, 0.7])
    map.drawparallels(np.arange(-90, 90, 10), color = [0.7, 0.7, 0.7])
    cbar = map.colorbar(cs, location = 'bottom', pad = "5%")
    cbar.set_label("%")
    plt.title(sub_id[j_sub] + " | member " + str(j_for).zfill(3)  + " | " + period_name + " mean")
    plt.savefig("../figs/" + sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_concentration_" + period_short_name + "-mean")
    print("    Monthly mean conc printed for " + sub_id[j_sub] + " " + str(j_for).zfill(3))
    plt.close("fig")


  # Monthly mean for the ensemble mean + spaghetti for forecasts
  fig = plt.figure("fig", figsize = (5, 5))
  sic_grandmean = np.mean(np.mean(data, axis = 0)[t1:t2], axis = 0)
  for j_for in np.arange(1, n_for[j_sub] + 1):
    print("      Contour for forecast # " + str(j_for).zfill(3))
    cl = map.contour(x, y, np.mean(data[j_for - 1, t1:t2, :, :], axis = 0), [15.0], latlon = False, colors = '#ffcccc', linewidths = 0.1, linestyles = "-")

  if sub_id[j_sub] == "ucl":
    cs = map.contourf(x[0:-1,0:-1], y[0:-1,0:-1], sic_grandmean[0:-1,0:-1], clevs, cmap = plt.cm.PuBu_r, latlon = False, extend = "neither")     
  else:
    cs = map.contourf(x, y, sic_grandmean, clevs, cmap = plt.cm.PuBu_r, latlon = False, extend = "neither")   

  # Put obs
  if plotobs:
    for j_obs in range(n_obs):
      map.contour(obs[j_obs][0], obs[j_obs][1], np.mean(obs[j_obs][2], axis = 0), [15.0], latlon = False, colors = 'y', linewidths = 1, linestyles = "-")
    
  map.fillcontinents(color = 'grey', lake_color = 'w'); map.drawcoastlines(linewidth = 1.0)
  map.drawmeridians(np.arange(0, 360, 30), color = [0.7, 0.7, 0.7])
  map.drawparallels(np.arange(-90, 90, 10), color = [0.7, 0.7, 0.7])
  cbar = map.colorbar(cs, location = 'bottom', pad = "5%")
  cbar.set_label("%")
  plt.title(sub_id[j_sub] + " | ens mean | " + period_name + " mean")
  plt.savefig("../figs/" + sub_id[j_sub] + "_ens-mean" + "_concentration_" + period_short_name + "-mean.png", dpi = 200)        
  print("  Monthly mean conc printed for ensemble mean")
  plt.close("fig")


  # Probability of sea ice concentration being more than 15% (for each day)
  for jt in np.arange(t1, t2):
    fig = plt.figure("fig", figsize = (5, 5))
    prob = 100.0 * np.sum(1.0 * (data[:, jt, :, :] > 15.0), axis = 0) / n_for[j_sub]
    if sub_id[j_sub] == "ucl":
      cs = map.contourf(x[0:-1,0:-1], y[0:-1,0:-1], prob[0:-1,0:-1], clevs, cmap = plt.cm.RdYlGn_r, latlon = False, extend = "neither")
    else:
      cs = map.contourf(x, y, prob, clevs, cmap = plt.cm.RdYlGn_r, latlon = False, extend = "neither")
    map.fillcontinents(color = 'grey', lake_color = 'w'); map.drawcoastlines(linewidth = 1.0)
    map.drawmeridians(np.arange(0, 360, 30), color = [0.7, 0.7, 0.7])
    map.drawparallels(np.arange(-90, 90, 10), color = [0.7, 0.7, 0.7])
    cbar = map.colorbar(cs, location = 'bottom', pad = "5%")
    cbar.set_label("%")

    # Plot contours from two observations
    if plotobs:
      for j_obs in range(n_obs):
        map.contour(obs[j_obs][0], obs[j_obs][1], obs[j_obs][2][jt, :, :], [15.0], latlon = False, colors = 'w', linewidths = 1, linestyles = "-")

    # Print line on the right
    plt.title(sub_id[j_sub] + " | prob > 15% | " + time[jt].strftime('%d %B %Y'))


    #x1 = map(-135, -40)[0]
    #y1 = map(-135, -40)[1]
    #x2 = map( 135, -40)[0]
    #y2 = map( 135, -40)[1]

    #plt.plot((x1, x2), (y1, y2), color = [0.7, 0.7, 0.7], lw = 2, zorder = 0)
    #plt.scatter(np.linspace(x1, x2, nt),     np.linspace(y1, y2, nt), 50, color = [0.7, 0.7, 0.7])
    #plt.scatter(np.linspace(x1, x2, nt)[jt], np.linspace(y1, y2, nt)[jt], 80, color = [0.1, 0.1, 0.1])
    #for day in [1, 10, 20, 28]:
    #  plt.text(np.linspace(x1, x2, nt)[day - 1], np.linspace(y1, y2, nt)[day - 1], " " + period_name + " " + str(day).zfill(2), rotation = 90, va = "bottom", ha = "center", color = [0.7, 0.7, 0.7])
    plt.savefig("../figs/" + sub_id[j_sub] + "_prob-15" + "_concentration_" + "d" + str(jt + 1).zfill(2) + ".png", dpi = 200)
    print("  Probability day " + str(jt + 1).zfill(2) + " printed")
    plt.close("fig")
  
