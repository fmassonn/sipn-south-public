#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : April 2019
# Data: SIPN South contributors

# Integrated ice edge errror (Goessling et al., GRL, 2016)

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

def iiee(sic_eva, sic_ref, cellarea, mask = 1, threshold = 15.0, lat = None, lon = None, plot = False):
  """
  sic_ref -> reference field [%]
  sic_eva -> evaluated field    [%]
  cellarea -> grid cell area [m2]
  mask -> 1 where data has to be included
  threshold -> sea ice edge imit

  lat, lon, plot: to plot what's going on

  returns iiee, that is, sum of grid cell
    areas where sic1 and sic2 disagree on the
    event "> 15%" or " < 15%" (units million km2)
  """

  if sic_ref.shape != sic_eva.shape:
    sys.exit("iiee: sic_ref and sic_eva have different shapes")

  nt, _, _ = sic_ref.shape
  
  overestim = np.array([np.sum(1.0 * (sic_eva[jt, :, :] >  threshold) * (sic_ref[jt, :, :] <= threshold)  * mask * cellarea) for jt in range(nt)]) / 1e12
  underestim= np.array([np.sum(1.0 * (sic_eva[jt, :, :] <= threshold) * (sic_ref[jt, :, :] > threshold)   * mask * cellarea) for jt in range(nt)]) / 1e12

  ref_area  = np.array([np.nansum(sic_ref[jt, :, :] / 100.0 * mask * cellarea) for jt in range(nt)]) / 1e12
  AEE = np.abs(overestim - underestim)
  ME  = 2.0 * np.minimum(overestim, underestim)

  IIEE = overestim + underestim
  NIIEE = 100.0 * IIEE /ref_area

  if np.max(np.abs((AEE + ME) - IIEE)) > 1e-14:
    print(np.abs((AEE + ME) - IIEE))
    sys.exit("ERROR")


  if plot:
    clevs = np.arange(0.0, 110.0, 10.0)
    from   mpl_toolkits.basemap import Basemap, addcyclic
    plt.figure(figsize = (6, 6))
    map = Basemap(projection = "spstere", boundinglat = - 50,lon_0 = 180, resolution = 'l')
    x, y = map(lon, lat)
    plt.subplot(1, 1, 1)
    cs = map.contourf(x, y, sic_ref[0, :, :], clevs, cmap = plt.cm.Blues_r, latlon = False, extend = "neither")
    map.fillcontinents(color = 'grey', lake_color = 'w'); map.drawcoastlines(linewidth = 1.0)
    map.drawmeridians(np.arange(0, 360, 30), color = [0.7, 0.7, 0.7])
    map.drawparallels(np.arange(-90, 90, 10), color = [0.7, 0.7, 0.7])
    cbar = map.colorbar(cs, location = 'bottom', pad = "5%")
    cbar.set_label("%")
    plt.title("REF")
    plt.savefig("../figs/map.png")
    plt.close("fig")

  
  return IIEE, NIIEE, AEE, ME, overestim, underestim

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

plt.figure(figsize = (6, 4))
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
    IIEE, NIIEE, AEE, ME, O, U = iiee(sic, sic_obs, cellarea, mask = 1.0 * (mask_obs == 100.0) * (latitude < 0), threshold = 15.0, lat = latitude, lon = longitude, plot = False)

    submission.append(IIEE)

    # Plot series, line thinner for large ensembles. Legend only if first member
    if j_for == list_for[j_sub][0]:
      mylab = info[j_sub][0] + " " + info[j_sub][3]
    else:
      mylab = "_nolegend_"

  mean = np.mean(np.array(submission), axis = 0)
  plt.plot(time, mean, color = col[j_sub], lw = 1.5, label = info[j_sub][0] + " " + info[j_sub][3])
  # Plot range as shading
  mymax = np.max(np.array(submission), axis = 0)
  mymin = np.min(np.array(submission), axis = 0)
  plt.fill_between(time, mymin, mymax, color = [1.0 * c for c in col[j_sub]], alpha = 0.5, lw = 0)

plt.title(period_name + " Integrated Ice Edge Error")
plt.ylabel("10$^6$ km$^2$")
import matplotlib.dates as mdates
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.grid()
#plt.legend()
plt.savefig("../figs/iiee.png", dpi = 300)
