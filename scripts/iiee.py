#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : April 2019
# Data: SIPN South contributors

# Integrated ice edge errror (Goessling et al., GRL, 2016)

# Requires all data to be on the same grid

# Imports
import pandas as pd
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import random
import matplotlib.dates as mdates

#import cartopy.crs as ccrs ; import cartopy.feature as cfeature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from   netCDF4 import Dataset
from   datetime import datetime


# Script parameters
myyear = "2021-2022"  # label with the year investigated (2017-2018, 2018-2019, ...)

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
                        # (Pythonic convention)
else:
  # Initialization date
  inidate = myyear[:4] + "1201"
  # Number of days in the forecast period
  ndays   = 90
  # Label for period that is forecasted
  period_name = "Dec-Jan-Feb " + myyear[:4] + "-" + myyear[5:]
  # Starting and ending time indices (Python conventions)
  t1, t2 = 63 - 1, 63 - 1 + 28
  target_period_name = "February"

# Time axis
time = pd.date_range(pd.to_datetime(inidate, format = "%Y%m%d"), 
                     periods = ndays).tolist()

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

def iiee(sic_eva, sic_ref, cellarea, mask = 1, threshold = 15.0, lat = None, 
         lon = None, plot = False):
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
  
  overestim = np.array([np.sum(1.0 * (sic_eva[jt, :, :] >  threshold) * \
                  (sic_ref[jt, :, :] <= threshold)  * mask * cellarea) \
                  for jt in range(nt)]) / 1e12
  underestim= np.array([np.sum(1.0 * (sic_eva[jt, :, :] <= threshold) * \
                  (sic_ref[jt, :, :] > threshold)   * mask * cellarea) \
                  for jt in range(nt)]) / 1e12

  ref_area  = np.array([np.nansum(sic_ref[jt, :, :] / 100.0 * mask * cellarea)\
                        for jt in range(nt)]) / 1e12
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
    map = Basemap(projection = "spstere", boundinglat = - 50, 
                  lon_0 = 180, resolution = 'l')
    x, y = map(lon, lat)
    plt.subplot(1, 1, 1)
    cs = map.contourf(x, y, sic_ref[0, :, :], clevs, cmap = plt.cm.Blues_r, 
                      latlon = False, extend = "neither")
    map.fillcontinents(color = 'grey', lake_color = 'w')
    map.drawcoastlines(linewidth = 1.0)
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

f = Dataset("../data/" + myyear + "/netcdf/regrid/" + \
            obs_name + "_000_concentration_2x2.nc")
sic_obs = f.variables["siconc"][:]
latitude = f.variables["latitude"][:]
longitude = f.variables["longitude"][:]
cellarea  = f.variables["areacello"][:]
mask_obs  = f.variables["sftof"][:]
nt = f.dimensions["time"].size
f.close()

fig, ax = plt.subplots(4, 1, figsize = (8, 16))
# Load model data
for j_sub in range(n_sub):
  print("Doing " + str(sub_id[j_sub]))

  thisIIEE  = list()
  thisNIIEE = list()
  thisAEE   = list()
  thisME    = list()

  for j_for in list_for[j_sub]:
    filein = "../data/" + myyear + "/netcdf/regrid/" + sub_id[j_sub] + "_" \
    + str(j_for).zfill(3) + "_concentration_2x2.nc"

    if not os.path.exists(filein):
      print(filein + " not found")

    else:
      # Open file, read geometric parameters if the first one 
      print("  Loading " + filein)
      f = Dataset(filein, mode = "r")
      sic = f.variables["siconc"][:]
      f.close()
      IIEE, NIIEE, AEE, ME, O, U = iiee(sic, sic_obs, cellarea, mask = 1.0 * \
             (mask_obs == 100.0) * (latitude < 0), threshold = 15.0, 
             lat = latitude, lon = longitude, plot = False)

      thisIIEE.append(IIEE)
      thisNIIEE.append(NIIEE)
      thisAEE.append(AEE)
      thisME.append(ME)

      # Plot series, line thinner for large ensembles. 
      # Legend only if first member
  mylab = info[j_sub][0] + " " + info[j_sub][3]
  
  
  # IIEE
  meanIIEE = np.mean(np.array(thisIIEE), axis = 0)
  ax[0].plot(time, meanIIEE, color = col[j_sub], lw = 1.5, label = mylab)
  # Add submission ID
  ax[0].text(time[0], meanIIEE[0], str(j_sub), fontsize = 4)
  # Plot range as shading
  mymaxIIEE = np.max(np.array(thisIIEE), axis = 0)
  myminIIEE = np.min(np.array(thisIIEE), axis = 0)
  ax[0].fill_between(time, myminIIEE, mymaxIIEE, color = col[j_sub], \
                  alpha = 0.5, lw = 0)
  ax[0].set_title(period_name + " Integrated Ice Edge Error")
  ax[0].set_ylabel("10$^6$ km$^2$")
  
  # NIIEE
  meanNIIEE = np.mean(np.array(thisNIIEE), axis = 0)
  ax[1].plot(time, meanNIIEE, color = col[j_sub], lw = 1.5, label = mylab)
  # Add submission ID
  ax[1].text(time[0], meanNIIEE[0], str(j_sub), fontsize = 4)
  # Plot range as shading
  mymaxNIIEE = np.max(np.array(thisNIIEE), axis = 0)
  myminNIIEE = np.min(np.array(thisNIIEE), axis = 0)
  ax[1].fill_between(time, myminNIIEE, mymaxNIIEE, color = col[j_sub], \
                  alpha = 0.5, lw = 0)
  ax[1].set_title(period_name + " Normalized Integrated Ice Edge Error")
  ax[1].set_ylabel("%")
  
  # AEE
  meanAEE = np.mean(np.array(thisAEE), axis = 0)
  ax[2].plot(time, meanAEE, color = col[j_sub], lw = 1.5, label = mylab)
  # Add submission ID
  ax[2].text(time[0], meanAEE[0], str(j_sub), fontsize = 4)
  # Plot range as shading
  mymaxAEE = np.max(np.array(thisAEE), axis = 0)
  myminAEE = np.min(np.array(thisAEE), axis = 0)
  ax[2].fill_between(time, myminAEE, mymaxAEE, color = col[j_sub], \
                  alpha = 0.5, lw = 0)
  ax[2].set_title(period_name + " Absolute Extent Error")
  ax[2].set_ylabel("10$^6$ km$^2$")
  

  # ME
  meanME = np.mean(np.array(thisME), axis = 0)
  ax[3].plot(time, meanME, color = col[j_sub], lw = 1.5, label = mylab)
  # Add submission ID
  ax[3].text(time[0], meanME[0], str(j_sub), fontsize = 4)
  # Plot range as shading
  mymaxME = np.max(np.array(thisME), axis = 0)
  myminME = np.min(np.array(thisME), axis = 0)
  ax[3].fill_between(time, myminME, mymaxME, color = col[j_sub], \
                  alpha = 0.5, lw = 0)
  ax[3].set_title(period_name + "Misplacement Error")
  ax[3].set_ylabel("10$^6$ km$^2$")

for a in ax:
  a.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
  a.set_xticks([time[j] for j in [0, 14, 31, 45, 62, 76, 89]])
  a.grid()
  a.set_axisbelow(True)
  a.legend(loc=(1.04,0))

fig.tight_layout()

for fmt in ["png",]:# "eps", "pdf"]:
    fig.savefig("../figs/iiee." + fmt, dpi = 300)
    

