# Author: F. Massonnet
# Date  : 10 oct 2022


# Purpose: Figure with daily integrated ice edge error

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import pandas            as pd
import matplotlib; matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
import matplotlib.colors
import numpy             as np
import os
import random

from datetime import datetime
from datetime import timedelta
from netCDF4  import Dataset
from datetime import date
from scipy    import stats

from seaice_commondiags import *

# Functions
# ---------

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

matplotlib.rcParams['font.family'] = "Arial Narrow"

plt.close("all")


# Open namelist information file
# ------------------------------
exec(open("./namelist.py").read())

# Script parameters
# -----------------

seasonId = 4 # Which season to look at
diagId   = 2 # We look at sea ice concentration

# ----
nDaysWeek = 7 # how many days in a week

nDays    = ((namelistOutlooks[seasonId][5]) - (namelistOutlooks[seasonId][4])).days + 1

daysAxis = [namelistOutlooks[seasonId][4] + timedelta(days = float(d)) for d in np.arange(nDays)]


fig, ax = plt.subplots(1, 1, figsize = (6, 4))


seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)

# Load observational data first
# -----------------------------
obs_name = "NSIDC-0081"
f = Dataset("../data/" + seasonName + "/netcdf/regrid/" + \
            obs_name + "_000_concentration_2x2.nc")
sic_obs = f.variables["siconc"][:]
latitude = f.variables["latitude"][:]
longitude = f.variables["longitude"][:]
cellarea  = f.variables["areacello"][:]
mask_obs  = f.variables["sftof"][:]
nt = f.dimensions["time"].size
f.close()

# Run through all forecasts

for j, n in enumerate(namelistContributions):
	print(n)
	thisName = n[0]
	thisNbForecasts	= n[seasonId + 1][diagId]
	thisType = n[-1] # Statistical or dynamical

	if thisNbForecasts > 0:
		thisIIEE  = list()
		thisNIIEE = list()
		thisAEE   = list()
		thisME    = list()

		for jFor in np.arange(thisNbForecasts):
			filein = "../data/" + seasonName + "/netcdf/regrid/" + thisName + "_" \
				  + str(jFor + 1).zfill(3) + "_concentration_2x2.nc"
	
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
		
		# Plot series
		if thisType == "s":
			thisColor = plt.cm.YlOrRd( int(128 + j / len(namelistContributions) * 128))
		elif thisType == "d":
			thisColor = plt.cm.PuBuGn( int(128 + j / len(namelistContributions) * 128))
		else:
			thisColor = "black"

		ax.plot(daysAxis, IIEE, color = thisColor, label = thisName)
		mymaxIIEE = np.max(np.array(thisIIEE), axis = 0)
		myminIIEE = np.min(np.array(thisIIEE), axis = 0)
		ax.fill_between(daysAxis, myminIIEE, mymaxIIEE, color = thisColor, \
                alpha = 0.2, lw = 0)

		# Add text to locate contributions
		ax.text(daysAxis[0] - timedelta(days  = 1), IIEE[0], thisName[0], ha = "right", fontsize = 4, color = [0.5, 0.5, 0.5], va = "center")

ax.legend()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
ax.set_ylabel("Million km$^2$")
ax.set_title("Integrated Ice Edge Error (" + seasonName + ")")
fig.savefig("../figs/fig6_paper.png", dpi = 300)
