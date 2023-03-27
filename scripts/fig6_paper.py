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
from matplotlib.lines import Line2D
import numpy             as np
import os
import random

from datetime import datetime
from datetime import timedelta
from netCDF4  import Dataset
from datetime import date
from scipy    import stats

from seaice_commondiags import *

thresholdSeaIcePresence = 15.0 # When a grid cell is said to have ice

# Note that there is some trickiness in the computation of the IIEE, to be well understood.
#
# The core idea of the IIEE is to compare two *contours* of sea ice defined by sea ice presence, defined as SIC > 15%

# How do we do that?

# 1/ We binarize the verifying SIC, i.e., we turn it to 0 and 1 based on some accepted threshold of sea ice presence, 15%

# 2/ We binarize the forecast SIC the same way. We do that for each ensemble member separately

# 3/ We average (in case of ensemble members) the fields in 2/ to obtain a probability of ice presence, SIP

# 4/ We compute the IIEE as the sum of grid cells where either the obs says there ice (> 15%) but the forecast says there is not (SIP < 50%), or the obs says there is no ice (<15%) and the forecast says there is (SIP > 50%)

# Note that if only one ensemble member is available, SIC > 15% <-> SIP=1 > 50% and SIC < 15% <-> SIP=0 < 50% so it's equivalent

# An example:

#			cell 1		cell2		cell3

# SIC fields
# Obs:     		10%		20%		100%

# Forecast member 1:	30%		70%		90%
# Forecast member 2:	50%		80%		20%
# Forecast member 3:	60%		0%		10%
# -----------------

# Binarized/sea ice presence fields:
# Obs:			0		1		1

# Forecast member 1:	1		1		1
# Forecast member 2:	1		1		1
# Forecast member 3:	1		0		0

# Forecast average:	1		0.66		0.66
# -----------------

# In all grid cells the forecast predicts more likely than not that there will be ice (>0.5 probability)
# In obs there is ice in 2 out of 3 grid cells
# --> The IIEE is the area of cell 1, where a disagreement exists.

# Note that if the forecast was only member 1, the IIEE would be also cell1 . If the forecast was only one member (member 3), the IIEE would be cell1 + cell2 + cell3


# --> These numbers may be different, see correspondence with H. Goessling (e-mail 10 Jan 2023):
# The following example should show that these seemingly inconsistent thresholds are actually consistent: Imagine the ensemble size is 1. In that case, the 50% SIP contour will be identical with the 15% SIC contour of the single ensemble member, so it obviously makes sense to compare that against the observed 15%-SIC contour. And it remains consistent for actual ensembles with more members. The 50%-SIP contour is simply the „median“ of all individual 15%-SIC contours

# Functions
# ---------

def iiee(sipForecast, sicVerif, cellarea, mask = 1, thresholdProbability = 0.5, thresholdSeaIcePresence = 15.0, timeDim = True):
  """
  Arguments:
    - sipForecast  -> probability of sea ice presence in forecast field (binary if one ensemble member, fractional else)
                      numpy array of values between 0.0 and 1.0
    - sicVerif     -> sea ice concentration in verification field, in %
                      numpy array of values between 0.0 and 100.0, same shape as sipForecast

    - cellarea     -> areas of grid cells in m2

    - mask         -> mask (1 or 0) to account e.g. for land and ocean distribution

    - thresholdProbability -> a value between 0 and 1 to define above what threshold sea ice is considered to be present in the forecast. Typically 0.5 ("more likely than not").

    - thresholdSeaIcePresence -> for obervations or one individual physical forecast, the threshold above which there is enough ice to declare ice presence. Typically 15%

    - timeDim -> Boolean to state if the first dimension of the arrays is a time dimension
   

  Returns: IIEE in million km2

  """

  if sipForecast.shape != sicVerif.shape:
    sys.exit("iiee: sipForecast and sicVerif have different shapes")

  if timeDim:
    nt, _, _ = sicVerif.shape
 
    IIEE = list()

    # Create a new array with nans where there are nans in obs, 1 where there is ice, and 0 where there is not
    presenceVerif = np.full(sicVerif.shape, np.nan)
    presenceVerif[sicVerif >= thresholdSeaIcePresence] = 1
    presenceVerif[sicVerif <  thresholdSeaIcePresence] = 0


    # Create a new arrray with nans where there are nans in forecast, 1 where there is ice and 0 else
    presenceForecast = np.full(sipForecast.shape, np.nan)
    presenceForecast[sipForecast >= thresholdProbability] = 1
    presenceForecast[sipForecast <  thresholdProbability] = 0
     
    # Sum these fields. Only where the sum is equal to 1, we have a disagreement
    disagree = ((presenceForecast + presenceVerif) == 1)

    for jt in np.arange(nt):
      if np.sum(np.isnan(presenceVerif[jt, :, :])) == presenceVerif[jt, :, :].size:
        total = np.nan
      else:
        total = np.sum(disagree[jt, :, :] * cellarea * mask) / 1e12
      
      IIEE.append(total)
  else:
    print("not coded yet")
    stop

  return IIEE

matplotlib.rcParams['font.family'] = "Arial Narrow"

plt.close("all")


# Open namelist information file
# ------------------------------
exec(open("./namelist.py").read())

# Script parameters
# -----------------

seasonId = 4 # Which season to look at
diagId   = 2 # We look at sea ice concentration

colorDict = {"statistical": "#008579", \
             "dynamical":   "#FF7200", \
             "climatology": "#000000", \
             "group forecast": "#0000FF", \
             "alternative verification": "#808080" }

plotOtherVerif = True # Whether to plot alternative dataset  IIEE to gauge obs uncertainty

# ----
nDays    = ((namelistOutlooks[seasonId][5]) - (namelistOutlooks[seasonId][4])).days + 1

daysAxis = [namelistOutlooks[seasonId][4] + timedelta(days = float(d)) for d in np.arange(nDays)]


fig, ax = plt.subplots(1, 1, figsize = (6, 4))


seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)

# Load observational data first
# -----------------------------
mainVerif = "NSIDC-0081"
#mainVerif = "OSI-401-b"
alternativeVerif = "OSI-401-b"
#alternativeVerif = "NSIDC-0081"
f = Dataset("../data/" + seasonName + "/netcdf/regrid/" + \
            mainVerif + "_000_concentration_2x2.nc")
sic_obs = f.variables["siconc"][:].data
latitude = f.variables["latitude"][:].data
longitude = f.variables["longitude"][:].data
cellarea  = f.variables["areacello"][:].data
mask_obs  = f.variables["sftof"][:].data
nt = f.dimensions["time"].size
f.close()

if plotOtherVerif:
	# If asked to plot second product, just add it to the namelist
	namelistContributions.append([alternativeVerif, [1, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0], "v"])

siProbGroupList = list() # for group forecast: list of arrays with probability of sea ice presence
# Run through all forecasts

for j, n in enumerate(namelistContributions):
	print(n)
	thisName = n[0]
	thisNbForecasts	= n[seasonId + 1][diagId]
	thisType = n[-1] # Statistical or dynamical

	if thisNbForecasts > 0:

		if thisName == "AWI-SDAP":
			# Special case for AWI-SDAP, that provides directly probability of sea ice presence basedon 15% threshold for SIC
			filein = "../data/" + seasonName + "/netcdf/regrid/" + thisName + "_" \
                                  + "001_probability_2x2.nc"

			f = Dataset(filein, mode = "r")
			siProb = f.variables["siprob"][:] / 100.0 # since provided as %
			f.close()


		else:
			siPresContribList = list()

			for jFor in np.arange(thisNbForecasts):

				if thisName == alternativeVerif:
					thisMember = "000"
				else:
					thisMember = str(jFor + 1).zfill(3)
				filein = "../data/" + seasonName + "/netcdf/regrid/" + thisName + "_" \
					  + thisMember + "_concentration_2x2.nc"
	
				f = Dataset(filein, mode = "r")
				sic = f.variables["siconc"][:]
				f.close()

				# Binarizing SIC based on 15% threshold to define sea ice presence
				siPres = 1.0 * (sic > thresholdSeaIcePresence)
	
				siPresContribList.append(siPres)
	
				del siPres
	
			# Now we are out of the loop of members and we have a list of sea ice presence (binary fields)
			# Let's average through members to obtain a probability of sea ice presence.
			siProb = np.mean(np.array(siPresContribList), axis = 0)
			
			del siPresContribList

		# Now we compute the IIEE for that contribution, based on the 0.5 threshold (Helge Goessling)

		thisIIEE = iiee(siProb, sic_obs, cellarea, mask = (latitude < 0), thresholdProbability = 0.5, thresholdSeaIcePresence = thresholdSeaIcePresence)

		# We also store the probability to compute later a group forecast
		if thisName != "climatology":
			siProbGroupList.append(siProb)
		
		# Plot series
		if thisType == "s":
			thisColor = colorDict["statistical"]
		elif thisType == "d":
			thisColor = colorDict["dynamical"]
		elif thisType == "b": # benchmark
			thisColor = colorDict["climatology"]
		elif thisType == "v": # other verification
			thisColor = colorDict["alternative verification"]

		ax.plot(daysAxis, thisIIEE, color = thisColor, lw = 2)

		#ax.text(daysAxis[-1] + timedelta(days  = 1), thisIIEE[-1], thisName, ha = "left", fontsize = 4, color = [0.5, 0.5, 0.5], va = "center")

		del siProb

# Compute group forecast
siProbGroup = np.mean(np.array(siProbGroupList), axis = 0)

iieeGroup = iiee(siProbGroup, sic_obs, cellarea, mask = (latitude < 0), thresholdProbability = 0.5, thresholdSeaIcePresence = thresholdSeaIcePresence)

ax.plot(daysAxis, iieeGroup, color = colorDict["group forecast"], lw = 3, linestyle = "-")


linesLegend = [Line2D([0], [0], color = colorDict[c], linewidth=3) for c in colorDict]
linesLabels = [c for c in colorDict]
ax.legend(linesLegend, linesLabels)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
ax.set_ylabel("Million km$^2$")
ax.set_title("Integrated Ice Edge Error (" + seasonName + ")")
fig.savefig("../figs/fig6_paper.png", dpi = 300)
