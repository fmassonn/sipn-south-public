# Author: F. Massonnet
# Date  : 10 Jan 2022


# Purpose: Figure with daly sea ice area for one particular season 
#          highlighting the bias at initial day and the excessive melt rates

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

matplotlib.rcParams['font.family'] = "Arial Narrow"

plt.close("all")


# Open namelist information file
# ------------------------------
exec(open("./namelist.py").read())

# Script parameters
# -----------------

seasonId = 5 # Which season to look at: 0 = 2017-2018, 1 = 2018- 2019, 2 = 2019-2020, 3 = 2020-2021, 4 = 2021-2022, 5 = 2022-2023
diagId   = 0 # Which diag (currently only works with total SIA so 0)

# ----
nDaysWeek = 7 # how many days in a week

nDays    = ((namelistOutlooks[seasonId][5]) - (namelistOutlooks[seasonId][4])).days + 1

daysAxis = [namelistOutlooks[seasonId][4] + timedelta(days = float(d)) for d in np.arange(nDays)]

fig, ax = plt.subplots(2, 2, figsize = (8, 6))

# Run through all forecasts
seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)

firstDyn = False
firstSta = False

# All data to make medians lter
groupDyn = list()
groupSta = list()

for j, n in enumerate(namelistContributions):
	print(n)
	thisName = n[0]
	thisNbForecasts	= n[seasonId + 1][diagId]
	thisType = n[-1] # Statistical or dynamical

	if thisNbForecasts > 0:
		thisList = list()
		for jFor in np.arange(thisNbForecasts):
			# Locate the file
			fileIn = "../data/" + seasonName + "/txt/" + thisName + "_" + str(jFor + 1).zfill(3) + "_total-area.txt"

			# Read file content
			csv = pd.read_csv(fileIn, header = None)
			seriesTmp = csv.iloc[0][:nDays]
			thisList.append(seriesTmp)

			#print(j / len(namelistContributions) * 255)
			
			#ax.plot(daysAxis, seriesTmp, color = plt.cm.gnuplot( int(j / len(namelistContributions) * 255)))
		# Compute ensemble mean
		thisMean = np.nanmean(np.array(thisList), axis = 0)
		thisMin = np.nanmax(np.array(thisList), axis = 0)
		thisMax = np.nanmin(np.array(thisList), axis = 0)
		if thisType == "s":
			if not firstSta:
				firstSta = True
				label = "statistical contributions"
			else:
				label = None
			thisColor = "#008579"
			ax[0, 0].plot(daysAxis, thisMean, color = thisColor, label = label, lw = 0.5)
			meltRate = thisMean[nDaysWeek:] - thisMean[:-nDaysWeek]
			ax[0, 1].plot(daysAxis[nDaysWeek:], meltRate, color = thisColor, label = label, lw = 0.5)

			groupSta.append(thisMean)

		elif thisType == "d":
			if not firstDyn:
				firstDyn = True
				label = "dynamical contributions"
			else:
				label = None

			thisColor = "#FF7200"
			ax[1, 0].plot(daysAxis, thisMean, color = thisColor, label = label, lw = 0.5)
			meltRate = thisMean[nDaysWeek:] - thisMean[:-nDaysWeek]
			ax[1, 1].plot(daysAxis[nDaysWeek:], meltRate, color = thisColor, label = label, lw = 0.5)
	
			groupDyn.append(thisMean)
		elif thisType == "b":
			thisColor = [0.2, 0.2, 0.2]
			[ax[jj, 0].plot(daysAxis, thisMean, color = thisColor, label = thisName, linestyle = ":", lw = 2, zorder = 1000) for jj in [0, 1]]
			# Plot melt rate
			meltRate = thisMean[nDaysWeek:] - thisMean[:-nDaysWeek]
			[ax[jj, 1].plot(daysAxis[nDaysWeek:], meltRate, color = thisColor, label = thisName, lw = 2, linestyle = ":", zorder = 1000) for jj in [0, 1]]
		else:
			stop("Type not known")

# Plot ensemble medians
medianSta = np.median(groupSta, axis = 0)
medianDyn = np.median(groupDyn, axis = 0)
ax[0, 0].plot(daysAxis, medianSta, color = "#008579", label = "Group mean", lw = 2)
ax[1, 0].plot(daysAxis, medianDyn, color = "#FF7200", label = "Group mean", lw = 2)

# Plot melt rates of ensemble medians (linear diagnostic so can be done on group means)

meltRateSta = medianSta[nDaysWeek:] - medianSta[:-nDaysWeek]
meltRateDyn = medianDyn[nDaysWeek:] - medianDyn[:-nDaysWeek]
ax[0, 1].plot(daysAxis[nDaysWeek:], meltRateSta, color = "#008579", label = "Group mean", lw = 2)
ax[1, 1].plot(daysAxis[nDaysWeek:], meltRateDyn, color = "#FF7200", label = "Group mean", lw = 2)


# Plot observational references
obsVerif = ["NSIDC-0081", "OSI-401-b"]
obsLine  = ["-", "-"]
firstObs = False

for j, obsname in enumerate(obsVerif):
	seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)

	fileIn = "../data/" + seasonName + "/txt/" + obsname + \
                         "_000_total-area.txt"

	csv = pd.read_csv(fileIn, header = None)
	seriesTmp = csv.iloc[0][:nDays]

	if not firstObs:
		firstObs = True
		label = "observational references"
	else:
		label = None

	ax[0, 0].plot(daysAxis, seriesTmp, label = label, linestyle = obsLine[j], color = "k")
	ax[1, 0].plot(daysAxis, seriesTmp, label = label, linestyle = obsLine[j], color = "k")

	meltRate = np.array(seriesTmp[nDaysWeek:]) - np.array(seriesTmp[:-nDaysWeek])


	ax[0, 1].plot(daysAxis[nDaysWeek:], meltRate, linestyle = obsLine[j], color = "k", lw = 2)
	ax[1, 1].plot(daysAxis[nDaysWeek:], meltRate, linestyle = obsLine[j], color = "k", lw = 2)

[a.legend(fontsize = 9) for a in ax.flatten()]

for a in ax.flatten():
	a.set_xlim(namelistOutlooks[seasonId][4], namelistOutlooks[seasonId][5])
	a.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
for j in [0, 1]:
  ax[j, 0].set_ylim(0.0, 15)
  ax[j, 1].set_ylim(-2.5, 1.0)
  ax[j, 1].plot(daysAxis, np.zeros(len(daysAxis)), color = [0.2, 0.2, 0.2], lw = 1)
  ax[j, 0].set_title("Daily mean Southern Ocean sea ice area, " + seasonName)
  ax[j, 0].set_ylabel("Million km$^2$")
  ax[j, 1].set_title("Weekly running melt rates, " + seasonName)
  ax[j, 1].set_ylabel("Million km$^2$/week")


fig.tight_layout()

plt.savefig("../figs/fig4_paper.png", dpi = 300)

