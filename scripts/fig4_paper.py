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

seasonId = 4 # Which season to look at
diagId   = 0 # Which diag (currently only works with total SIA so 0)

# ----
nDaysWeek = 7 # how many days in a week

nDays    = ((namelistOutlooks[seasonId][5]) - (namelistOutlooks[seasonId][4])).days + 1

daysAxis = [namelistOutlooks[seasonId][4] + timedelta(days = float(d)) for d in np.arange(nDays)]

fig, ax = plt.subplots(1, 2, figsize = (10, 4))

# Run through all forecasts
seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)

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
		# Compute ensemble median
		thisMedian = np.median(np.array(thisList), axis = 0)
		thisMin = np.max(np.array(thisList), axis = 0)
		thisMax = np.min(np.array(thisList), axis = 0)
		if thisType == "s":
			thisColor = plt.cm.YlOrRd( int(128 + j / len(namelistContributions) * 128))
		elif thisType == "d":
			thisColor = plt.cm.PuBuGn( int(128 + j / len(namelistContributions) * 128))
		else:
			stop("Type not known")
		ax[0].plot(daysAxis, thisMedian, color = thisColor, label = thisName + " (" +  thisType + ")")

		# Display range
		#ax.fill_between(daysAxis, thisMin, thisMax, color = thisColor, alpha = 0.5, lw = 0)

		# Plot melt rate
		meltRate = thisMedian[nDaysWeek:] - thisMedian[:-nDaysWeek]
		ax[1].plot(daysAxis[nDaysWeek:], meltRate, color = thisColor, label = thisName)




# Plot observational references
obsVerif = ["NSIDC-0081", "OSI-401-b"]
obsLine  = ["--", ":"]
for j, obsname in enumerate(obsVerif):
	seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)

	fileIn = "../data/" + seasonName + "/txt/" + obsname + \
                         "_000_total-area.txt"

	csv = pd.read_csv(fileIn, header = None)
	seriesTmp = csv.iloc[0][:nDays]

	ax[0].plot(daysAxis, seriesTmp, label = obsVerif[j], linestyle = obsLine[j], color = "k")

	meltRate = np.array(seriesTmp[nDaysWeek:]) - np.array(seriesTmp[:-nDaysWeek])


	ax[1].plot(daysAxis[nDaysWeek:], meltRate, linestyle = obsLine[j], color = "k", lw = 2)

ax[1].plot((namelistOutlooks[seasonId][4], namelistOutlooks[seasonId][5]), (0, 0), color = "grey")
ax[0].legend(ncol = 2, fontsize = 9)

for a in ax:
	a.set_xlim(namelistOutlooks[seasonId][4], namelistOutlooks[seasonId][5])
	a.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
ax[0].set_title("Daily mean Southern Ocean sea ice area, " + seasonName)
ax[0].set_ylabel("Million km$^2$")
ax[1].set_title("Weekly running melt rates, " + seasonName)
ax[1].set_ylabel("Million km$^2$/week")

fig.tight_layout()

plt.savefig("./fig4_paper.png", dpi = 300)

