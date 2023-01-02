# Author: F. Massonnet
# Date  : 04 Oct 2022


# Purpose: Figure with PDFs and CDFs and CRPSs of February mean sea ice area

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

nDays    = ((namelistOutlooks[seasonId][5]) - (namelistOutlooks[seasonId][4])).days + 1
seasonName = str(namelistOutlooks[seasonId][0].year) + "-" + str(namelistOutlooks[seasonId][1].year)



fig, ax = plt.subplots(1, 2, figsize = (6, 5))

# Get value from observations
# Plot observational references
obsVerif = ["NSIDC-0081", "OSI-401-b"]
obsLine  = ["--", ":"]
listVerif = list()
for j, obsname in enumerate(obsVerif):
	fileIn = "../data/" + seasonName + "/txt/" + obsname + \
                         "_000_total-area.txt"
	csv = pd.read_csv(fileIn, header = None)
	seriesTmp = csv.iloc[0][:nDays]
	t1 = (namelistOutlooks[seasonId][2] -  namelistOutlooks[seasonId][4]).days
	t2 = (namelistOutlooks[seasonId][3] -  namelistOutlooks[seasonId][4]).days + 1
	
	thisValue = np.mean(seriesTmp[t1: t2])
	listVerif.append(thisValue)

	ax[0].plot((thisValue, thisValue), (-1e9, 1e9), label = obsVerif[j], linestyle = obsLine[j], color = "k", lw = 1, zorder = -100)

verifValue = listVerif[0]
# Run through all forecasts

# Prepare lists for plotting
listInfo = list()
for j, n in enumerate(namelistContributions):
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

			t1 = (namelistOutlooks[seasonId][2] -  namelistOutlooks[seasonId][4]).days
			t2 = (namelistOutlooks[seasonId][3] -  namelistOutlooks[seasonId][4]).days + 1

			thisValue = np.mean(seriesTmp[t1: t2])
			thisList.append(thisValue)
		if thisType == "s":
			thisColor = plt.cm.YlOrRd_r( int(j / len(namelistContributions) * 128))
		elif thisType == "d":
			thisColor = plt.cm.PuBuGn_r( int(j / len(namelistContributions) * 128))
		else:
			thisColor = "black"
		
		# Compute CRPS
		step = 0.01
		x = np.arange(0, 100.0, step) # x-axis
		cdfVerif = (x > verifValue) * 1
		cdfForecast = np.sum([(x > f) * 1 / thisNbForecasts for f in thisList], axis = 0)


		thisCRPS = np.sum((cdfForecast - cdfVerif) ** 2 * step)

		listInfo.append([thisName, thisCRPS, thisList, thisType, thisColor])

		# Plot PDF to check
		figTmp, axTmp = plt.subplots(2, 1, figsize = (4, 8))
		axTmp[0].plot(x, cdfForecast, label = thisName)
		axTmp[0].plot(x, cdfVerif, color = "k", label = "OBS")
		axTmp[1].set_title("cumulative distribution functions")
		axTmp[0].legend()
		axTmp[0].set_xlim(0.0, 4.0)
		axTmp[0].set_ylim(0.0, 1.1)

		axTmp[1].plot(x, (cdfForecast - cdfVerif) ** 2)
		axTmp[1].set_title("(forecast CDF - obs CDF)$^2$")
		axTmp[1].set_xlim(0.0, 4.0)
		axTmp[1].set_ylim(0.0, 1.1)
		figTmp.savefig(thisName + ".png")
		plt.close(figTmp)


# Print result

print(seasonName)
for l in listInfo:
	print("{:<16} {:} {:<16}".format(l[0], l[3], str(np.round(l[1], 2))))

# Sort by CRPS
sortedList = sorted(listInfo, key = lambda x: x[1])


# Plots
for j, s in enumerate(sortedList):
	thisName   = s[0]
	thisNbForecasts = len(s[2])
	thisCRPS   = s[1]
	theseAreas = s[2]
	thisType   = s[3] # Statistical or dynamical
	thisColor =  s[4]

	#Plot dots
	ax[0].scatter(theseAreas, j * np.ones(thisNbForecasts), 3, marker = "x", color = thisColor, lw = 0.5)
	ax[0].text(-0.5, j, thisName + " (" + thisType + ")", ha = "right", va = "center")


	ax[0].set_yticklabels("")
	ax[0].set_title("February mean sea ice area\nforecast distributions (" + seasonName + ")")
	ax[0].set_xlabel("Million km$^2$")
	# Plot PDF
	if thisNbForecasts > 1:
		xpdf = np.linspace(0, 10, 1000)
		kernel = stats.gaussian_kde(theseAreas)
		pdf = kernel(xpdf).T
		# Normalize PDF max for visual
		pdf /= np.max(pdf)
		ax[0].fill_between(xpdf, j, j + 0.5 * pdf, color = thisColor, alpha = 0.2, lw = 0)

	# Plot CRPS
	ax[1].fill_between((0, thisCRPS), (j - 0.45, j - 0.45), (j + 0.45, j + 0.45), color = thisColor, alpha = 0.5, lw = 0)
	


ax[1].set_yticklabels("")
ax[1].set_title("Continuous rank probability score")
ax[1].set_xlabel("(Million km$^2$)$^2$")
for a in ax:
	a.set_ylim(-0.5, len(sortedList))
ax[0].set_xlabel("Million km$^2$")
ax[0].set_xlim(0, 4)
ax[0].legend(fontsize = 8)
fig.tight_layout()
figName = "../figs/fig5_paper.png"
fig.savefig(figName, dpi = 300)
print("Figure " + figName + " printed")
