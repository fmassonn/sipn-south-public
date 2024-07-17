#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:29:28 2022

@author: massonnetf
"""


# CRPS and rank histograms for SIPN South

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

from colInterpolatOr import *

matplotlib.rcParams['font.family'] = "Arial Narrow"

# Script parameters
# -----------------

# label with the year investigated (2017-2018, 2018-2019, ...)
myyear = "2023-2024"   
stringPeriod=myyear[:4] + "1201-" + myyear[5:9] + "0228"


# Name of observational products. Evaluation will be done with dataset 1   
obs = ["NSIDC-0081", "OSI-401-b"]

# Resolution for saving figures
dpi = 300

# End script parameters
# ---------------------

# Function CRPS
def computeCRPS(forecast, verifValue, step = 0.01, minVal = 0.0, maxVal = 100.0, isRange = False):
	""" 	forecast is a list of ensemble forecasts
		verifValue is the verifying observation
		step is the spacing to estimate the CDFs
		minVal and maxVal are the domain limits of the variables
		isRange is False when the PDF must be reconstructed by assuming that all
		values are realizations from an unknown PDF, while it is True when the PDF
		must be a rectangular shape starting at the lowest member and ending at the highest member
		(like Nico Sun)
	"""
	thisNbForecasts = len(forecast)
	x = np.arange(minVal, maxVal, step) # x-axis
	cdfVerif = (x > verifValue) * 1
	if not isRange:
		cdfForecast = np.sum([(x > f) * 1 / thisNbForecasts for f in forecast], axis = 0)
	else:
		cdfForecast = 1 * (x > np.max(forecast)) + (x - np.min(forecast)) / (np.max(forecast) - np.min(forecast)) * (x <= np.max(forecast)) * (x > np.min(forecast))
	print(cdfForecast)
	CRPS = np.sum((cdfForecast - cdfVerif) ** 2 * step)

	return x, CRPS

# Read or create meta-data
# ------------------------

if myyear == "2017-2018":
  # Initialization date
  inidate = "20180201"
  # Number of days in the forecast period
  ndays   = 28 
  # Label for period that is forecasted
  period_name = "February 2018"
  # Starting and ending time indices (Python conventions)
  t1, t2 = 0, 0 + ndays

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


# Load namelist
exec(open("./namelist_" + myyear + ".py").read())

# Create time axis
time = [datetime.strptime(inidate, "%Y%m%d") + timedelta(days = x) for x in 
                    range(0, ndays)]
nt = len(time)

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# Range of forecasts for each submission
range_for = [info[j_sub][1] for j_sub in range(n_sub)]

# Number of forecasts for each submission
n_for    = [len(l) for l in range_for]

# Colors for each submission
# Create colors
sourceColors = ["#000075", "#645ccc", "#4363d8", "#1898e0", "#00b2ed", "#00bb62", \
               "#8bcd45", "#dbe622", "#f9c410", \
               "#f89e13", "#fb4c27", "#fb4865", \
               "#d24493", "#8f57bf", "#911eb4"]
col = colInterpolatOr(sourceColors, nTarget = n_sub, colorCode = "HEX")

# Change clim's color to black
for j_sub in range(n_sub):
  if sub_id[j_sub] == "climatology":
    col[j_sub] = "black"
    orderSub = 1000


# Load the data in
# The daily total areas will be stored in a list. Each item of the list
# will correspond to one submission and will be a 2-D numpy array
# of dimensions {time, number of ensemble forecasts}
data = list()


for j_sub in range(n_sub):
    print("Reading " + sub_id[j_sub])
    # Create empty numpy array
    sub_data = list()
    # Read in data for that submission
    # Note that j_for refers to the forecast index in non-Pythonic
    # conventions
    for j_for in range_for[j_sub]:
      filein = "../data/" + myyear + "/txt/" + sub_id[j_sub] + "_" + \
               str(j_for).zfill(3) + "_" + stringPeriod + "_total-area.txt"
      # Read the CSV file
      csv = pd.read_csv(filein, header = None)
      series = csv.iloc[0][:]
      # Append that series to the contribution data
      if len(series) != nt:
        print("WARNING: INPUT SERIES TOO LONG, CROPPING")
        print(sub_id[j_sub])
        stop()
      
      # Compute time-mean  
        
      siaMean = np.mean(series[t1:t2])
      sub_data.append(siaMean)
      del series, csv
    
    # Store numpy array
    data.append(sub_data)
    # Delete temporary data for that submission
    del sub_data

# Repeat with observations

# A list of 1-D numpy arrays, each of dimensions {time}
data_obs= list()

for obsname in obs:
    filein = "../data/" + myyear + "/txt/" + obsname + \
    "_000_total-area.txt"
    # Read the CSV file
    csv = pd.read_csv(filein, header = None)
    series = csv.iloc[0][:]

    siaMeanObs = np.mean(series[t1:t2])
    data_obs.append(siaMeanObs)
    
    del series, csv,

# Now compute CRPS

CRPS  = np.full(n_sub + 1, np.nan) # + 1 to accommodate the multi model mean

for j_sub in range(n_sub):

  thisList = data[j_sub]
  obsRef = data_obs[0]

  if sub_id[j_sub] == "NicoSun":
    isRange = True
  else:
    isRange = False

  x, thisCRPS = computeCRPS(thisList, obsRef, step = 0.01, minVal = 0.0, maxVal = 100.0, isRange = isRange)

  CRPS[j_sub] = thisCRPS


  siaForecasts = list(np.sort(data[j_sub]))
  
  # Plot cdf
  fig, ax = plt.subplots(1, 1, figsize = (4, 5))
  ax.scatter(siaForecasts, np.zeros(len(siaForecasts)), 50, marker = "x", color = "blue")
  # Draw CDF
  listXpoints = [0] + siaForecasts + [100]
  listCDF     = [0] + [i / len(siaForecasts) for i in range(1, len(siaForecasts) + 1 )] + [1]
  [ax.plot((listXpoints[i], listXpoints[i + 1]), (listCDF[i], listCDF[i]), 'b-') for i in range(len(siaForecasts) + 1)]
  [ax.plot((listXpoints[i + 1], listXpoints[i + 1]), (listCDF[i], listCDF[i + 1]), "b-") for i in range(len(siaForecasts) + 1) ]

  # Plot obs
  ax.plot((obsRef, obsRef), (-1000, 1000), "r-")
  ax.set_ylim(-0.05, 1.05)
  ax.set_xlim(-0.05, 2.5)
  
  ax.set_title(sub_id[j_sub])
  
  fig.savefig("../figs/CDF_" + sub_id[j_sub] + ".png", dpi = 300)


# Create the multi-model-mean forecast
mmForecast = [np.mean(d) for j_sub, d in enumerate(data) if sub_id[j_sub] == "climatology"]
# Compute CRPS
_, mmCRPS = computeCRPS(mmForecast, obsRef, step = 0.01, minVal = 0.0, maxVal = 100.0)
print("MM forecast: " + str(mmCRPS)) 
CRPS[n_sub] = mmCRPS
sub_id += ["group forecast"]
col  += ["blue"]
  
zipped = zip(sub_id, CRPS, col, sub_id)

zipped_sorted = sorted(zipped, key = lambda x: x[1])

# Plot result

fig, ax = plt.subplots(1, 1, figsize = (3, 4))
for j, z in enumerate(zipped_sorted):
    
    ax.fill_between((0, z[1]), (n_sub + 1 - j, n_sub + 1 - j), color = z[2], alpha = 1.0)
    ax.text(z[1], n_sub + 1 - j - 0.5, "  " + z[3], color = z[2], ha = 'left', va = "center")
ax.set_xlim(0.0, 2.5)
ax.set_ylim(-0.05, n_sub + 1)
ax.set_title("Continuous rank probability score\nfor " + target_period_name + " " + str(myyear[5:]) + " total sea ice area")
ax.set_xlabel("(Million km$^2$)$^2$")
ax.set_yticklabels("")

fig.tight_layout()
fig.savefig("../figs/CRPS.png", dpi = 300)
