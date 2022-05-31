# Author: F. Massonnet
# Date  : 10 Jan 2022


# Purpose: Figure with Continuous Rank Probability Scores for each region
# and per submission

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import csv
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


# Record all ensemble forecast values

nMemb = 50
t1 = 62
t2 = 90

# Load forecast

siaForecasts = list()

for jMemb in np.arange(1, nMemb + 1):
  fileIn = "/Users/massonnetf/git/sipn-south-public/data/2021-2022/txt/cmcc_" \
      + str(jMemb).zfill(3) + "_total-area.txt"
  
  csvFile = open(fileIn,)
  csvReader = csv.reader(csvFile, delimiter = ",")
    
  for row in csvReader:
    dataTmp = [float(d) for d in row]
    
    # Make mean over period
    timeMeanTmp = np.mean(dataTmp[t1:t2])
    siaForecasts.append(timeMeanTmp)
    
  csvFile.close()
  
# Sort sia forecasts (important for CRPS)
siaForecasts.sort()
  
# Load obs
fileIn = "/Users/massonnetf/git/sipn-south-public/data/2021-2022/txt/NSIDC-0081_000_total-area.txt"
csvFile = open(fileIn)
csvReader = csv.reader(csvFile, delimiter = ",")
    
for row in csvReader:
  dataTmp = [float(d) for d in row]
    
# Make mean over period
obsRef = np.mean(dataTmp[t1:t2])
csvFile.close()
print(obsRef)



# Compute CRPS
fig, ax = plt.subplots(1, 1, figsize = (4, 3))
ax.scatter(siaForecasts, np.zeros(len(siaForecasts)), 50, marker = "x", color = "blue")
# Draw CDF
listXpoints = [0] + siaForecasts + [100]
listCDF     = [0] + [i / len(siaForecasts) for i in range(1, len(siaForecasts) + 1 )] + [1]
[ax.plot((listXpoints[i], listXpoints[i + 1]), (listCDF[i], listCDF[i]), 'b-') for i in range(len(siaForecasts) + 1)]
[ax.plot((listXpoints[i + 1], listXpoints[i + 1]), (listCDF[i], listCDF[i + 1]), "b-") for i in range(len(siaForecasts) + 1) ]
ax.set_title("Cumulative forecast distribution for\nucl sea ice area February 2022")


# Plot obs
ax.plot((obsRef, obsRef), (-1000, 1000), "r-")

ax.set_xlim(1.2, 2.7)
ax.set_ylim(0.0, 1.0)
fig.tight_layout()
fig.savefig("./CRPS.png", dpi = 300)



