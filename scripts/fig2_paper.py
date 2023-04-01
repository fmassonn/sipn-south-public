# Author: F. Massonnet
# Date  : 10 Jan 2022


# Purpose: Overview figure with meta-data about submissions

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

#from seaice_commondiags import *

matplotlib.rcParams['font.family'] = "Arial Narrow"

plt.close("all")



# Open Namelist file
# ------------------
exec(open("./namelist.py").read())


# Fetch starting dates and their number
startDates  = [n[0] for n in namelistOutlooks]
endDates    = [n[1] for n in namelistOutlooks]

nStartDates = len(startDates)

# The list of years defining the start dates
startYears   = [s.year for s in startDates]
endYears     = [e.year for e in endDates]

# The name of seasons (ex: 2017-2018)
nameSeasons = [str(s) + "-" + str(e) for s, e in zip(startYears, endYears)]

# Number of identified contributors to SIPN South
nContributors = len(namelistContributions)

# How many contributors (at least one file submitted) as a function of
# start date?

nCount = list()
nFiles = list()

for s in range(nStartDates):
    # We are in a start date. Let's create two counters, one that 
    # accumulates how many of the identified contributors have indeed
    # submitted at least one file, and one that accumulates the total
    # number of files
    nCountTmp = 0
    nFilesTmp = 0
    
    for c in range(nContributors):
        if max(namelistContributions[c][s + 1]) > 0:
            nCountTmp += 1
        
        nFilesTmp += sum(namelistContributions[c][s + 1])

    nCount.append(nCountTmp)
    nFiles.append(nFilesTmp)
    

  
fig, ax1 = plt.subplots(1, 1, figsize = (4, 3), dpi = 300)

# Instantiate new axis system based on ax2
color1 = "#4189DD"
ax1.bar(startYears, nCount, color = color1, label = "Nb. contributing groups", alpha = 0.9)
ax1.tick_params(axis = "y", color = color1, labelcolor = color1)
ax1.set_ylim(bottom = 0, top = 25)

ax2 = ax1.twinx()

color2 = "#381D59"
ax2.plot(startYears, nFiles, "-s", color = color2, label = "Nb. files contributed")
ax2.set_ylim(bottom = 0, top = 1000)
ax1.set_yticks(np.arange(0, 25, 5))
ax2.tick_params(axis = "y", color = color2, labelcolor = color2)

ax1.set_xticks(startYears)
ax1.set_xticklabels(nameSeasons, rotation = 20)


ax2.set_xlim(2016, 2023)

ax1.set_title("Evolution of input statistics to SIPN South")
ax1.grid(axis = "y")
ax1.set_axisbelow(True)
fig.legend(ncol = 1, bbox_to_anchor=(0.04, 0.78, 0.5, .102))

# Legend
fig.tight_layout()
fig.savefig("../figs/fig2_paper.png", dpi = 300)
