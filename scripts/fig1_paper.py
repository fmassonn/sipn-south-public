# Author: F. Massonnet
# Date  : 10 Jan 2022


# Purpose: Overview figure with meta-data

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import pandas            as pd
from   matplotlib import font_manager
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

plt.close("all")


# Change font globally
# --------------------
#font_dirs  = ["/System/Library/Fonts/", ]
#font_files = font_manager.findSystemFonts(fontpaths = font_dirs)
#font_list  = font_manager.createFontList(font_files)
#font_manager.fontManager.ttflist.extend(font_list)
matplotlib.rcParams['font.family'] = "Arial Narrow"


# Open Namelist file
exec(open("./namelist_contributors.py").read())


nYears = yearEnd - yearStart + 1

years = np.arange(yearStart, yearEnd + 1)

# How many contributors (at least one file submitted) over the years?
# 
nCount = np.full(nYears, 0)
nFiles = np.full(nYears, 0)

for j, y in enumerate(years):
    nCount[j] = len([i for i in info if np.max(i[j + 1]) > 0])
    
    nFiles[j] = np.sum([i[j + 1] for i in info])


  
fig, ax1 = plt.subplots(1, 1, figsize = (5, 3), dpi = 300)

# Instantiate new axis system based on ax2
color1 = "#4189DD"
ax1.bar(years, nCount, color = color1, label = "Nb. contributing groups", alpha = 0.9)
ax1.tick_params(axis = "y", color = color1, labelcolor = color1)
ax1.set_ylim(bottom = 0, top = 20)

ax2 = ax1.twinx()

color2 = "#381D59"
ax2.plot(years, nFiles, "-bs", color = color2, label = "Nb. files contributed")
ax2.set_ylim(bottom = 0, top = 800)
ax1.set_yticks(np.arange(0, 25, 5))
ax2.tick_params(axis = "y", color = color2, labelcolor = color2)

myLabels = [str(y) + "-" + str(y +1) for y in years]
ax1.set_xticks(years)
ax1.set_xticklabels(myLabels, rotation = 20)


ax2.set_xlim(2016, 2022)

ax1.set_title("Evolution of input statistics to SIPN South")
ax1.grid(axis = "y")
ax1.set_axisbelow(True)
fig.legend(ncol = 1, bbox_to_anchor=(0.04, 0.78, 0.5, .102))
# Legend

fig.savefig("../figs/fig1_paper.png")