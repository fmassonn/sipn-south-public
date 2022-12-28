#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 10:57:54 2021

@author: massonnetf
"""

# Fig. 1 of the SIPN paper with anomalies of sea ice extent from the OSISAF SH sea ice index

from netCDF4 import Dataset
from datetime import date, datetime, timedelta
import matplotlib
matplotlib.rcParams['font.family'] = "Arial Narrow"

import matplotlib.pyplot as plt
import numpy as np

# To update:
# wget https://thredds.met.no/thredds/fileServer/osisaf/met.no/ice/index/v2p1/sh/osisaf_sh_sie_monthly.nc

filein = "./osisaf_sh_sie_monthly.nc"

f = Dataset(filein, mode = "r")
dateRefString = f.variables["time"].units.split(" ")[2]
time = f.variables["time"][:]
sie  = f.variables["sie"][:]
fillV = f.variables["sie"]._FillValue

f.close()

dateRef = datetime.strptime(dateRefString, "%Y-%m-%d")
dates = [dateRef + timedelta(days = float(d)) for d in time]


# Screen out masked data
sie = sie.data
sie[sie == fillV] = np.nan

dates = [d for j, d in enumerate(dates) if ~np.isnan(sie[j])]
sie   = [s for j, s in enumerate(sie)   if ~np.isnan(sie[j])]

# Compute seasonal cycle
yearb, yeare = 1981, 2010

cycle = np.full(12, np.nan)
for m in np.arange(1, 12 + 1):
    cycle[m - 1] = np.nanmean(np.array([e for e, d in zip(sie, dates) if \
                         d.month == m \
                             and d.year >= yearb \
                                 and d.year <= yeare]))

# Compute anomalies
# First "tile" the cycle
cycleTile = np.array([cycle[d.month - 1] for d in dates])
sie_ano = sie - cycleTile


# Trend
p = np.polyfit(np.arange(len(dates)), sie_ano, 1)
fit = np.arange(len(dates)) * p[0] + p[1]

fig, ax = plt.subplots(1, 1, figsize = (6, 2.5), dpi = 300)

ymin, ymax = -1.5, 1.5
for j, d in enumerate(dates):
    if not np.ma.is_masked(sie_ano[j]):
      color = plt.cm.RdBu(int((sie_ano[j]- ymin) * 255 / (ymax - ymin)))[:3]
      ax.bar(d, sie_ano[j], color = color, width = 50)

#ax.plot(dates, fit, "b--")
ax.plot((dates[0], dates[-1]), (0, 0), "k", lw = 1)

ax.set_title("Antarctic monthly sea ice extent\nanomalies relative to " + \
             str(yearb) + "-" + str(yeare) + " average")
ax.set_ylabel("10$^6$ km$^2$")
#ax.text(dates[-1] + timedelta(days = 930), -2.3, "Data: OSISAF sea ice index ; Figure: @FMassonnet", va = 'bottom', rotation = 90, fontsize = 6)
fig.tight_layout()
figName = "../figs/fig1_paper.png"
fig.savefig(figName)
print("Figure " + figName + " printed")
