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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import os

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
figName = "../figs/fig1a_paper.png"
fig.savefig(figName)
print("Figure " + figName + " printed")

# Print map Antarctica

# Load the 2023 data SH
file = "/Users/massonnetf/CLIMDATA/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/siconc_SIday_OSI-401-b_r1i1p1_20221201-20230228_sh.nc"
f = Dataset(file, mode = "r")
sic = f.variables["siconc"][:]

sicFeb2023 = np.mean(sic[-28:, :, :], axis = 0)
    # The time indices in the daily file
sicFeb2023[sicFeb2023 < 15] = np.nan
lon1= f.variables["longitude"][:]

lat1 = f.variables["latitude"][:]

f.close()

# Mean
file = "/Users/massonnetf/CLIMDATA/obs/ice/siconc//OSI-SAF/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_197901-201512_sh.nc"

f = Dataset(file, mode = "r")
sic = f.variables["siconc"][:]
sicFebMean = np.mean(sic[(1979-1979)*12+1:(2015-1979)*12+1+1:12],axis = 0)
#sicFebMean[sicFebMean < 15] = np.nan
lat2= f.variables["lat"][:]
lon2 = f.variables["lon"][:]
lon2[lon2<0]+=360 # Tweak to avoid an horrible hexagon appear


f.close()



# Compute a circle in axes coordinates,
# which we can use as a boundary for the map.
# https://scitools.org.uk/cartopy/docs/latest/gallery/lines_and_polygons/

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


fig, ax = plt.subplots(figsize=(3, 3))
plt.axis('off')
myProj = ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0)
ax = plt.axes(projection= \
            myProj)
xlims = [-180,180]
ylims = [-90,-50]
ax.set_extent(xlims+ylims, crs=ccrs.PlateCarree())
ax.set_boundary(circle, transform=ax.transAxes)
ax.stock_img()
ax.coastlines("110m", linewidth=0.5, color="black")
cs = ax.pcolormesh(lon1, lat1, sicFeb2023, \
                   transform=ccrs.PlateCarree(), cmap = plt.cm.Greys_r
                   )
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.contour(lon2, lat2, sicFebMean, [15.1, 15.2], \
           transform=ccrs.PlateCarree(), colors = "#A2DBEB", linewidths = 3)

ax.text(-52, -75, "2023", transform=ccrs.PlateCarree(), fontsize = 10, \
        ha = "center" )
ax.text(-165, -65, "1979-2015", transform=ccrs.PlateCarree(), fontsize = 10, \
        ha = "center", color =  "#A2DBEB", fontweight = "bold" )


# Draw sectors
exec(open("./sectorsDefinitions.py").read())
for j, s in enumerate(sectors):
	if s[0] != "Southern Ocean":
		print(s[0])
		ax.plot((s[1] + 360, s[1] + 360), (-90, 0), transform=ccrs.PlateCarree(), color = [0.3, 0.3, 0.3])
		if s[2] < s[1]:
			meanLon = (s[1] + s[2] + 360) / 2
		else:
			meanLon = (s[1] + s[2]) / 2

		print(meanLon)
		if meanLon < -90 or meanLon > 90:
			angleText = 180 - meanLon
		else:
			angleText = -meanLon
		ax.text(meanLon, -48, s[0], transform = ccrs.PlateCarree(), rotation = angleText, ha = "center", va = "center", fontsize = 7)

figName = "../figs/fig1b_paper.png"
fig.savefig(figName, dpi = 300)
print("Figure " + figName + " printed")


# Merge figures
os.system("convert +append ../figs/fig1a_paper.png ../figs/fig1b_paper.png ../figs/fig1.png")
os.system("convert ../figs/fig1.png ../figs/fig1.JPEG")
figName = "../figs/fig1_paper.png"
print("Figure " + figName + " printed")
