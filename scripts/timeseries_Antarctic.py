#!/usr/bin/python
#
# Time series of sea ice in February


from   netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import wget
import os
import csv
import pandas as pd
import calendar

# Target month (normal convention: 1 = January)

month = 2
hemi = "south"
diag = "area"

# -----
hemi_region = {"south": "Antarctic",
               "north": "Arctic"   ,
              }
rootdir = "ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/" + hemi +"/monthly/data/"
file    = hemi[0].upper() + "_" + str(month).zfill(2) + "_extent_v3.0.csv"

if os.path.exists(file):
    print("File already exists")
else:
    wget.download(rootdir + file)


# Read the CSV file
df = pd.read_csv(file, sep = ",", skipinitialspace = True)




time = np.array(df["year"])
sie  = np.array(df[diag])

# Trend fit
a = np.polyfit(time, sie, 1)
fit = a[0] * time + a[1]
res = sie - fit
std = np.std(res)



fig, ax = plt.subplots(1, 1, figsize = (6, 3), dpi = 300)

ax.plot(time, sie, color = [0.2, 0.6, 1.0], lw = 3, ls = "-", zorder = -1)
ax.scatter(time, sie, 5, color = [0.1, 0.3, 0.5], marker = "s")
ax.plot(time, fit, color = [0.1, 0.3, 0.5], lw = 2, linestyle = "-")
ax.scatter(time[-1], sie[-1], 200, marker = '*', color = [0.1, 0.3, 0.5],
           zorder = 1000)
ax.fill_between(time,  fit - 2 * std, fit + 2 * std, linewidth = 0, 
                edgecolor = None, facecolor = [0.8, 0.9, 1.0], zorder = -3)
ax.fill_between(time,  fit - 1 * std, fit + 1 * std, linewidth = 0, 
                edgecolor = None, facecolor = [0.6, 0.8, 1.0], zorder = -2)

ax.set_ylim(0.0, 1.2 * np.max(sie))
ax.set_xlim(time[0] - 1, time[-1] + 1)
ax.set_axisbelow(True)
ax.set_ylabel("10$^6$ km$^2$")
ax.set_title(calendar.month_name[month] + " " + hemi_region[hemi] + " sea ice " + diag)
ax.grid()

# Sorted data
sorted_data = [x for x in sorted(zip(sie,time))]
print(sorted_data)

plt.tight_layout()
plt.savefig("../figs/figTimeSeries.png", dpi = 400)

