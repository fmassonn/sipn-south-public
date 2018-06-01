#!/usr/bin/python
#
# Time series of sea ice

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

f = Dataset("../data/netcdf/NSIDC-G02135_000_extent.nc")
sie = f.variables["siextents"][:]
f.close()

time = np.arange(1979, 2018 + 1)
sie  = sie[2 + 12::12]

# Trend fit
a = np.polyfit(time, sie, 1)
fit = a[0] * time + a[1]
res = sie - fit
std = np.std(res)



plt.figure(figsize = (5, 3))
plt.plot(time, sie, color = [0.2, 0.6, 1.0], lw = 3)
#plt.plot(time, fit, color = [0.2, 0.6, 1.0], lw = 1, linestyle = "--")
plt.scatter(time[-1], sie[-1], 100, marker = '*', color = [0.2, 0.35, 1.0], zorder = 1000)
plt.fill_between(time,  fit - 2 * std, fit + 2 * std, linewidth = 0, edgecolor = None, facecolor = [0.8, 0.9, 1.0], zorder = 0)
plt.fill_between(time,  fit - 1 * std, fit + 1 * std, linewidth = 0, edgecolor = None, facecolor = [0.6, 0.8, 1.0], zorder = 0)

plt.ylim(0, 6.0)
plt.xlim(1978, 2020)
plt.ylabel("10$^6$ km$^2$")
plt.title("February Antarctic sea ice extent")
plt.tight_layout()
plt.grid()
plt.savefig("../figs/sie.png", dpi = 400)

