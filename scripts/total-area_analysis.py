#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : December 2017
# Data: SIPN South contributors

# Imports
import pandas as pd
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import Dataset
from datetime import date

# Load namelist
exec(open("./namelist.py").read())

# Time axis
time = pd.date_range(pd.to_datetime("20180201", format = "%Y%m%d"), periods = 28).tolist()

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# Number of forecasts for each submission
n_for = [info[j_sub][1] for j_sub in range(n_sub)]

# Colors
col = [info[j_sub][2] for j_sub in range(n_sub)]

plt.figure(figsize = (6, 4))

for j_sub in range(n_sub):
  #plt.text(29, 5.0 - j_sub * 0.3, info[j_sub][0] + " " + info[j_sub][3], color = col[j_sub])
  for j_for in range(n_for[j_sub]):
    # Total area
    filein = "../data/txt/" + sub_id[j_sub] + "_" + str(j_for + 1).zfill(3) + "_total-area.txt"
    # Read the CSV file
    csv = pd.read_csv(filein, header = None)
    series = csv.iloc[0][:]
    # Plot series, line thinner for large ensembles. Legend only if first member
     
    if j_for == 0:
      mylab = info[j_sub][0] + " " + info[j_sub][3]
    else:
      mylab = "_nolegend_"

    plt.plot(np.arange(len(time)) + 1, series, color = col[j_sub], lw = (5.0 / n_for[j_sub]) ** 0.2, label = mylab)

# Plot observations
# -----------------
filein = "../data/txt/NSIDC-0081_000_total-area.txt"
csv = pd.read_csv(filein, header = None)
series = csv.iloc[0][:]
plt.plot(np.arange(len(time)) + 1, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = "--", label = "OBS (NSIDC-0081)")

filein = "../data/txt/OSI-401-b_000_total-area.txt"
csv = pd.read_csv(filein, header = None)
series = csv.iloc[0][:]
plt.plot(np.arange(len(time)) + 1, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = ":", label = "OBS (OSI-401-b)")

plt.legend(fontsize = 8)

plt.title("Forecast of February 2018 total\nAntarctic sea ice area")
plt.xticks([1, 5, 10, 15, 20, 25, 28], ["01", "05", "10", "15", "20", "25", "28"])
plt.xlabel("Day in February 2018")
plt.ylim(0.0, 5)
plt.xlim(0.0, 46)
plt.ylabel("10$^6$ km$^2$")
plt.grid()
plt.tight_layout()
plt.savefig("../figs/fig1.png", dpi = 500)
print("Figure ../figs/fig1.png printed")
