#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : December 2017
#         December 2018: update to account for new data arrival
#                        each year
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

plt.close("all")

myyear = "2018-2019"  # label with the year investigated (2017-2018, 2018-2019, ...)
plotobs = False       # Add obs as reference or not (False if forecast mode)

# Load namelist
exec(open("./namelist_" + myyear + ".py").read())

if myyear == "2017-2018":
  inidate = "20180201"
  ndays   = 28 #number of days of the forecast period
  period_name = "February 2018"

elif myyear == "2018-2019":
  inidate = "20181201"
  ndays   = 90
  period_name = "Dec-Jan-Feb 2018-2019"

# Time axis
time = pd.date_range(pd.to_datetime(inidate, format = "%Y%m%d"), periods = ndays).tolist()

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# Number of forecasts for each submission
n_for = [info[j_sub][1] for j_sub in range(n_sub)]

# Colors
col = [info[j_sub][2] for j_sub in range(n_sub)]

# Time series
plt.figure("fig1", figsize = (6, 4))
# Minimum 
plt.figure("fig2", figsize = (6, 4))
# Monthly means
plt.figure("fig3", figsize = (6, 4))

for j_sub in range(n_sub):
  for j_for in range(n_for[j_sub]):
    # Total area
    filein = "../data/" + myyear + "/txt/" + sub_id[j_sub] + "_" + str(j_for + 1).zfill(3) + "_total-area.txt"
    # Read the CSV file
    csv = pd.read_csv(filein, header = None)
    series = csv.iloc[0][:]
    # Plot series, line thinner for large ensembles. Legend only if first member
     
    if j_for == 0:
      mylab = info[j_sub][0] + " " + info[j_sub][3]
    else:
      mylab = "_nolegend_"

    plt.figure("fig1")
    plt.plot(time, series, color = col[j_sub], lw = (5.0 / n_for[j_sub]) ** 0.2, label = mylab)


    # Record when the minimum of the series occurs
    # --------------------------------------------
    # Fit a quadratic polynomial to filter out weather variability
    tt = np.arange(0.0, len(series), 0.1)
    daymin = 1 + tt[np.argmin(np.polyval(np.polyfit(range(len(series)), series, 2), tt))]

    plt.figure("fig2")
    if info[j_sub][3] != "(monthly)":
      plt.scatter(daymin, n_sub - j_sub, color = col[j_sub], label = mylab)
      #plt.text(46, n_sub - j_sub, info[j_sub][0] + " " + info[j_sub][3], color = col[j_sub], ha = "right")

    plt.figure("fig3")
    plt.scatter(np.mean(series), n_sub - j_sub, color = col[j_sub], label = mylab)

# Plot observations
# -----------------
if plotobs:

  filein = "../data/" + myyear + "/txt/NSIDC-0081_000_total-area.txt"
  csv = pd.read_csv(filein, header = None)
  series = csv.iloc[0][:]
  plt.figure("fig1")
  plt.plot(np.arange(len(time)) + 1, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = "--", label = "OBS (NSIDC-0081)")
  plt.figure("fig2")
  tt = np.arange(0.0, len(series), 0.1)
  daymin = 1 + tt[np.argmin(np.polyval(np.polyfit(range(len(series)), series, 2), tt))]
  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS (NSIDC-0081)")
  plt.figure("fig3")
  plt.plot((np.mean(series), np.mean(series)), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS (NSIDC-0081)")
  
  
  
  filein = "../data/" + myyear + "/txt/OSI-401-b_000_total-area.txt"
  csv = pd.read_csv(filein, header = None)
  series = csv.iloc[0][:]
  plt.figure("fig1")
  plt.plot(np.arange(len(time)) + 1, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = ":", label = "OBS (OSI-401-b)")
  plt.legend(fontsize = 8)
  plt.figure("fig2")
  tt = np.arange(0.0, len(series), 0.1)
  daymin = 1 + tt[np.argmin(np.polyval(np.polyfit(range(len(series)), series, 2), tt))]
  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS (OSI-401-b)")
  plt.legend(fontsize = 8)
  plt.figure("fig3")
  plt.plot((np.mean(series), np.mean(series)), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS (OSI-401-b)")
  plt.legend(fontsize = 8)
  
plt.figure("fig1")
plt.title(period_name + " total Antarctic sea ice area")
#plt.xticks([1, 5, 10, 15, 20, 25, 28], ["01", "05", "10", "15", "20", "25", "28"])
#plt.xlabel("Day in February 2018")
plt.ylim(0.0, 10)

#"plt.xlim(0.0, 46)
import matplotlib.dates as mdates
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.ylabel("10$^6$ km$^2$")
plt.grid()
plt.tight_layout()
plt.savefig("../figs/fig1.png", dpi = 500)
print("Figure ../figs/fig1.png printed")
  
  
  #plt.figure("fig2")
  #plt.title("When did the 2018 minimum of Antarctic sea ice area occur?")
  #plt.xlim(0.0, 46)
  #plt.xticks([1, 5, 10, 15, 20, 25, 28], ["<01", "05", "10", "15", "20", "25", ">28"])
  #plt.xlabel("Day in February 2018")
  #plt.yticks([],[])
  #plt.tight_layout()
  #plt.savefig("../figs/fig2.png", dpi = 500)
  #print("Figure ../figs/fig2.png printed")
  
  
  #plt.figure("fig3")
  #plt.title("Monthly mean Antarctic sea ice area")
  #plt.xlim(0.0, 5.0)
  #plt.xticks([1, 5, 10, 15, 20, 25, 28], ["<01", "05", "10", "15", "20", "25", ">28"])
  #plt.xlabel("10$^6$ km$^2$")
  #plt.yticks([],[])
  #plt.tight_layout()
  #plt.savefig("../figs/fig3.png", dpi = 500)
  #print("Figure ../figs/fig3.png printed")
  
