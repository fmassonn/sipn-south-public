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
plotobs = True        # Add obs as reference or not (False if forecast mode)

# Load namelist
exec(open("./namelist_" + myyear + ".py").read())

if myyear == "2017-2018":
  inidate = "20180201"
  ndays   = 28 #number of days of the forecast period
  period_name = "February 2018"
  t1, t2 = 0, 0 + ndays # time indices defining the period for diagnostics
                        # (Pythonic convention)

elif myyear == "2018-2019":
  inidate = "20181201"
  ndays   = 90
  period_name = "Dec-Jan-Feb 2018-2019"
  t1, t2 = 63 - 1, 63 - 1 + 28    # Period to compute means and identify minimum
                                  # nb of days since first day

# Time axis
time = pd.date_range(pd.to_datetime(inidate, format = "%Y%m%d"), periods = ndays).tolist()

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# List of forecasts for each submission
list_for = [info[j_sub][1] for j_sub in range(n_sub)]
n_for    = [len(l) for l in list_for] # Nb of forecasts
# Colors
col = [info[j_sub][2] for j_sub in range(n_sub)]

# Time series
plt.figure("fig1", figsize = (6, 4))
# Minimum 
plt.figure("fig2", figsize = (6, 4))
# Monthly means
plt.figure("fig3", figsize = (6, 4))

for j_sub in range(n_sub):
  for j_for in list_for[j_sub]:
    # Total area
    filein = "../data/" + myyear + "/txt/" + sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_total-area.txt"
    # Read the CSV file
    csv = pd.read_csv(filein, header = None)
    series = csv.iloc[0][:]
    # Plot series, line thinner for large ensembles. Legend only if first member
     
    if j_for == list_for[j_sub][0]:
      mylab = info[j_sub][0] + " " + info[j_sub][3]
    else:
      mylab = "_nolegend_"

    plt.figure("fig1")
    plt.plot(time, series, color = col[j_sub], lw = (5.0 / n_for[j_sub]) ** 0.2, label = mylab)


    # Record when the minimum of the series occurs
    # --------------------------------------------
    # Fit a quadratic polynomial to filter out weather variability
    tt = np.arange(0.0, len(series[t1:t2]), 1)
    daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]

    plt.figure("fig2")
    if info[j_sub][3] != "(monthly)":
      plt.scatter(daymin, n_sub - j_sub, color = col[j_sub], label = mylab)
      #plt.text(46, n_sub - j_sub, info[j_sub][0] + " " + info[j_sub][3], color = col[j_sub], ha = "right")

    plt.figure("fig3")
    plt.scatter(np.mean(series), n_sub - j_sub, color = col[j_sub], label = mylab)

# Plot observations
# -----------------
if plotobs:

  filein = "../data/2017-2018/txt/NSIDC-0081_000_total-area.txt"
  csv = pd.read_csv(filein, header = None)
  series = csv.iloc[0][:]
  plt.figure("fig1")
  plt.plot(time, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = "--", label = "OBS 2017-2018 (NSIDC-0081)")
  #Also plot (hard-coded for now...) the current state
  mydata = [8.4169,8.2069,7.9816,7.8211,7.6011,7.2880,7.1140,6.9764,6.7289,6.4839,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
  plt.plot(time, mydata, color = [1.0, 0.2, 0.2], lw = 3.0, linestyle = "--", label = "OBS 2018-2019 (NSIDC-0081)")

  plt.figure("fig2")
  tt = np.arange(0.0, len(series[t1:t2]), 1)
  daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS 2017-2018 (NSIDC-0081)")
  plt.figure("fig3")
  plt.plot((np.mean(series), np.mean(series)), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS 2017-2018 (NSIDC-0081)")
  
  
  
  filein = "../data/2017-2018/txt/OSI-401-b_000_total-area.txt"
  csv = pd.read_csv(filein, header = None)
  series = csv.iloc[0][:]
  plt.figure("fig1")
  plt.plot(time, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = ":", label = "OBS 2017-2018 (OSI-401-b)")
  #Also plot (hard-coded for now...) the current state
  mydata = [8.9798,8.7607,8.4994,8.2861,8.0624,7.7428,7.5236,7.2674,7.0243,6.7677,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
  plt.plot(time, mydata, color = [1.0, 0.2, 0.2], lw = 3.0, linestyle = ":", label = "OBS 2018-2019 (OSI-401-b)")

  plt.legend(fontsize = 8)
  plt.figure("fig2")
  tt = np.arange(0.0, len(series[t1:t2]), 1)
  daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS 2017-2018 (OSI-401-b)")
  plt.legend(fontsize = 8)
  plt.figure("fig3")
  plt.plot((np.mean(series), np.mean(series)), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS 2017-2018 (OSI-401-b)")
  plt.legend(fontsize = 8)
  
plt.figure("fig1")
plt.title(period_name + " total Antarctic sea ice area")
#plt.xticks([1, 5, 10, 15, 20, 25, 28], ["01", "05", "10", "15", "20", "25", "28"])
#plt.xlabel("Day in February 2018")
plt.ylim(0.0, 13)

import matplotlib.dates as mdates
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.legend(fontsize = 'small')
plt.ylabel("10$^6$ km$^2$")
plt.grid()
plt.tight_layout()
plt.savefig("../figs/fig1.png", dpi = 500)
print("Figure ../figs/fig1.png printed")
  
  
plt.figure("fig2")
plt.legend()
plt.title("When does the minimum of Antarctic sea ice area occur?")
plt.legend(fontsize = 'small')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.grid()
plt.yticks([],[])
plt.tight_layout()
plt.savefig("../figs/fig2.png", dpi = 500)
print("Figure ../figs/fig2.png printed")
  
  
plt.figure("fig3")
plt.title("Monthly mean Antarctic sea ice area")
plt.xlabel("10$^6$ km$^2$")
plt.yticks([],[])
plt.tight_layout()
plt.savefig("../figs/fig3.png", dpi = 500)
print("Figure ../figs/fig3.png printed")
  
