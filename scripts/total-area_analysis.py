#!/usr/bin/python
# 
# Author: F. Massonnet
# Date  : December 2017
#         December 2018: update to account for new data arrival
#                        each year
#         February 2019: update to plot uncertainty better

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import pandas            as pd
import matplotlib.pyplot as plt
import numpy             as np
import os
import random

from datetime import datetime
from datetime import timedelta
from netCDF4  import Dataset
from datetime import date

plt.close("all")

# Script parameters
myyear = "2018-2019"  # label with the year investigated (2017-2018, 2018-2019, ...)
plotobs = True        # Add obs as reference or not (False if forecast mode)

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

# Load namelist
exec(open("./namelist_" + myyear + ".py").read())

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
plt.figure("fig2", figsize = (6, 5))
# Monthly means
plt.figure("fig3", figsize = (6, 5))

for j_sub in range(n_sub):
  print("Processing " + sub_id[j_sub])
  # Will have list of forecasted areas
  submission = list() 
  for j_for in list_for[j_sub]:
    # Total area
    filein = "../data/" + myyear + "/txt/" + sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_total-area.txt"
    # Read the CSV file
    csv = pd.read_csv(filein, header = None)
    series = csv.iloc[0][:]
    # Append that series to the contribution data
    submission.append(series)

    # Plot series, line thinner for large ensembles. Legend only if first member
    if j_for == list_for[j_sub][0]:
      mylab = info[j_sub][0] + " " + info[j_sub][3]
    else:
      mylab = "_nolegend_"

    #plt.figure("fig1")
    #plt.plot(time, series, color = col[j_sub], lw = (5.0 / n_for[j_sub]) ** 0.2, label = mylab)


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

  # Plot ensemble mean
  mean = np.mean(np.array(submission), axis = 0)
  plt.figure("fig1")
  plt.plot(time, mean, color = col[j_sub], lw = 1.5, label = info[j_sub][0] + " " + info[j_sub][3])
  # Plot range as shading
  mymax = np.max(np.array(submission), axis = 0)
  mymin = np.min(np.array(submission), axis = 0)
  plt.fill_between(time, mymin, mymax, color = [c * 1.0 for c in col[j_sub]], alpha = 0.5, lw = 0)
    

# Plot observations
# -----------------
if plotobs:

  filein = "../data/2018-2019/txt/NSIDC-0081_000_total-area.txt"
  csv = pd.read_csv(filein, header = None)
  series = csv.iloc[0][:]
  # It can be that the obs is not as long as the forecast, when the verification period is not over yet
  # --> complete with NaN.
  obscomplete = (len(time) == np.sum((1 - np.isnan(series)) * 1))
  if not obscomplete:
      series = series.append(pd.Series([np.nan for i in range(len(time) - len(series))]), ignore_index = True)
  plt.figure("fig1")
  plt.plot(time, series, color = [0.1, 0.1, 0.1], lw = 1.5, linestyle = "--", label = "OBS 2018-2019 (NSIDC-0081)")

  plt.figure("fig2")
  tt = np.arange(0.0, len(series[t1:t2]), 1)
  # Don't calculate minimum if series not yet completed
  if obscomplete:
      daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
  else:
      daymin = np.nan
  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS 2018-2019 (NSIDC-0081)")
  plt.figure("fig3")
  plt.plot((np.mean(series), np.mean(series)), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS 2018-2019 (NSIDC-0081)")
  
  
  
  filein = "../data/2018-2019/txt/OSI-401-b_000_total-area.txt"
  csv = pd.read_csv(filein, header = None)
  series = csv.iloc[0][:]
  # It can be that the obs is not as long as the forecast, when the verification period is not over yet
  # --> complete with NaN.
  obscomplete = (len(time) == np.sum((1 - np.isnan(series)) * 1))
  if not obscomplete:
      series = series.append(pd.Series([np.nan for i in range(len(time) - len(series))]), ignore_index = True)
  plt.figure("fig1")
  plt.plot(time, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = ":", label = "OBS 2018-2019 (OSI-401-b)")

  plt.legend(fontsize = 8)
  plt.figure("fig2")
  tt = np.arange(0.0, len(series[t1:t2]), 1)
  if obscomplete:
      daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
  else:
      daymin = np.nan

  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS 2018-2019 (OSI-401-b)")
  plt.legend(fontsize = 8)
  plt.figure("fig3")
  plt.plot((np.mean(series), np.mean(series)), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS 2018-2019 (OSI-401-b)")
  plt.legend(fontsize = 8)
  
plt.figure("fig1")
plt.title(period_name + " total Antarctic sea ice area")

#plt.xticks([1, 5, 10, 15, 20, 25, 28], ["01", "05", "10", "15", "20", "25", "28"])
#plt.xlabel("Day in February 2018")
plt.ylim(0.0, 13)

import matplotlib.dates as mdates
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.legend(fontsize = 7)
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
  
