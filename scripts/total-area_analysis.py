# Author: F. Massonnet
# Date  : December 2017
#         December 2018: update to account for new data arrival
#                        each year
#         February 2019: update to plot uncertainty better
#         December 2019: update for 2019-2020 forecast
#         April    2020: update to read data once and store it in arrays

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import pandas            as pd
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
import numpy             as np
import os
import random

from datetime import datetime
from datetime import timedelta
from netCDF4  import Dataset
from datetime import date
from scipy    import stats

plt.close("all")

# Script parameters
# -----------------

# label with the year investigated (2017-2018, 2018-2019, ...)
myyear = "2019-2020"   

# Add obs as reference or not (False if forecast mode)
plotobs = True   
# Name of observational products      
obs = ["NSIDC-0081", "OSI-401-b"]
# line styles to be used
lst = ["--",         ":"]

# Resolution for saving figures
dpi = 300

# End script parameters
# ---------------------


# Read or create meta-data
# ------------------------

if myyear == "2017-2018":
  # Initialization date
  inidate = "20180201"
  # Number of days in the forecast period
  ndays   = 28 
  # Label for period that is forecasted
  period_name = "February 2018"
  # Starting and ending time indices (Python conventions)
  t1, t2 = 0, 0 + ndays

elif myyear == "2018-2019":
  # Initialization date
  inidate = "20181201"
  # Number of days in the forecast period
  ndays   = 90
  # Label for period that is forecasted
  period_name = "Dec-Jan-Feb 2018-2019"
  # Starting and ending time indices (Python conventions)  
  t1, t2 = 63 - 1, 63 - 1 + 28    
                                  
elif myyear == "2019-2020":
  # Initialization date
  inidate = "20191201"
  # Number of days in the forecast period
  ndays   = 90
  # Label for period that is forecasted
  period_name = "Dec-Jan-Feb 2019-2020"
  # Starting and ending time indices (Python conventions)  
  t1, t2 = 63 - 1, 63 - 1 + 28


# Load namelist
exec(open("./namelist_" + myyear + ".py").read())

# Create time axis
#time = pd.date_range(pd.to_datetime(inidate, format = "%Y%m%d"), 
#                     periods = ndays).tolist()
time = [datetime.strptime(inidate, "%Y%m%d") + timedelta(days = x) for x in 
                    range(0, ndays)]
nt = len(time)

# Number of submissions
n_sub = len(info)

# Submission ID
sub_id = [info[j_sub][0] for j_sub in range(n_sub)]

# Range of forecasts for each submission
range_for = [info[j_sub][1] for j_sub in range(n_sub)]

# Number of forecasts for each submission
n_for    = [len(l) for l in range_for]

# Colors for each submission
col = [info[j_sub][2] for j_sub in range(n_sub)]


# Store the raw data
# ------------------
# The daily total areas will be stored in a list. Each item of the list
# will correspond to one submission and will be a 2-D numpy array
# of dimensions {time, number of ensemble forecasts}
data = list()

for j_sub in range(n_sub):
    print("Reading " + sub_id[j_sub])
    # Create empty numpy array
    sub_data = np.full((nt, n_for[j_sub]), np.nan)
  
    # Read in data for that submission
    # Note that j_for refers to the forecast index in non-Pythonic
    # conventions
    for j_for in range_for[j_sub]:
      filein = "../data/" + myyear + "/txt/" + sub_id[j_sub] + "_" + \
               str(j_for).zfill(3) + "_total-area.txt"
      # Read the CSV file
      csv = pd.read_csv(filein, header = None)
      series = csv.iloc[0][:]
      # Append that series to the contribution data
      sub_data[:, j_for - 1] = series
      del series, csv
    
    # Store numpy array
    data.append(sub_data)
    # Delete temporary data for that submission
    del sub_data

# Repeat with observations, if needed
if plotobs:
    # A list of 1-D numpy arrays, each of dimensions {time}
    data_obs= list()
    
    for obsname in obs:
        filein = "../data/" + myyear + "/txt/" + obsname + \
        "_000_total-area.txt"
        # Read the CSV file
        csv = pd.read_csv(filein, header = None)
        series = csv.iloc[0][:]
        # It might be that the obs is not as long as the forecast, 
        # when the verification period is not over yet
        # --> complete the data with NaNs
        obscomplete = (nt == np.sum((1 - np.isnan(series)) * 1))
        if not obscomplete:
            series = series.append(pd.Series([np.nan for i in \
                     range(nt - len(series))]), ignore_index = True)
    
        data_obs.append(series)
        
        del series, csv, obscomplete
        
# Process the raw data to produce diagnostics
# and plot them
# -------------------------------------------  
        
# Figure 1: daily ensemble median and range
fig, ax = plt.subplots(figsize = (6, 4), dpi = dpi)

for j_sub in range(n_sub):
    median = np.median(data[j_sub], axis = 1)

    plt.plot(time, median, color = col[j_sub], lw = 1.5, 
             label = info[j_sub][0] + " " + info[j_sub][3])
    
    # Plot range as shading
    mymax = np.max(data[j_sub], axis = 1)
    mymin = np.min(data[j_sub], axis = 1)
    plt.fill_between(time, mymin, mymax, 
                     color = [c * 1.0 for c in col[j_sub]], 
                     alpha = 0.2, lw = 0)
    
# Plot observations if required
if plotobs:
    for j_obs, obsname in enumerate(obs):
        plt.plot(time, data_obs[j_obs], color = [0.1, 0.1, 0.1], lw = 1.5, \
                 linestyle = lst[j_obs], label = "OBS " + myyear + \
                 "\n" + obsname)
# Figure polishing
plt.title(period_name + " total Antarctic sea ice area")
plt.xticks([time[j] for j in [0, 14, 31, 45, 62, 76, 89]])
plt.ylim(0.0, 14)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.legend(fontsize = 7)
plt.ylabel("10$^6$ km$^2$")
plt.grid()
plt.tight_layout()

for fmt in ["png", "eps", "pdf"]:
    plt.savefig("../figs/fig1." + fmt)
    print("Figure ../figs/fig1." + fmt + " printed")


plt.close(fig)

# Figure 2: minimum of time series
fig, ax = plt.subplots(figsize = (5, 6), dpi = dpi)

for j_sub in range(n_sub):
    # List that will have the day of minimum for each member
    daymin = list()
    for j_for in range_for[j_sub]:

        # Pythonic convention requires to start one index earlier
        series = data[j_sub][:, j_for - 1]
        # Create dummy time index to locate the minimum. 
        # This index is expressed as seconds since first date analyzed
    
        #tt = np.arange(time[t1].timestamp(), time[t2 - 1].timestamp())
        tt = np.arange(t1, t2, 1 / 24)
        # The raw data is fitted by a quadratic polynomial,
        # and the day at which the minimum is achieved is recorded
        
        coeffs = np.polyfit(np.arange(t1, t2), series[t1:t2], 2)
        # Minimum of ax^2 + bx + c occurs at  - b / 2a
        daymin.append(time[0] + timedelta(days = -coeffs[1] / (2 * coeffs[0])))

        del series, coeffs
        
    # Plot all days of minimum
    ax.scatter(daymin, np.full(n_for[j_sub], n_sub - j_sub), 
               15, color = col[j_sub], 
               label = info[j_sub][0] + " " + info[j_sub][3],
               edgecolor = "white", lw = 0.2)
        
    # Plot associated PDF
    scale = 5e5
    if n_for[j_sub] >= 3:
        daymin_sec = np.array([d.timestamp() for d in daymin])
        #xpdf = np.linspace(time[0].timestamp(), time[-1].timestamp() + , 10000)#
        xpdf = np.linspace(np.min(daymin_sec) - 10 * 86400, \
                           np.max(daymin_sec) + 10 * 86400, 10000)
        kernel = stats.gaussian_kde(daymin_sec)
        pdf = kernel(xpdf).T
        ax.fill_between([datetime.fromtimestamp(x) for x in xpdf],  \
                 np.full(len(pdf), n_sub - j_sub),
                 n_sub - j_sub + scale * pdf, 
                 color = col[j_sub], alpha = 0.2, lw = 0)

#  
# Figure polishing
plt.title("When does the minimum of Antarctic sea ice area occur?")
plt.legend(fontsize = "small")
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.xticks([time[j] for j in [62, 71, 81, 89]])
plt.xlim(time[t1] - timedelta(days = 10), time[t2 - 1] + timedelta(days = 20))

plt.grid()
plt.yticks([],[])
plt.tight_layout()
for fmt in ["png", "eps", "pdf"]:
    plt.savefig("../figs/fig2." + fmt)
    print("Figure ../figs/fig2." + fmt + " printed")

stop()
## Record when the minimum of the series occurs
#    # --------------------------------------------
#    # Fit a quadratic polynomial to filter out weather variability
#daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
#
#plt.figure("fig2")
#if info[j_sub][3] != "(monthly)":
#  plt.scatter(daymin, n_sub - j_sub, 50 / n_for[j_sub], color = col[j_sub], 
#               label = mylab, marker = ".")
#  #plt.text(46, n_sub - j_sub, info[j_sub][0] + " " + info[j_sub][3], color = col[j_sub], ha = "right")
#
#plt.figure("fig3")
#plt.scatter(np.mean(series[t1:t2]), n_sub - j_sub, 50 / n_for[j_sub], 
#            color = col[j_sub], 
#            label = mylab)
#
#
#    
#
#
#
#
#
#
## Monthly means
#plt.figure("fig3", figsize = (6, 5))
#
#for j_sub in range(n_sub):
#  print("Processing " + sub_id[j_sub])
#  # Will have list of forecasted areas
#  submission = list() 
#  for j_for in range_for[j_sub]:
#    # Total area
#    filein = "../data/" + myyear + "/txt/" + sub_id[j_sub] + "_" + str(j_for).zfill(3) + "_total-area.txt"
#    # Read the CSV file
#    csv = pd.read_csv(filein, header = None)
#    series = csv.iloc[0][:]
#    # Append that series to the contribution data
#    submission.append(series)
#
#    # Plot series, line thinner for large ensembles. Legend only if first member
#    if j_for == range_for[j_sub][0]:
#      mylab = info[j_sub][0] + " " + info[j_sub][3]
#    else:
#      mylab = "_nolegend_"
#
#    #plt.figure("fig1")
#    #plt.plot(time, series, color = col[j_sub], lw = (5.0 / n_for[j_sub]) ** 0.2, label = mylab)
#
#
#    
#  # Plot PDF

#
#  # Plot ensemble mean
#  mean = np.mean(np.array(submission), axis = 0)
#  plt.figure("fig1")
#  plt.plot(time, mean, color = col[j_sub], lw = 1.5, label = info[j_sub][0] + " " + info[j_sub][3])
#  # Plot range as shading
#  mymax = np.max(np.array(submission), axis = 0)
#  mymin = np.min(np.array(submission), axis = 0)
#  plt.fill_between(time, mymin, mymax, color = [c * 1.0 for c in col[j_sub]], alpha = 0.5, lw = 0)
#    
#
## Plot observations
## -----------------
#if plotobs:
#
#  filein = "../data/" + myyear + "/txt/NSIDC-0081_000_total-area.txt"
#  csv = pd.read_csv(filein, header = None)
#  series = csv.iloc[0][:]
#  # It can be that the obs is not as long as the forecast, when the verification period is not over yet
#  # --> complete with NaN.
#  obscomplete = (len(time) == np.sum((1 - np.isnan(series)) * 1))
#  if not obscomplete:
#      series = series.append(pd.Series([np.nan for i in range(len(time) - len(series))]), ignore_index = True)
#  plt.figure("fig1")
#  plt.plot(time, series, color = [0.1, 0.1, 0.1], lw = 1.5, linestyle = "--", label = "OBS " + myyear + "\n(NSIDC-0081)")
#
#  plt.figure("fig2")
#  tt = np.arange(0.0, len(series[t1:t2]), 1)
#  # Don't calculate minimum if series not yet completed
#  if obscomplete:
#      daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
#  else:
#      daymin = np.nan
#  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS " + myyear + "\n(NSIDC-0081)")
#  plt.figure("fig3")
#  plt.plot((np.mean(series[t1:t2]), np.mean(series[t1:t2])), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = "--", label = "OBS " + myyear + " (NSIDC-0081)")
#  
#  
#  
#  filein = "../data/" + myyear + "/txt/OSI-401-b_000_total-area.txt"
#  csv = pd.read_csv(filein, header = None)
#  series = csv.iloc[0][:]
#  # It can be that the obs is not as long as the forecast, when the verification period is not over yet
#  # --> complete with NaN.
#  obscomplete = (len(time) == np.sum((1 - np.isnan(series)) * 1))
#  if not obscomplete:
#      series = series.append(pd.Series([np.nan for i in range(len(time) - len(series))]), ignore_index = True)
#  plt.figure("fig1")
#  plt.plot(time, series, color = [0.1, 0.1, 0.1], lw = 3.0, linestyle = ":", label = "OBS " + myyear + "\n(OSI-401-b)")
#
#  plt.legend(fontsize = 8)
#  plt.figure("fig2")
#  tt = np.arange(0.0, len(series[t1:t2]), 1)
#  if obscomplete:
#      daymin = time[t1:t2][np.argmin(np.polyval(np.polyfit(range(len(series[t1:t2])), series[t1:t2], 2), tt))]
#  else:
#      daymin = np.nan
#
#  plt.plot((daymin, daymin), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS " + myyear + "\n(OSI-401-b)")
#  plt.legend(fontsize = 8)
#  plt.figure("fig3")
#  plt.plot((np.mean(series[t1:t2]), np.mean(series[t1:t2])), (0, n_sub), color = [0.1, 0.1, 0.1], linestyle = ":", label = "OBS " + myyear + " (OSI-401-b)")
#  plt.legend(fontsize = 8)
#  
#plt.figure("fig1")
#plt.title(period_name + " total Antarctic sea ice area")
#
##plt.xticks([1, 5, 10, 15, 20, 25, 28], ["01", "05", "10", "15", "20", "25", "28"])
##plt.xlabel("Day in February 2020")
#plt.xticks([time[j] for j in np.arange(0, len(time), 15)])
#plt.ylim(0.0, 13)
#
#import matplotlib.dates as mdates
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
#plt.legend(fontsize = 7)
#plt.ylabel("10$^6$ km$^2$")
#plt.grid()
#plt.tight_layout()
#plt.savefig("../figs/fig1.png", dpi = 500)
#print("Figure ../figs/fig1.png printed")
#plt.savefig("../figs/fig1.eps")
#plt.savefig("../figs/fig1.pdf")
#
#  
#plt.figure("fig3")
#plt.title("Monthly mean Antarctic sea ice area")
#plt.xlabel("10$^6$ km$^2$")
#plt.yticks([],[])
#plt.tight_layout()
#plt.savefig("../figs/fig3.png", dpi = 500)
#print("Figure ../figs/fig3.png printed")
#plt.savefig("../figs/fig3.pdf")
#plt.savefig("../figs/fig3.eps")
#  
