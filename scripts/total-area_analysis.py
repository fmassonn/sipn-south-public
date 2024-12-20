# Author: F. Massonnet
# Date  : December 2017
#         December 2018: update to account for new data arrival
#                        each year
#         February 2019: update to plot uncertainty better
#         December 2019: update for 2019-2020 forecast
#         April    2020: update to read data once and store it in arrays


# Purpose: large-scale sea ice diagnostics from SIPN South data

# Data: SIPN South contributors

# Imports and clean-up
# --------------------
import pandas            as pd
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

from colInterpolatOr import *

plt.close("all")

# Script parameters
# -----------------

# label with the year investigated (2017-2018, 2018-2019, ...)
myyear = "2023-2024"   

# Add obs as reference or not (False if forecast mode)
plotobs = True

# Are we after the period to be forecasted? (to know if need to plot verif)
postseason = True

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
elif myyear == "TOP2022":
  # Initialization date
  inidate = myyear[-4:] + "0501"
  # Number of days in the forecast period
  ndays   = 123
  # Label for period that is forecasted
  period_name = "May-Jun-Jul-Aug " + myyear[-4:] 
  # Starting and ending time indices (Python conventions)
  t1, t2 = 93 - 1, 93 - 1 + 31
  target_period_name = "August"
else:
  # Initialization date
  inidate = myyear[:4] + "1201"
  # Number of days in the forecast period
  ndays   = 90
  enddate = str(int(myyear[:4]) + 1) + "0228"
  # Label for period that is forecasted
  period_name = "Dec-Jan-Feb " + myyear[:4] + "-" + myyear[5:]
  # Starting and ending time indices (Python conventions)  
  t1, t2 = 63 - 1, 63 - 1 + 28
  target_period_name = "February"


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

# Create colors
sourceColors = ["#000075", "#645ccc", "#4363d8", "#1898e0", "#00b2ed", "#00bb62", \
               "#8bcd45", "#dbe622", "#f9c410", \
               "#f89e13", "#fb4c27", "#fb4865", \
               "#d24493", "#8f57bf", "#911eb4"]
col = colInterpolatOr(sourceColors, nTarget = n_sub, colorCode = "HEX")

# Change clim's color to black
for j_sub in range(n_sub):
  if sub_id[j_sub] == "climatology":
    col[j_sub] = "black"
    orderSub = 1000

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
               str(j_for).zfill(3) + "_" + inidate + "-" + enddate + "_total-area.txt"
      # Read the CSV file
      csv = pd.read_csv(filein, header = None)
      series = csv.iloc[0][:]
      # Append that series to the contribution data
      if len(series) != nt:
        print(sub_id[j_sub] + " --> WARNING: INPUT SERIES TOO LONG, CROPPING")
        series = series[:nt]
      sub_data[:, j_for - 1] = series
      
      del series, csv
    
    # Store numpy array
    data.append(sub_data)
    # Delete temporary data for that submission
    del sub_data

# Add the multi-model ensemble forecast
mmef = np.array([np.mean(d, axis=1) for j, d in enumerate(data) if sub_id[j] != "climatology"])



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
        
# Figure 1: daily ensemble mean and range
fig, ax = plt.subplots(figsize = (6, 4), dpi = dpi)

for j_sub in range(n_sub):
    mean = np.nanmean(data[j_sub], axis = 1)

    if info[j_sub][0] == "climatology":
      zorder = 1000
    else:
      zorder = 0

    plt.plot(time, mean, color = col[j_sub], lw = 1.0, 
             label = info[j_sub][0] + " " + info[j_sub][2], zorder = zorder)
    
    # Plot range as shading
    mymax = np.nanmax(data[j_sub], axis = 1)
    mymin = np.nanmin(data[j_sub], axis = 1)
    
    if info[j_sub][0] == "climatology":
	    plt.fill_between(time, mymin, mymax, 
	                     color = col[j_sub], 
	                     alpha = 0.2, lw = 0)

# Plot model ensemble median
plt.plot(time, np.median(mmef, axis = 0), color = "blue", linestyle = "--",lw = 2, label = "group forecast")    

# Plot observations if required
if plotobs:
    for j_obs, obsname in enumerate(obs):
        plt.plot(time, data_obs[j_obs], color = [0.1, 0.1, 0.1], lw = 1.5, \
                 linestyle = lst[j_obs], label = "OBS " + obsname)
# Figure polishing
plt.title(period_name + " Antarctic sea ice area")
plt.xticks([time[j] for j in [0, 14, 31, 45, 62, 76, 89]])
plt.ylim(0.0, 14)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.legend(fontsize = 6, ncol = 2)
plt.ylabel("10$^6$ km$^2$")
plt.grid()
plt.tight_layout()

for fmt in ["png" , "pdf"]:
    plt.savefig("../figs/fig1." + fmt, dpi = dpi)
    print("Figure ../figs/fig1." + fmt + " printed")


plt.close(fig)


# Figure 2: monthly means 
fig, ax = plt.subplots(figsize = (6, 6), dpi = dpi)

for j_sub in range(n_sub):
    print(sub_id[j_sub])
    if np.sum(np.isnan(data[j_sub][t1:t2, :])) > 5:
       print("STOP: too many Nans")
       print(data[j_sub][t1:t2, :])
       stop()
    monmean = np.nanmean(data[j_sub][t1:t2, :], axis = 0)
    
    ax.scatter(monmean, np.full(n_for[j_sub], n_sub - j_sub), 
               15, color = col[j_sub], 
               label = info[j_sub][0] + " " + info[j_sub][2],
               edgecolor = "white", lw = 0.2)
        
    # Plot associated PDF
    scale = 1.0
    if n_for[j_sub] >= 3:
        
        xpdf = np.linspace(0, 10, 1000)
        
        kernel = stats.gaussian_kde(monmean)
        
        if sub_id[j_sub] == "NicoSun": # He provides the range
            pdf = 1.0 * (xpdf > np.min(monmean)) * (xpdf < np.max(monmean)) 
        else:
            pdf = kernel(xpdf).T
        
        ax.fill_between(xpdf, n_sub - j_sub , n_sub - j_sub + 0.5 * pdf * scale, 
                 color = col[j_sub], alpha = 0.2, lw = 0)

# Plot mmef
xpdf = np.linspace(0, 10, 1000)
mmef_monmean = np.nanmean(mmef[:, t1:t2], axis = 1)
kernel = stats.gaussian_kde(mmef_monmean)
pdf = kernel(xpdf).T
ax.scatter(mmef_monmean, np.full(len(mmef_monmean), 0), 15, color = "blue", \
           label = "model ensemble", lw = 0.2, edgecolor = "white")
ax.fill_between(xpdf, 0 , 0 + 0.5 * pdf * scale,
                 color = "blue", alpha = 0.4, lw = 0)

# Plot observations if required
if plotobs and postseason:
    monmean_obs = list()
    
    for j_obs, obsname in enumerate(obs):
        series = data_obs[j_obs]
        mean_tmp= np.mean(series[t1:t2])
        monmean_obs.append(mean_tmp)
        
        ax.plot((mean_tmp, mean_tmp), (-1e9, 1e9), color = [0.1, 0.1, 0.1], lw = 1.5, \
                 linestyle = lst[j_obs], label = "OBS " + obsname,
                 zorder = -100)

        del mean_tmp
        
# Figure polishing
ax.set_axisbelow(True)
ax.set_title(str(myyear[5:]) + " " + target_period_name + " mean Antarctic sea ice area")
ax.legend(loc = "upper center", ncol = 4, fontsize = 6)
ax.set_xlabel("Million km$^2$")
ax.set_xlim(0, 4)
ax.set_ylim(-1.0, n_sub + 6)
ax.grid()
ax.set_yticks([])
plt.tight_layout()
for fmt in ["png", "pdf"]:
    plt.savefig("../figs/fig2." + fmt, dpi = dpi)
    print("Figure ../figs/fig2." + fmt + " printed")  
