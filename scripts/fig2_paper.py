# Author: F. Massonnet
# Date  : 10 Jan 2022


# Purpose: Figure with time series of SIPN South contributions

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

from seaice_commondiags import *


plt.close("all")


# Open sector information file
# ----------------------------
exec(open("./namelist.py").read())
nSectors = len(sectors)
# Script parameters
# -----------------

# Years for the long-term record and the climatology
yearbClim = 1979
yeareClim = 2015

targetMonth = 1 # February (Python convention)

fileObsLong = "/Users/massonnetf/CLIMDATA/obs/ice/siconc/OSI-SAF" + \
              "/OSI-450/processed/native/siconc_SImon_OSI-450_r1i1p1_"   + \
              str(yearbClim) + "01-" + str(yeareClim) + "12_sh.nc"

yearsClim = np.arange(yearbClim, yeareClim  + 1)



# Open file
# ---------
f = Dataset(fileObsLong)
sic = f.variables["siconc"][:]
lat = f.variables["lat"][:]
lon = f.variables["lon"][:]
lon[lon > 180] = lon[lon > 180] - 360
f.close()

cellarea = 25000 ** 2

# Mean sea ice area for the target month

sia = compute_area(sic[targetMonth::12], cellarea, lat < 0.0)



fig, ax = plt.subplots(3, 2, dpi = 300, figsize = (8, 6) )


for s, a in zip(sectors, ax.flatten()):
    print(s[0])
    lonW, lonE = s[1], s[2]
    
    # Create mask
    if lonE > lonW:
        mask = (lat < 0.0) * (lon > lonW) * (lon <= lonE)
    else:
        mask = (lat < 0.0) * ( (lon > lonW) + (lon <= lonE))
    # Compute area
    sia_sector = compute_area(sic[targetMonth::12], cellarea, mask)
    
    # Fit a PDF
    x_pdf = np.linspace(0, 5, 1000)
    kernel0 = stats.gaussian_kde(sia_sector)
    pdf = kernel0(x_pdf).T
    
    
    a.plot(yearsClim, sia_sector)#label = "OSI-450 (OBS)")
    a.plot(2015 + pdf, x_pdf)
    
    a.set_title(s[0])

    a.set_xlim(1978, 2023) 
    a.set_ylim(bottom = 0, top = np.max(sia_sector) * 1.2)
    a.grid()


# Add the two verification products
# ---------------------------------


# Fetch starting dates and their number
startDates  = [n[0] for n in namelistOutlooks]
endDates    = [n[1] for n in namelistOutlooks]

nStartDates = len(startDates)

# The list of years defining the start dates
startYears   = [s.year for s in startDates]
endYears     = [e.year for e in endDates]

# The name of seasons (ex: 2017-2018)
nameSeasons = [str(s) + "-" + str(e) for s, e in zip(startYears, endYears)]

# Number of identified contributors to SIPN South
nContributors = len(namelistContributions)



# 
obsVerif = ["NSIDC-0081", "OSI-401-b"]

for s, season in enumerate(nameSeasons):
    
    for obsname in obsVerif:
        filein = "../data/" + season + "/txt/" + obsname + \
         "_000_total-area.txt"

        # Read the CSV file
        csv = pd.read_csv(filein, header = None)
        series = csv.iloc[0][:]
        
        # if more than 2 days missing, do nothing
        
        if np.sum(np.isnan(series[-28:])) > 2:
            pass
        else:
            obsMean = np.nanmean(series[-28:])


        # Plot it
        ax[2, 1].scatter(endYears[s], obsMean, 5, marker = "s", color = "black")
        

# # Complete with the SIPN South forecast data
# # ----------------------------------------------


data = [list() for _ in range(nStartDates)]
# This list has as many items as there are start dates / outlooks

# Now we populate each of the item by going from one outlook to the next
for jStartDates in range(nStartDates):
    
    # We create a new list with as many items as there are diagnostics
    listTmp = [list() for _ in range()]
    
    


# # The idea is to populate a 2-item list. The first item is for total area
# # and the second is for regional area. The first item is an array with dim-
# # ensions: nSub x nFor x nDay
# # where
# nSub = len(info) # Number of groups contributing
# nFor = np.max([i[1:] for i in info]) # Max number of ensemble forecasts
# nDay = 90 # number of days in forecasting season
# nYear= yearEnd - yearStart + 1 # Number of hindcasts
# years = np.arange(yearStart, yearEnd + 1)

# itemTotalArea =    np.full((nSub, nFor, nYear, nDay), np.nan)

# for jSub, i in enumerate(info):
#     for jYear, y in enumerate(years):
        
#         # Total Area
#         if i[jYear + 1][0] == 0:
#             # If there is no forecast of total area for that year
#             pass
#         else:
#             # Load the files
#             for jFor in np.arange(i[jYear + 1][0]):
#                 filein = "../data/" + str(y) + "-" + str(y + 1) + "/txt/" + i[0] + "_" + \
#                     str(jFor + 1).zfill(3) + "_total-area.txt"
#                 # Read the CSV file
#                 csv = pd.read_csv(filein, header = None)
#                 series = csv.iloc[0][:]
#                 # Append that series to the contribution data
#                 if y == 2017:
#                     itemTotalArea[jSub, jFor, jYear, -28:] = series[:28]
#                 else:
#                     itemTotalArea[jSub, jFor, jYear, :] = series[:90]

# # Model mean, Ensemble mean
# grandMean = np.mean(np.nanmean(itemTotalArea, axis = (0, 1))[:, -28:], axis = 1)
# ax[2, 1].plot(years, grandMean, "-s")

#fig.tight_layout()

for a in ax.flatten():
    a.set_xlim(2015, 2023)
    a.set_axisbelow(True)

fig.legend()
fig.tight_layout()

fig.savefig("../figs/fig2_paper.png", dpi = 300)




stop()




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

else:
  # Initialization date
  inidate = myyear[:4] + "1201"
  # Number of days in the forecast period
  ndays   = 90
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

# Colors for each submission
col = [info[j_sub][2] for j_sub in range(n_sub)]

# Convert to RGB if necessary
for j_sub in range(n_sub):
    if type(col[j_sub]) is str:
        col[j_sub] = matplotlib.colors.to_rgb(col[j_sub])

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
      if len(series) != nt:
        print("WARNING: INPUT SERIES TOO LONG, CROPPING")
        series = series[:nt]
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
                 linestyle = lst[j_obs], label = "OBS " + obsname)
# Figure polishing
plt.title(period_name + " total Antarctic sea ice area")
plt.xticks([time[j] for j in [0, 14, 31, 45, 62, 76, 89]])
plt.ylim(0.0, 14)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
plt.legend(fontsize = 7)
plt.ylabel("10$^6$ km$^2$")
plt.grid()
plt.tight_layout()

for fmt in ["png" , "pdf"]:
    plt.savefig("../figs/fig1." + fmt, dpi = dpi)
    print("Figure ../figs/fig1." + fmt + " printed")


plt.close(fig)

# Figure 2: minimum of time series
fig, ax = plt.subplots(figsize = (6, 4), dpi = dpi)

for j_sub in range(n_sub):
    # List that will have the day of minimum for each member
    daymin = list()
    for j_for in range_for[j_sub]:

        # Pythonic convention requires to start one index earlier
        series = data[j_sub][:, j_for - 1]
        # Create dummy time index to locate the minimum. 
        # This index is expressed as seconds since first date analyzed
    
        # Create dummy time axis (hourly resolution) to locate the minimum
        # after the quadratic fit has been performed
        tt = np.arange(t1, t2, 1 / 24)
        # The raw data is fitted by a quadratic polynomial,
        # and the day at which the minimum is achieved is recorded
        
        coeffs = np.polyfit(np.arange(t1, t2), series[t1:t2], 2)
        
        if np.max(np.abs(coeffs)) == 0.0:
            # In this case likely the model goes to zero --> set first date of 
            # zero
            daymin.append(time[0] + timedelta(days = float(np.where(series == 0)[0][0])))
        else:
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
        xpdf = np.linspace(np.min(daymin_sec) - 10 * 86400, \
                           np.max(daymin_sec) + 10 * 86400, 10000)
        kernel = stats.gaussian_kde(daymin_sec)
        pdf = kernel(xpdf).T
        ax.fill_between([datetime.fromtimestamp(x) for x in xpdf],  \
                 np.full(len(pdf), n_sub - j_sub),
                 n_sub - j_sub + scale * pdf, 
                 color = col[j_sub], alpha = 0.2, lw = 0)

# Plot grey area right to end of forecasting period
ax.fill_between((time[-1], time[-1] + timedelta(days= 31), time[-1] + 
                 timedelta(days= 31), time[-1]),
                (-1e9, -1e9, 1e9, 1e9), color = [0.9, 0.9, 0.9], 
                zorder = -1000)

# Plot observations if required and if complete
if plotobs and postseason:
    daymin_obs = list()
    
    for j_obs, obsname in enumerate(obs):
        series = data_obs[j_obs]
        # Create dummy time axis (hourly resolution) to locate the minimum
        # after the quadratic fit has been performed
        tt = np.arange(t1, t2, 1 / 24)
        # The raw data is fitted by a quadratic polynomial,
        # and the day at which the minimum is achieved is recorded
        myTime   = np.arange(t1, t2)
        mySeries = series[t1:t2]


        # It happens (in 2021 at least) that there are NaNs because
        # missing data, which need to be taken out then
        myTime = myTime[~np.isnan(mySeries)]
        mySeries = mySeries[~np.isnan(mySeries)]    

        coeffs = np.polyfit(myTime, mySeries, 2)
        
        # Minimum of ax^2 + bx + c occurs at  - b / 2a
        daymin_obs_tmp = time[0] + \
                         timedelta(days = -coeffs[1] / (2 * coeffs[0]))
        daymin_obs.append(daymin_obs_tmp)
        
        ax.plot((daymin_obs_tmp, daymin_obs_tmp), (-1e9, 1e9), 
                color = [0.1, 0.1, 0.1], lw = 1.5, \
                 linestyle = lst[j_obs], label = "OBS " + obsname, 
                 zorder = -100)
        
        del daymin_obs_tmp
    del daymin_obs

# Figure polishing
ax.set_axisbelow(True)
ax.set_title("When will the minimum of "  + myyear[5:]+ " Antarctic sea ice area occur?")
ax.legend(loc = "upper left", ncol = 2, fontsize = 7)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
ax.set_xticks([time[j] for j in [62, 71, 81, 89]])
ax.set_xlim(time[t1] - timedelta(days = 10), time[t2 - 1] + \
            timedelta(days = 20))
ax.set_ylim(0.0, n_sub + 1)
ax.grid()
ax.set_yticks([],[])
plt.tight_layout()
for fmt in ["png", "pdf"]:
    plt.savefig("../figs/fig2." + fmt, dpi = dpi)
    print("Figure ../figs/fig2." + fmt + " printed")




# Figure 3: monthly means 
fig, ax = plt.subplots(figsize = (5, 5), dpi = dpi)

for j_sub in range(n_sub):
    # List that will have the day of minimum for each member
    monmean = np.mean(data[j_sub][t1:t2, :], axis = 0)
    
    ax.scatter(monmean, np.full(n_for[j_sub], n_sub - j_sub), 
               15, color = col[j_sub], 
               label = info[j_sub][0] + " " + info[j_sub][3],
               edgecolor = "white", lw = 0.2)
        
    # Plot associated PDF
    scale = 0.5
    if n_for[j_sub] >= 3:
        
        xpdf = np.linspace(0, 10, 1000)
        
        kernel = stats.gaussian_kde(monmean)
        
        pdf = kernel(xpdf).T
        
        ax.fill_between(xpdf, n_sub - j_sub , n_sub - j_sub + 0.5 * pdf, 
                 color = col[j_sub], alpha = 0.2, lw = 0)

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
ax.set_title(target_period_name + " mean sea ice area")
ax.legend(loc = "upper center", ncol = 4, fontsize = 7)
ax.set_xlabel("Million km$^2$")
ax.set_xlim(0, 4)
ax.set_ylim(0.0, n_sub + 6)
ax.grid()
ax.set_yticks([],[])
plt.tight_layout()
for fmt in ["png", "pdf"]:
    plt.savefig("../figs/fig3." + fmt, dpi = dpi)
    print("Figure ../figs/fig3." + fmt + " printed")  
