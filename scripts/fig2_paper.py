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
        

# Complete with the SIPN South forecast data
# ----------------------------------------------

# We are creating a list with as many elements as there are startDates
#
myData = list()

# myData : 1 / Which startDate
#          2 / Which sector
#          3 / Which contributor
#          4 / Which forecast

for i in range(nStartDates):
    print("Loading " + nameSeasons[i])
    
    iList = list()
    for j in range(nSectors):
        jList = list()
        
        for k in range(nContributors):
            nameContributor = namelistContributions[k][0]
            kList = list()
            
            
            if sectors[j][0] == "Southern Ocean": # total area
                nForecasts = namelistContributions[k][i + 1][0] # 0 since total area, i+ 1 since first item is name
            
                for l in range(nForecasts):
                                
                    # Fetch the relevant file
                    filein = "../data/" + nameSeasons[i] + "/txt/" + nameContributor + "_" + \
                     str(l + 1).zfill(3) + "_total-area.txt"
                    
                    
                    csv = pd.read_csv(filein, header = None)
                    seriesTmp = csv.iloc[0][:]
                    
                    kList.append(seriesTmp)
                    
            else:
                
                # Read regional file
                nForecasts = namelistContributions[k][i + 1][1] # 1 since regional area, i+ 1 since first item is name
                
                
                lonMin = sectors[j][1]
                lonMax = sectors[j][2]
                
                # Recast everything into 0, 360
                if lonMin < 0:
                    lonMin += 360
                if lonMax < 0:
                    lonMax +=360
                

                # Fetch relevant file
                for l in range(nForecasts):
                     filein = "../data/" + nameSeasons[i] + "/txt/" + nameContributor + "_" + \
                     str(l + 1).zfill(3) + "_regional-area.txt"
                     
                     csv = pd.read_csv(filein, header = None)
                     arrayTmp = np.array(csv.iloc[:][:])
                     
                     # Sum over relevant rows to accumulate regional average
                     # The data comes by 10° increments
                     
                     if lonMin > lonMax: # This can happen for e.g. weddell sea

                         row0 = int(lonMax / 10.0)
                         row1 = int(lonMin / 10.0)
                         seriesTmp = np.sum(arrayTmp, axis = 0) - np.sum(arrayTmp[row0:row1, :], axis = 0)
                         #if i == 1 and j == 0 and k == 12 and l ==6:
                         #    stop()
                     else:
                         row0 = int(lonMin / 10.0)
                         row1 = int(lonMax / 10.0)
                         seriesTmp = np.sum(arrayTmp[row0:row1, :], axis = 0)
                         
                     kList.append(seriesTmp)
                         
            
            jList.append(kList)
        
        iList.append(jList)
        
    myData.append(iList)
    
    



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

# =============================================================================
# for a in ax.flatten():
#     a.set_xlim(2015, 2023)
#     a.set_axisbelow(True)
# 
# fig.legend()
# fig.tight_layout()
# 
# fig.savefig("../figs/fig2_paper.png", dpi = 300)
# 
# 
# 
# 
# stop()
# =============================================================================

