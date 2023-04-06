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

matplotlib.rcParams['font.family'] = "Arial Narrow"

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



fig, ax = plt.subplots(3, 2, dpi = 300, figsize = (8, 8) )


# counter for putting letters
jLett = 0

for s, a in zip(sectors, ax.flatten()):

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
    
    
    # Rescale for viz purposes.
    pdf /= np.max(pdf) 
    
    myClim = np.median(sia_sector)
    
    #a.plot(yearsClim, sia_sector)
    #a.fill_between(2016 + pdf, x_pdf, facecolor = "#E16E79", alpha = 0.4, edgecolor = "#E16E79")

    a.scatter((2016, 2016), (myClim, myClim), 30, marker = "s", color = "#D02433")
    a.plot((2016, 2016), (np.percentile(sia_sector, 10), np.percentile(sia_sector, 90)), color = "#D02433", lw = 2)
    a.scatter(np.ones(len(sia_sector)) * 2016, sia_sector, 5, marker = ".", color = "#D02433")
    a.plot((2016 - 0.1, 2016 + 0.1), (np.percentile(sia_sector, 90), np.percentile(sia_sector, 90)), color = "#D02433", lw = 2)
    a.plot((2016 - 0.1, 2016 + 0.1), (np.percentile(sia_sector, 10), np.percentile(sia_sector, 10)), color = "#D02433", lw = 2)

    myClim = np.median(sia_sector)

    a.plot((2016, 2030), (myClim, myClim), color = "#D02433", linestyle = ":", lw = 1)

    
    a.set_title("(" + "abcdefghijkl"[jLett] + ") " + s[0])

    jLett += 1

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
        
        for j, (_, a) in enumerate(zip(sectors, ax.flatten())):
            if sectors[j][0] == "Southern Ocean":
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
                    a.scatter(endYears[s], obsMean, 15, marker = "x", color = "black", lw = 1)
                    
            else:
                filein = "../data/" + season + "/txt/" + obsname + \
                         "_000_regional-area.txt"
                csv = pd.read_csv(filein, header = None)
                arrayTmp = np.array(csv.iloc[:][:])
                
                lonMin = sectors[j][1]
                lonMax = sectors[j][2]
                
                # Recast everything into 0, 360
                if lonMin < 0:
                    lonMin += 360
                if lonMax < 0:
                    lonMax +=360
                

                     
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
  
                # Time (monthly) mean:
                t1 = (namelistOutlooks[s][2] -  namelistOutlooks[s][4]).days
                t2 = (namelistOutlooks[s][3] -  namelistOutlooks[s][4]).days + 1
                    
                # Do the time monthly-mean over target month
                if np.sum(np.isnan(seriesTmp[-28:])) > 2:
                    obsMean = np.nan
                else:
                    obsMean = np.nanmean(seriesTmp[t1:t2])
                print(sectors[j][0] + " " + season + ": " + str(obsMean))
           
                a.scatter(endYears[s], obsMean, 15, marker = "x", color = "black")



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
    
    
# Plot the grand distribution statistics. For each start date, and sector, we 
# want to end up with a distribution of monthly-mean February sea ice area
# forecasts. The number of contributors varies from one year to the next,
# so we need a list of list of arrays.

sipnAreas = list()

for i in range(nStartDates):
    iList = list()
    for j in range(nSectors):
        jList = list()
        
        
        for k in range(nContributors):
            
            thisInput = myData[i][j][k]
            
            if len(thisInput) == 0: # If no contribution for that year
                pass
            else:
                nForecasts = len(myData[i][j][k])
                arrayTmp = np.full(nForecasts, np.nan)
                for l in range(nForecasts):
                    seriesTmp = myData[i][j][k][l]
                    t1 = (namelistOutlooks[i][2] -  namelistOutlooks[i][4]).days
                    t2 = (namelistOutlooks[i][3] -  namelistOutlooks[i][4]).days + 1
                    
                    # Do the time monthly-mean over target month
                    arrayTmp[l] = np.mean(seriesTmp[t1:t2])
                    
 
                    
                # Ensemble mean
                jList.append(np.mean(arrayTmp))

        iList.append(jList)
        
    sipnAreas.append(iList)


# Plot
for i in range(nStartDates):
    for j, (_, a) in enumerate(zip(sectors, ax.flatten())):
        
        #a.scatter(endYears[i], 
        a.scatter(endYears[i] + 0.3, np.median(sipnAreas[i][j]), 30, color = "#3878DB", marker = "s")
        a.plot((endYears[i] + 0.3, endYears[i] + 0.3), (np.percentile(sipnAreas[i][j], 10), np.percentile(sipnAreas[i][j], 90)), \
               lw = 2, color = "#3878DB")
        a.plot((endYears[i] + 0.2, endYears[i]+ 0.4), (np.percentile(sipnAreas[i][j], 10), np.percentile(sipnAreas[i][j], 10)), \
                 color = "#3878DB", lw = 2)
        a.plot((endYears[i] + 0.2, endYears[i]+ 0.4), (np.percentile(sipnAreas[i][j], 90), np.percentile(sipnAreas[i][j], 90)), \
                 color = "#3878DB", lw = 2)
        a.scatter(np.ones(len(sipnAreas[i][j])) * (endYears[i] + 0.3), sipnAreas[i][j], 5, color = "#3878DB", marker = ".")
        # Plot PDF
        #box = a.boxplot( sipnAreas[i][j], positions = [endYears[i] + 0.3], patch_artist=True, \
        #       whis = 1.5)
        #[b.set_facecolor("#4189DD") for b in box["boxes"]]
        #[b.set_markeredgecolor("#381D59") for b in box ["fliers"]]
        #box["boxes"][0].set_facecolor("#4189DD")
        #box["fliers"][0].set_markeredgecolor("#381D59")
        #box["fliers"][0].set_markeredgecolor("red")
        
        #box["medians"][0].set_markeredgecolor("white")
        #colbox = "#364EB9"
        #plt.setp(box["medians"], color = colbox)
        #for element in ['boxes', 'whiskers', 'fliers', 'means', 'caps']:
        #    plt.setp(box[element], color = colbox)
        #for patch in box['boxes']:
        #    patch.set(facecolor = "#228FCF")
   
maxS = [2, 0.5, 1.0, 1.2    , 1.0, 4]
    
for j, a in enumerate(ax.flatten()):
     a.set_axisbelow(True)
     a.set_ylabel("Million km$^2$")
     
     a.set_xticks([2016] + endYears)
     a.set_xticklabels([str(yearbClim) + "-" + str(yeareClim) + "\nobs. climatology"] + endYears, rotation = 20)

     a.set_yticks(np.arange(0, 5, 0.5))     
     a.set_xlim(2015, 2024) 
 
     myMax = maxS[j]
     a.set_ylim(-0.0, myMax * 1.1)

     print(sectors[j])
     if sectors[j][0] == "Southern Ocean":
       a.set_yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
     a.grid()
     
plt.suptitle("Observed and SIPN South forecast February mean sea ice areas")

fig.tight_layout()
figName = "../figs/fig3_paper.png"
fig.savefig(figName, dpi = 300)
print("Figure " + figName + " printed") 

figName = "../figs/fig3_paper.JPEG"
fig.savefig(figName, dpi = 300)


