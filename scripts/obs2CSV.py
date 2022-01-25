#!/usr/bin/python
#
# Script to convert NetCDF observational references of SIC to
# SIPN South compliant format (CSV), i.e. text files with
# areas for the period 1 Dec year -- 28 Feb year + 1

# NSIDC 0081 and OSI-401b
#
# Author - F. Massonnet
# Date   - March 5, 2018
#        - January 2022: update to work on 2015-now data

# Imports, modules, etc.
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import date, timedelta
import os


# The (machine-dependent) location of CLIMDATA folder
data_dir = os.environ["TECLIM_CLIMATE_DATA"]

# Name of obs to process
#        OBS LABEL      OBS DIRECTORY
obs = [ 
        ["OSI-401-b" , data_dir + "/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/"], \
        ["NSIDC-0081", data_dir + "/obs/ice/siconc/NSIDC/NSIDC-0081/processed/native/"], \
      ]

# Date parameters
# Years for which a CSV file need to be created
# The years always refer to the year of the month of December of the summer season

yearb = 2015
yeare = 2021




# =================

# Function to compute sea ice area from sea ice concentation
# -----
def compute_area(concentration, cellarea, mask = 1):
  """ Input: - sea ice concentration in %
               numpy array. if 3-D, time is assumed to be 1st
             - cellarea: array of grid cell areas (sq. meters)
             - mask (1 on ocean, 0 on continent)

      Output: Sea ice area in the region defined by the mask
  """
  import sys
  import numpy as np

  if np.max(concentration) < 10.0:
    sys.exit("(compute_area): concentration seems to not be in percent")

  if len(concentration.shape) == 3:
    nt, ny, nx = concentration.shape
    are = np.asarray([np.sum( concentration[jt, :, :] / 100.0 * cellarea * mask) / 1e12 for jt in range(nt)])
  elif len(concentration.shape) == 2:
    are = np.sum( concentration / 100.0 * cellarea * mask) / 1e12
  else:
    sys.exit("(compute_area): concentration has not 2 nor 3 dimensions")

  return are
# ------
# ------


for year in np.arange(yearb, yeare + 1):
    yearp1 = year + 1
    print(str(year) + "-" + str(yearp1))
    target = str(year) + "-" + str(yearp1)    # Where to place the output file (../data/$target/txt)
    
    d0 = date(1850,   1, 1 )   # Zero-time reference of the input file
    d1 = date(year,  12, 1 )   # Start investigated period
    d2 = date(yearp1, 2, 28)   # End investigated period (include)
    

    daterange = [d1 + timedelta(days=x) for x in range((d2-d1).days + 1)]

    # For internal check
    plt.figure(figsize = (4, 4))

    for j_obs in range(len(obs)):
        
        print(obs[j_obs][0])
        # Input file, following CMIP conventions
        filein = obs[j_obs][1] + "siconc_SIday_" + obs[j_obs][0] + "_r1i1p1_" + str(yearb) + "0101-" + str(yeare + 1) + "1231_sh.nc"
        print("Source file: " + filein)

        f = Dataset(filein, mode = "r")
        siconc = f.variables["siconc"][:]
        time   = f.variables["time"][:]
        cellarea = f.variables["areacello"][:]
        sftof    = f.variables["sftof"][:]
        lat      = f.variables["latitude"][:]
        lon      = f.variables["longitude"][:]
        # Re-range longitude to [0, 360.0]
        lon[lon < 0.0] = lon[lon < 0.0] + 360.0
        f.close()
        
        # Subset to the forecasting season
        # --------------------------------
        t1 = (d1 - d0).days - time[0]
        t2 = (d2 - d0).days - time[0]
        
        # Compute sea ice area for that period
        # ------------------------------------
        areatot = compute_area(siconc[t1:t2 + 1, :, :], cellarea, mask = 1.0 * (lat < 0.0)) # + 1 because of Python indexing convention
        print(areatot)
        # Save as CSV file
        # ----------------

	#
        # Total area
        print("  Doing global")
        with open("../data/" + target + "/txt/" + obs[j_obs][0] + "_000" + "_total-area.txt", "w") as file:
            file.write(",".join(["{0:.4f}".format(a) for a in areatot]))  
            file.write("\n")
        
        # Per longitude
        print("  Doing regional")
        with open("../data/" + target + "/txt/" + obs[j_obs][0] + "_000" + "_regional-area.txt", "w") as file:
            # Per longitude bin
            for j_bin in np.arange(36):
              # print(j_bin)
              area = compute_area(siconc[t1:t2 + 1, :, :], cellarea, mask = 1.0 * (lat < 0) * (sftof == 100.0) * (lon >= j_bin * 10.0) * (lon < (j_bin + 1) * 10.0))
              file.write(",".join(["{0:.4f}".format(a) for a in area]))  # + 1 as python does not take the last bit
              file.write("\n")
        
        # Plot for internal check
        # -----------------------
        plt.plot(daterange, areatot, label = obs[j_obs][0])

plt.legend()
plt.ylim(0.0, 17.0)
plt.savefig("../figs/obs.png", dpi = 300)
    
