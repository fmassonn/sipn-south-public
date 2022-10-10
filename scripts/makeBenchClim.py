#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:15:47 2022

@author: massonnetf
"""

# Creates benchmark climatology forecasts by retrieving the 
# observations of previous years


import numpy as np
from datetime import datetime
from datetime import timedelta
import os
from netCDF4 import Dataset


# Computation of areas given sea ice concentration
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


# Source dir
climDir = "/Users/massonnetf/CLIMDATA/"

# Year of the 1st of December starting the summer season
# Last year to be processed. That's the year of the forecast minus one
#
# !!! ATTENTION !!!!
# The year is referenced to the 1st of December of the summer season, so
# 2020 --> 2020-2021 = the last season before 2021-2022
target = "2021-2022"
yeare = int(target[:4]) - 1

# Year to start the collection of climatology. ATTENTION to the convention
yearb = yeare - 30 + 1 



# Member ID to be incremented each year
jMemb = 1

for year in np.arange(yearb, yeare + 1):
    print(year)
    
    #  Create date axis
    d0 = datetime(year, 12, 1    )
    d1 = datetime(year + 1, 2, 28)
    nd = (d1 - d0).days + 1
    
    dates = [d0 + timedelta(days = int(d)) for d in np.arange(0, nd)]
    
    
    listSia = list()
    # Loop through files
    for d in dates:
        thisYear = d.year
        thisMonth= d.month
        thisDay  = d.day
        
        dirData = climDir + "/obs/ice/siconc/OSI-SAF/OSI-450/raw/"
        
        if thisYear < 2015:
          dirData = climDir + "/obs/ice/siconc/OSI-SAF/OSI-450/raw/"
          fileName = "ice_conc_sh_ease2-250_cdr-v2p0_" + str(thisYear) + str(thisMonth).zfill(2) + str(thisDay).zfill(2) + "1200.nc"
          cellArea = (25e3) ** 2
        else:
          dirData = climDir + "/obs/ice/siconc/OSI-SAF/OSI-401-b/raw/"
          fileName = "ice_conc_sh_polstere-100_multi_" + str(thisYear) + str(thisMonth).zfill(2) + str(thisDay).zfill(2) + "1200.nc"
          cellArea = (10e3) ** 2
        if os.path.exists(dirData + fileName):
            
            f = Dataset(dirData + fileName, mode = "r")
            sic = f.variables["ice_conc"][:]

            listSia.append(compute_area(sic, cellArea)[0])
            f.close()
        else:
           print("File not found: " + dirData + fileName)
           listSia.append(np.nan)
            
    # Write to txt file
    with open("../data/" + target + "/txt/benchClim_" + str(jMemb).zfill(3) + "_total-area.txt", "w") as file:
        file.write(",".join(["{0:.4f}".format(a) for a in listSia]))  
        file.write("\n")
        
        
    jMemb += 1
