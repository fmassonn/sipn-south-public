#!/usr/bin/python3
#
# Script to convert NetCDF observational references of SIC to
# SIPN South compliant format (CSV)
# NSIDC 0081 and OSI-401-b
#
# Author - F. Massonnet
# Date   - March 5, 2018


def convert(target):
  """
  expName = "2021-2022", or "TOP2022", ...
  """

  # Imports, modules, etc.
  import numpy as np
  from netCDF4 import Dataset
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  from datetime import date, timedelta
  import os
  
  
  # Name of obs to process
  #        OBS LABEL      OBS DIRECTORY
  obs = [ 
          ["OSI-401-b" , "../data/" + target +  "/netcdf/"], \
          ["NSIDC-0081", "../data/" + target +  "/netcdf/"], \
        ]
  
  # ========================================
  
  
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
  
  
  # For internal check
  plt.figure(figsize = (4, 4))
  
  for j_obs in range(len(obs)):
      print(obs[j_obs][0])
      # Input file, following CMIP conventions
      filein = obs[j_obs][1] + obs[j_obs][0] + "_000_concentration.nc"
      print(filein)
  
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
      
      # Compute sea ice area for that period
      # ------------------------------------
      areatot = compute_area(siconc[:, :, :], cellarea, mask = 1.0 * (lat < 0.0)) # + 1 because of Python indexing convention
      print(areatot)
      # Save as CSV file
      # ----------------
      # Total area
      with open("../data/" + target + "/txt/" + obs[j_obs][0] + "_000" + "_total-area.txt", "w") as file:
          file.write(",".join(["{0:.4f}".format(a) for a in areatot]))  
          file.write("\n")
      
      # Per longitude
      with open("../data/" + target + "/txt/" + obs[j_obs][0] + "_000" + "_regional-area.txt", "w") as file:
          # Per longitude bin
          for j_bin in np.arange(36):
            print(j_bin)
            area = compute_area(siconc[:, :, :], cellarea, mask = 1.0 * (lat < 0) * (sftof == 100.0) * (lon >= j_bin * 10.0) * (lon < (j_bin + 1) * 10.0))
            file.write(",".join(["{0:.4f}".format(a) for a in area]))  # + 1 as python does not take the last bit
            file.write("\n")
      
  
  plt.legend()
  plt.ylim(0.0, 17.0)
  plt.savefig("../figs/obs.png", dpi = 300)


if __name__ == "__main__":
  import sys
  if len(sys.argv) != (1 + 1):
    print("No argument given")
    sys.exit
  else:
    target = sys.argv[1]

    convert(target)
