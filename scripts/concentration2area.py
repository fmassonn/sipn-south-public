#!/usr/bin/python
from datetime import date
from netCDF4 import Dataset
import numpy as np

# From SIPN-South compliant sea ice concentration fields, produces 
# sea ice areas (circumpolar + regional) text files

subId = "barreiraOrig"
myyear = "2022-2023"
nmemb = 2


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

for jmemb in np.arange(1, nmemb + 1):

  jjmemb = str(jmemb).zfill(3)

  filein = "../data/" + myyear + "/netcdf/" + subId + "_" + jjmemb + "_concentration.nc"

  # Open file
  f = Dataset(filein, mode = "r")
  siconc = f.variables["siconc"][:]
  lat    = f.variables["latitude"][:]
  lon    = f.variables["longitude"][:]
  sftof  = f.variables["sftof"][:]
  areacello=f.variables["areacello"][:]
  f.close()
  ny, nx = sftof.shape


  lon[lon < 0.0] = lon[lon < 0.0] + 360.0

  area = compute_area(siconc, areacello, 1.0 * (lat < 0) * (sftof / 100.0))

  with open("../data/" + myyear + "/txt/" + subId + "_" + jjmemb + "_total-area.txt", "w") as file:
    file.write(",".join(["{0:.4f}".format(a) for a in area]))  # + 1 as python does not take the last bit
    file.write("\n")

  del area

  with open("../data/" + myyear + "/txt/" + subId + "_" + jjmemb + "_regional-area.txt", "w") as file:
    # Per longitude bin
    for j_bin in np.arange(36):
      print(j_bin)
      area = compute_area(siconc, areacello, 1.0 * (lat < 0) * (sftof / 100.0) * (lon >= j_bin * 10.0) * (lon < (j_bin + 1) * 10.0))
      file.write(",".join(["{0:.4f}".format(a) for a in area]))  # + 1 as python does not take the last bit
      file.write("\n")
