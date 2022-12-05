#!/usr/bin/python
# Common sea ice diags
# Francois Massonnet & Martin Vancoppenolle
#
import numpy as np
import sys
import scipy.stats
sys.path.insert(1, "~/git/ClimateData/Visualization/")
from   ts_analyses import *

exec(open("./ts_analyses.py").read())

def compute_volume(avgthickness, cellarea, mask = 1):
  """ Input: - avgthickness: sea ice or snow volume per unit cell area, in meters
                             numpy array. if 3-D, time is assumed to be 1st
             - cellarea: array of grid cell areas (sq. meters)
             - mask (1 on ocean, 0 on continent)

      Output: Sea ice or snow volume in the region defined by the mask
  """

  import sys
  import numpy as np

  if np.max(mask) != 1.0 or np.min(mask) < 0.0:
    sys.exit("(compute_volume): mask not between 0 and 1")

  if np.max(avgthickness) > 20.0:
    print("(compute_volume): W A R N I N G: large sea ice thickness")
    print("(compute_volume): np.max(avgthickness) = " + str(np.max(avgthickness)))

  if len(avgthickness.shape) == 3:
    nt, ny, nx = avgthickness.shape
    vol = np.asarray([np.sum(avgthickness[jt, :, :] * cellarea * mask) / 1e12 for jt in range(nt)])
  elif len(avgthickness.shape) == 2:
    vol = np.sum(avgthickness * cellarea * mask) / 1e12
  else:
    sys.exit("(compute_volume): avgthickness has not 2 nor 3 dimensions")

  return vol

def compute_extent(concentration, cellarea, threshold = 15.0, mask = 1):
  """ Input: - sea ice concentration in %
               numpy array. if 3-D, time is assumed to be 1st
             - Threshold over which to consider cell as icy
             - cellarea: array of grid cell areas (sq. meters)
             - mask (1 on ocean, 0 on continent)

      Output: Sea ice extent in the region defined by the mask
  """
  import sys
  import numpy as np
  
  if np.max(concentration) < 10.0:
    sys.exit("(compute_extent): concentration seems to not be in percent")

  if len(concentration.shape) == 3:
    nt, ny, nx = concentration.shape
    ext = np.asarray([np.sum( (concentration[jt, :, :] > threshold) * cellarea * mask) / 1e12 for jt in range(nt)])
  elif len(concentration.shape) == 2:
    ext = np.sum( (concentration > threshold) * cellarea * mask) / 1e12
  else:
    sys.exit("(compute_extent): concentration has not 2 nor 3 dimensions")

  return ext
  
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

def mask_ice(concentration, threshold = 15.0):
   """ Input  - sea ice concentration in fractional units: time [optional], y, x
       Output - mask: 1 if concentration reached threshold or more at least once
   """
   import numpy as np
   import sys

   if np.max(concentration) < 10.0:
    sys.exit("(mask_ice): concentration seems to not be in percent")

   if len(concentration.shape) == 2:
     mask = 1.0 * (concentration > threshold)
   elif len(concentration.shape) == 3:
     mask = np.max(1.0 * (concentration > threshold), axis = 0)
   else:
     sys.exit("(mask_ice): input concentration has incorrect dimensions")

   return mask

def compute_ke(sispeed, concentration, cellarea, mask = 1):
  """ Input: - sispeed: numpy array of sea ice speed in m/s
                        if 3-D, one assumes that the first dim is time
             - concentration: sea ice concentration in fractional units
                        if 3-D, one assumes that the first dim is time
             - cellarea: grid cell area in m^2
             - mask (1 on ocean, 0 on continent)

      Output: Sea ice kinetic energy per unit mass, weighted by area
              Weighting by area ensures to not assign unrealistic importance
              to ice in the MIZ 
              Units are J / kg
  """

  import sys
  import numpy as np

  if len(sispeed.shape) != len(concentration.shape):
    sys.exit("(compute_ke) Sea ice speed and concentration arrays don't have matching number of dimensions")

  if len(sispeed.shape) == 3:
    nt, ny, nx = sispeed.shape
    ke = np.asarray([np.sum( (0.5 * sispeed[jt, :, :] ** 2) * concentration[jt, :, :] * cellarea * mask) / \
          np.sum(concentration[jt, :, :] * cellarea * mask) for jt in range(nt)])
  elif len(sispeed.shape) == 2:
    ny, nx = sispeed.shape
    ke = np.sum( (0.5 * sispeed ** 2) * concentration * cellarea * mask) / np.sum(concentration * cellarea * mask)
 
  return ke

def negative_seaice_feedback(volume, period, order = 1):
  """ Function to estimate the negative ice-thickness ice growth feedback
      and its significance.

      INPUTS 
        (1) volume = 1-D numpy array containing time series of sea ice volume
        (2) period = period of the signal (period expressed
                     in time steps: 12 for monthly, 365 for daily...)
        (3) order  = order of the polynomial detrending (integer >=0)

      OUTPUTS
        (1) Feedback parameter expressed as the regression
            between dV on V_min (see (3))
        (2) Correlation between those two and the p-value under the null hypothesis
            of no correlation between dV and V_min
        (3) [V_min, dV]: detrended time series of annual minimum of sea ice volume,
            detrended series of wintertime volume production
  """

  if len(volume.shape) != 1:
    sys.exit("(negative_seaice_feedback) volume is not 1-D")

  nt = len(volume)
  if  nt // period != 1.0 * nt / period:
    sys.exit("(negative_seaice_feedback) length of volume series is not multiple of period")

  # 1. Locate the minima for each year
  imin = [t + np.nanargmin(volume[t:t + period]) for t in np.arange(0, nt, period)]
  
  # 2. Locate the maxima for each year
  imax = [t + np.nanargmax(volume[t:t + period]) for t in np.arange(0, nt, period)]

  # 3. Detrend series. A one-year shift is introduced to make sure we 
  #    compute volume production *after* the summer minimum
  Vmin = detrend(volume[imin[:-1]]                   , order = order)
  dV   = detrend(volume[imax[1:]] - volume[imin[:-1]], order = order) 

  # 4. Compute diagnostics
  # If all Vmins are zero or all dVs are zero, return Nan (pathological case)
  if np.max(Vmin) == 0.0 or np.max(dV == 0.0):
    nf   = np.nan
    r    = np.nan
    pval = np.nan
  else:
    r  = np.corrcoef(Vmin, dV)[0, 1]
    N = len(Vmin)
    tstat = r / np.sqrt((1 - r ** 2) / (N - 2))  # The t-statistic.
                                                 # Under the null hypothesis of no correlation,
                                                 # tstat follows a student's law with  N - 2 dof.
    pval = 1.0 - scipy.stats.t.cdf(np.abs(tstat), N - 2)

    if pval > 0.05:
      print("(negative_seaice_feedback) W A R N I N G")
      print("                           Check the scatterplot of dV versus V_min, it is most")
      print("                           likely suspicious, and the feedback factor likely meaningless")
      print("p-value: " + str(pval))

    try:
      nf = np.polyfit(Vmin, dV, 1)[0]
    except ValueError:
      print("(negative_seaice_feedback) PROBLEM, series badly conditioned")
      print("Input volume: ")
      print(volume)
      print("Vmin: ")
      print(Vmin)
      print("dV: ")
      print(dV)
      sys.exit()
  # 5. Return
  return [nf, [r, pval], [Vmin, dV]]

def positive_seaice_feedback(volume, area, period, order = 1):
  """ Function to estimate the negative ice-thickness ice growth feedback
      and its significance.
      INPUTS
        (1) Volume = 1-D Numpy array containing time series of sea ice volume
                     Units in 10^3 km^3
        (2) Area   = 1-D Numpy array containing time series of sea ice area
                     Units in 10^6 km^2

                     Submitting consistent units between volume and area
                     is important as the feedback returned must have units 1 / m

        (3) Period = period of the signal (period expressed
                     in time steps: 12 for monthly, 365 for daily...)
        (4) Order is the order of the polynomial detrending (integer >=0)

      OUTPUTS
        (1) Feedback parameter expressed as the regression
            between dA on dV between max and min volume
        (2) Correlation between those two and the p-value under the null hypothesis
            of no correlation between dA and dV
        (3) [dV, dA]: detrended time series of summer loss of volume
                      and area
  """
  import numpy as np
  exec(open("./ts_analyses.py").read())
  import sys

  if len(volume.shape) != 1 or len(area.shape) != 1:
    sys.exit("(positive_seaice_feedback) volume or area is not 1-D")

  nt = len(volume)
  if len(area) != nt:
    sys.exit("(positive_seaice_feedback): area and volume have not the same length")

  if  nt // period != 1.0 * nt / period:
    sys.exit("(positive_seaice_feedback) length of time series is not multiple of period")

  # 1. Locate the maximum of ice volume for each year
  imax = [t + np.nanargmax(volume[t:t + period]) for t in np.arange(0, nt, period)]

  # 2. Locate the minimum of ice volume for each year
  imin = [t + np.nanargmin(volume[t:t + period]) for t in np.arange(0, nt, period)]

  # Compute area loss and volume loss for each year
  dV_raw = np.array([volume[imin[k]] - volume[imax[k]] for k in range(len(imax))])
  dA_raw = np.array([area[imin[k]]   - area[imax[k]]   for k in range(len(imax))])

  # 3. Detrend series.
  dV = detrend(dV_raw, order = order)
  dA = detrend(dA_raw, order = order)

  # 4. Compute diagnostics
  if np.max(dV) == 0.0 or np.max(dA == 0.0):
    pf = np.nan
    r = np.nan
  else:
    r  = np.corrcoef(dV, dA)[0, 1]
    N = len(dV)
    tstat = r / np.sqrt((1 - r ** 2) / (N - 2))  # The t-statistic.
                                                 # Under the null hypothesis of no correlation,
                                                 # tstat follows a student's law with  N - 2 dof.
    pval = 1.0 - scipy.stats.t.cdf(np.abs(tstat), N - 2)

    if pval > 0.05:
      print("(positive_seaice_feedback) W A R N I N G")
      print("                           Check the scatterplot of dV versus V_min, it is most")
      print("                           likely suspicious, and the feedback factor likely meaningless")
      print("p-value: " + str(pval))


    try:
      pf = np.polyfit(dV, dA, 1)[0]
    except ValueError:
      print("(positive_seaice_feedback) PROBLEM, series badly conditioned")
      print("Input volume: ")
      print(volume)
      print("Input area: ")
      print(area)
      print("dV: ")
      print(dV)
      print("dA: ")
      print(dA)
      sys.exit()

  # 5. Return
  return [pf, [r, pval], [dV, dA]]

def make_scicex_mask(lon, lat):
#  """ 
#  Input: lon and lat 
#  Output: mask with SCICEX box
#  """
  import numpy as np
  import matplotlib.path as mplPath

  if len(lon.shape) == 1 or len(lat.shape) == 1:
    lon, lat = np.meshgrid(lon, lat) # Expand if regular grid

  if np.max(lon) > 200.0:
    lon[lon > 180.0] = lon[lon > 180.0] - 360.0 # Reset to [-180, 180] if not in this range
  
  ny, nx = lon.shape
  
  mask = np.zeros((ny, nx))
  R = 6370000.0   # Earth's radius

  # The SCICEX box coordinate (Rothrock et al., JGR, 2008; Table 2)
  p = [ \
             [-15.0,  87.0 ], \
             [-60.0,  86.58], \
             [-130.0, 80.0 ], \
             [-141.0, 80.0 ], \
             [-141.0, 70.0 ], \
             [-155.0, 72.0 ], \
             [175.0,  75.50], \
             [172.0,  78.50], \
             [163.0,  80.50], \
             [126.0, 78.50 ], \
             [110.0, 84.33 ], \
             [80.0,  84.42 ], \
             [57.0,  85.17 ], \
             [33.0,  83.83 ], \
             [8.0,   84.08 ], \
      ] 
  # Convert to Lambert azimuthal equivalent (Rothrock's eq 1)
  rho = 2 * R * np.sin((45.0 - 0.5 * lat) * np.pi / 180.0)
  x   = rho *   np.cos((lon - 35.0) * np.pi / 180.0) / 1000.0
  y   = rho *   np.sin((lon - 35.0) * np.pi / 180.0) / 1000.0
 
  # Convert points to x-y
  pxy = list()
  for j in range(len(p)):
    rhotmp = 2 * R * np.sin((45.0 - 0.5 * p[j][1]) * np.pi / 180.0)
    pxy.append([rhotmp * np.cos((p[j][0] - 35.0) * np.pi / 180.0) / 1000.0, rhotmp * np.sin((p[j][0] - 35.0) * np.pi / 180.0) / 1000.0])

  bbPath = mplPath.Path(pxy)
  for jy in np.arange(ny):
    for jx in np.arange(nx):
      if bbPath.contains_point([x[jy, jx], y[jy, jx]]):
        mask[jy, jx] = 1.0

  return mask 
