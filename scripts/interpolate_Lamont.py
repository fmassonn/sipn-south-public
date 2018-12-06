#!/usr/bin/python


# Script to interpolate the monthly sea ice concentration
# data from monthly to daily (Lamont Forecast, Yuan)
#
# For each grid cell, three data points are given for 
# December, January and February average
# We are going to look for the quadratic function of time
# whose time average over December, January and February
# match the given values

# Our function is f(t) = x_1*t^2 + x_2*t + x_3 (t = 1, ... 90 as DJF has 90 days)
# Its time-average over interval t1, t2 is:
# <f> = 1 / nd * x_1 * (t2^3 - t1^3) / 3 + x_2 * (t2^2 - t1^2) / 2 + x_3 * (t2 - t1)
#     = A(t1, t2, nd) * x_1 + B(t1, t2, nd) * x_2 + C(t1, t2, nd) * x_3
#     = SIC_12
# where nd is the number of days of averaging
#
# So we need to solve a basic 3x3 system of equations
# 
# BX = D with B = [A(t11,t12) B(t11,t12) C(t11,t12) 
#                  A(t21,t22) B(t21,t22) C(t21,t22)
#                  A(t31,t32) B(t31,t32) C(t31,t32)]
# and D = [ SIC_1
#           SIC_2
#           SIC_3]

import numpy as np
from netCDF4 import Dataset

def A(t1,t2,nd):
  return 1.0 / 3.0 * (t2 ** 3 - t1 ** 3) / nd
def B(t1,t2,nd):
  return 1.0 / 2.0 * (t2 ** 2 - t1 ** 2) / nd
def C(t1,t2,nd):
  return 1.0 * (t2 - t1) /nd

t11 = 1.0
t12 = 31.0
t21 = 32.0
t22 = 62.0
t31 = 63.0
t32 = 90.0

B = np.array([[A(t11,t12,31), B(t11, t12,31), C(t11,t12,31)], \
              [A(t21,t22,31), B(t21, t22,31), C(t21,t22,31)], \
              [A(t31,t32,28), B(t31, t32,28), C(t31,t32,28)], \
             ])

Binv = np.linalg.inv(B)

# Read in the data
f = Dataset("../data/2018-2019/netcdf/Lamont_001_concentration_orig.nc", mode = "r")
sic = f.variables["siconc"][:]
lat = f.variables["latitude"][:]
lon = f.variables["longitude"][:]

nt, ny, nx = sic.shape
f.close()

time = np.arange(90)

sic_out = np.empty((90, ny, nx))
sic_out[:] = np.nan

# Mask & cell area
sftof     = np.empty((ny, nx))
areacello = np.empty((ny, nx))

R = 6378000.0 #Earth Radius

for jy in np.arange(ny):
  for jx in np.arange(nx):
    D = np.array([[sic[0, jy, jx]], [sic[1, jy, jx]], [sic[2, jy, jx]]])
    X = np.matmul(Binv, D)

    tmp = X[0] * time ** 2 + X[1] * time + X[2]

    tmp[tmp > 100.0] = 100.0
    tmp[tmp < 0.0  ] = 0.0

    sic_out[:, jy, jx] = tmp

    # No information available on mask, so I take as land those points
    # with ice for all three months
    if np.min(sic[:, jy, jx]) == 100.0:
      sftof[jy, jx] = 0.0
    else:
      sftof[jy, jx] = 100.0

    # regular grid, 0.5 degree spacing longitude, 2.0 latitude

    areacello[jy, jx] = R * np.cos(2.0 * np.pi / 360.0 * lat[jy]) *  2.0 * np.pi / 360.0 * 2.0  * \
                        R *                                         (2.0 * np.pi / 360.0 * 0.5) 

f = Dataset("../data/2018-2019/netcdf/Lamont_001_concentration.nc", mode = "w")
t_d   = f.createDimension("time", None)
lon_d = f.createDimension("longitude", len(lon))
lat_d = f.createDimension("latitude", len(lat))

mysic = f.createVariable("siconc", "f4", ("time", "latitude", "longitude",))
mylon = f.createVariable("longitude", "f4", ("longitude",))
mylat = f.createVariable("latitude", "f4", ("latitude",))
mytime= f.createVariable("time"   , "i4", ("time",))
mymask= f.createVariable("sftof", "i4", ("latitude", "longitude",))
myarea = f.createVariable("areacello", "f4", ("latitude", "longitude",))

mysic.units = "%"
mylon.units = "degrees_east"
mylat.units = "degrees_north"
mytime.units= "days since 2018-12-01"
mymask.units= "%"
myarea.units = "m2"

mylat[:] = lat
mylon[:] = lon
mytime[:] = time
mysic[:] = sic_out
mymask[:] = sftof
myarea[:]  = areacello[:]

f.close()
