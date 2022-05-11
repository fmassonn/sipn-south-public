#!/usr/bin/python


# Script to interpolate the monthly sea ice concentration
# data from monthly to daily (Lamont Forecast, Yuan)
#
# For each grid cell, three or four data points are given by Yuan for 
# (December, January and February averages or May, June, July, August)
#
# We are going to look for the polynomial (quadratic or cubic) function of time
# whose time average over the given months
# match the given values

# Our function is f(t) = [(x_0) * t^3] + x_1*t^2 + x_2*t + x_3 (t = 1, ... 90 or 123 since DJF or MJJA
# has 90 days or 123 days)
#
# The example below is for the 3-month but can be extende to 4 months:

# Its time-average over interval t1, t2 is:
# <f> = 1 / nd * [x_1 * (t2^3 - t1^3) / 3 + 
#                 x_2 * (t2^2 - t1^2) / 2 + 
#                 x_3 * (t2 - t1)]
#     = A(t1, t2, nd) * x_1 + B(t1, t2, nd) * x_2 + C(t1, t2, nd) * x_3
#     = SIC_12
# where nd is the number of days of averaging
#
# So we need to solve a basic 3x3 system of equations
# 
# DX = E with D = [A(t11,t12) B(t11,t12) C(t11,t12) 
#                  A(t21,t22) B(t21,t22) C(t21,t22)
#                  A(t31,t32) B(t31,t32) C(t31,t32)]
# and D = [ SIC_12
#           SIC_34
#           SIC_56
#         ]

import numpy as np
from netCDF4 import Dataset

myyear = "TOP2022"

if myyear == "TOP2022":
  nmonth = 4
  ndout = 123
else:
  nmonth = 3
  ndout = 90

if nmonth == 3:

    def A(t1, t2, nd):
      return 1.0 / 3.0 * (t2 ** 3 - t1 ** 3) / nd
    def B(t1, t2, nd):
      return 1.0 / 2.0 * (t2 ** 2 - t1 ** 2) / nd
    def C(t1, t2, nd):
      return 1.0 * (t2 - t1) /nd

    t11 = 1.0
    t12 = 31.0
    t21 = 32.0
    t22 = 62.0
    t31 = 63.0
    t32 = 90.0

    Btmp = np.array([[A(t11,t12,31), B(t11, t12,31), C(t11,t12,31)], \
                     [A(t21,t22,31), B(t21, t22,31), C(t21,t22,31)], \
                     [A(t31,t32,28), B(t31, t32,28), C(t31,t32,28)], \
                 ])

    Binv = np.linalg.inv(Btmp)
elif nmonth == 4:
    def A(t1, t2, nd):
      return 1.0 / 4.0 * (t2 ** 4 - t1 ** 4) / nd
    def B(t1, t2, nd):
      return 1.0 / 3.0 * (t2 ** 3 - t1 ** 3) / nd
    def C(t1, t2, nd):
      return 1.0 / 2.0 * (t2 ** 2 - t1 ** 2) / nd
    def D(t1, t2, nd):
      return 1.0 * (t2 - t1) /nd

    t11 = 1.0
    t12 = 31.0
    t21 = 32.0
    t22 = 61.0
    t31 = 62.0
    t32 = 92.0
    t41 = 93.0
    t42 = 123.0

    Btmp = np.array([[A(t11,t12,31), B(t11, t12,31), C(t11,t12,31), D(t11, t12, 31)], \
                     [A(t21,t22,30), B(t21, t22,30), C(t21,t22,30), D(t21, t22, 30)], \
                     [A(t31,t32,31), B(t31, t32,31), C(t31,t32,31), D(t31, t32, 31)], \
                     [A(t41,t42,31), B(t41, t42,31), C(t41,t42,31), D(t41, t42, 31)], \
                 ])

    Binv = np.linalg.inv(Btmp)
# Read in the data
f = Dataset("../data/" + myyear + "/netcdf/Lamont_001_concentration_orig.nc", mode = "r")
sic = f.variables["siconc"][:]
lat = f.variables["latitude"][:]
lon = f.variables["longitude"][:]

nt, ny, nx = sic.shape
f.close()

time = np.arange(ndout)

sic_out = np.empty((ndout, ny, nx))
sic_out[:] = np.nan

# Mask & cell area
sftof     = np.empty((ny, nx))
areacello = np.empty((ny, nx))

R = 6378000.0 #Earth Radius

for jy in np.arange(ny):
  for jx in np.arange(nx):

    if myyear == "TOP2022":
      D = np.array([[sic[0, jy, jx]], [sic[1, jy, jx]], [sic[2, jy, jx]], [sic[3, jy, jx]]])
    else:
      D = np.array([[sic[0, jy, jx]], [sic[1, jy, jx]], [sic[2, jy, jx]]])

    X = np.matmul(Binv, D)

    if myyear == "TOP2022":
      tmp = X[0] * time ** 3 + X[1] * time ** 2 + X[2] * time + X[3]
    else:
      tmp = X[0] * time ** 2 + X[1] * time      + X[0]

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

f = Dataset("../data/" + myyear + "/netcdf/Lamont_001_concentration.nc", mode = "w")
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
if myyear == "TOP2022":
  mytime.units= "days since " + myyear[-4:] + "-05-01"
else:
  mytime.units= "days since " + myyear[:4] + "-12-01"
mymask.units= "%"
myarea.units = "m2"

mylat[:] = lat
mylon[:] = lon
mytime[:] = time
mysic[:] = sic_out
mymask[:] = sftof
myarea[:]  = areacello[:]

f.close()



# Interpolation of areas calculated by X. Yuan herself (text files).
# ------------------------------------------------------------------
import pandas as pd
# Total area
# ----------
# Read the CSV file
filein = "../data/" + myyear + "/txt/Lamont_001_total-area_orig.txt"
csv = pd.read_csv(filein, header = None)
series = csv.iloc[0][:]

if myyear == "TOP2022":
  D = np.array([[series[0]], [series[1]], [series[2]], [series[3]] ])
else:
  D = np.array([[series[0]], [series[1]], [series[2]]              ])

X = np.matmul(Binv, D)

if myyear == "TOP2022":
  newseries = X[0] * time ** 3 + X[1] * time ** 2 + X[2] * time + X[3]
else:
  newseries = X[0] * time ** 2 + X[1] * time + X[2]


# Write
with open("../data/" + myyear + "/txt/Lamont_001_total-area.txt", "w") as file:
  file.write(",".join(["{0:.4f}".format(a) for a in newseries]))  # + 1 as python does not take the last bit
  file.write("\n")


# Regional areas
# --------------
filein = "../data/" + myyear + "/txt/Lamont_001_regional-area_orig.txt"

csv = pd.read_csv(filein, header = None)
with open("../data/" + myyear + "/txt/Lamont_001_regional-area.txt", "w") as file2:
  for j in np.arange(36):
    series = csv.iloc[j][:]
    if myyear == "TOP2022":
      D = np.array([[series[0]], [series[1]], [series[2]], [series[3]] ])
    else:
      D = np.array([[series[0]], [series[1]], [series[2]] ])

    X = np.matmul(Binv, D)

    if myyear == "TOP2022":
      newseries = X[0] * time ** 3 + X[1] * time ** 2 + X[2] * time + X[3]
    else:
      newseries = X[0] * time ** 2 + X[1] * time + X[2]
  
    file2.write(",".join(["{0:.4f}".format(a) for a in newseries]))  # + 1 as python does not take the last bit
    file2.write("\n")

