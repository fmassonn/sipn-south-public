#!/usr/bin/python
#
# Francois Massonnet francois.massonnet@uclouvain.be
# 
# Converts near-real time OSISAF sea ice concentration data
# (osi-401-b) to TECLIM compliant format
#
import numpy as np
from netCDF4 import Dataset
from datetime import date, timedelta
import os

# The (machine-dependent) location of CLIMDATA folder
data_dir = "../data/"


yearb=2021        # First year to process

# =========================


# ATTENTION, script is supposed to take end year different from
# first year to make NCO command extract work automatically

yeare=yearb + 1   # Last year to process
hemi="sh"    # Hemisphere to proces ("nh" or "sh")

# Where to read the data from
rootdir = data_dir + "/obs/ice/siconc/OSI-SAF/OSI-401-b/raw/"
# Where to write the data to
outdir =  data_dir + "/obs/ice/siconc/OSI-SAF/OSI-401-b/processed/native/"


# Check existence of paths
if not os.path.exists(rootdir):
    sys.exit("Input directory does not exist")

if not os.path.exists(outdir):
    os.makedirs(outdir)

# ========
# Input grid dimensions
if hemi == "sh":
    ny, nx = 830, 790
elif hemi == "nh":
    pass
else:
    sys.exit("Hemishpere unknown")

# Create time axis
d1 = date(yearb, 1, 1)  # Start investigated period
d2 = date(yeare, 12, 31) # End investigated period (included)
daterange = [d1 + timedelta(days=x) for x in range((d2-d1).days + 1)]
nt = len(daterange)

# Create numpy arrays with the output variables: siconc, longitude, latitude, time
siconc       = np.empty((nt, ny, nx))
siconc[:]    = np.nan
longitude    = np.empty((ny, nx))
longitude[:] = np.nan
latitude     = np.empty((ny, nx))
latitude[:]  = np.nan
time         = np.empty((nt))
time[:]      = np.nan
cellarea     = np.empty((ny, nx))
cellarea[:]  = np.nan
mask         = np.empty((ny, nx))
mask[:]      = np.nan

jt = 0 # Time counter to write the data
read_geom = True # boolean to read the geometric data only once

for day in daterange:
    print(day)
    # Check file existence
    filein = rootdir + "/ice_conc_" + hemi + "_polstere-100_multi_" + day.strftime("%Y%m%d") + "1200.nc"
    if not os.path.exists(filein):
        print("File " + filein + " not found")
    else:
        f = Dataset(filein, mode = "r")
        sic = f.variables["ice_conc"][:]
        siconc[jt, :, :] = sic[:]

        if read_geom:
            lat = f.variables["lat"][:]
            lon = f.variables["lon"][:]
            msk = f.variables["status_flag"][:]
            fillval = f.variables["ice_conc"]._FillValue

            latitude[:] = lat
            longitude[:] = lon
            mask[:] = 100.0 * (msk < 100.0)
            cellarea[:] = 10.0 * 1000.0 * 10.0 * 1000.0 # grid resolution is 10km by 10km
            read_geom = False
        f.close()
    jt += 1


# Save the data
# -------------

# 5. Save as NetCDF
# -----------------
date_ref =  date(1850, 1, 1) # zero-time reference

fileout = outdir + "siconc_SIday_OSI-401-b_r1i1p1_" + str(yearb) + "0101-" + str(yeare) + "1231_" + hemi + ".nc"
# Create file
f = Dataset(fileout, mode = "w")
# Create dimensions
time = f.createDimension('time', None)
y   = f.createDimension('y', ny)
x   = f.createDimension('x', nx)
# Create variables
times = f.createVariable('time', np.int32, ('time',))
times[:] = np.arange((d1 - date_ref).days, (d2 - date_ref).days + 1) # + 1 because of Python indexing
times.units = "days since " + str(date_ref.year) + "-" + str(date_ref.month) + "-" + str(date_ref.day)

s = f.createVariable('siconc', np.float32, ('time', 'y', 'x'), fill_value = fillval)
s.units = "%"
s.missing_value = fillval
s[:] = siconc

lons = f.createVariable('longitude', np.float32, ('y', 'x'))
lons.units = "degrees East"
lons[:] = lon

lats = f.createVariable('latitude', np.float32, ('y', 'x'))
lats.units = "degrees North"
lats[:] = lat

ca = f.createVariable('areacello', np.float32, ('y', 'x'))
ca.units = "m2"
ca[:] = cellarea

sf  = f.createVariable('sftof', np.float32, ('y', 'x'))
sf.units = "%"
sf[:] = mask

# Close
f.close()
print("")
print("File produced: " + fileout)
