#!/usr/bin/python

import csv
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

myyear = "2019-2020"

R = 6378000.0

dlon, dlat = 1.0, 1.0

lat_out = np.arange(-90.0 + dlat / 2.0, 90.0   , dlat)
lon_out = np.arange(-180.0 + dlon / 2.0, 180.0 , dlon)

lon_out, lat_out = np.meshgrid(lon_out, lat_out)
ny, nx = lon_out.shape

areacello = R * np.cos( (lat_out * 2.0 * np.pi / 360.0 )) * (dlon * 2.0 * np.pi / 360.0)  * R * (dlat * 2.0 * np.pi / 360.0)

root = "../data/" + myyear + "/barreira/"
files = ["SHN_ " + str(j) + " concentracion.txt" for j in range(1, 90 + 1)]

nt = len(files)


for jt in range(nt):
    print(jt)

    data = list()
    
    if os.path.exists(root + files[jt]):

        with open (root + files[jt], "r") as f:
            csvreader = csv.reader(f)
            next(csvreader)
            for row in csvreader:
                data.append(np.array([float(r) for r in row[0].split()]))
        
        data = np.array(data)
        
        lat = data[:, 0]
        lon = data[:, 1] 
        are = data[:, 2]
        sic = data[:, 3]
        
        # Recenter longitude to [-180, 180]
        lon[lon > 180.0] = lon[lon > 180.0] - 360.0
        
        # Initialize array if first time step
        if jt == 0:
            sic_out = np.full((nt, ny, nx), np.nan)
        
        for jy in np.arange(ny):
            for jx in np.arange(nx):
                # Locate all points that fall in that box
                lonmin = lon_out[jy, jx] - dlon / 2.0
                lonmax = lon_out[jy, jx] + dlon / 2.0
                latmin = lat_out[jy, jx] - dlat / 2.0
                latmax = lat_out[jy, jx] + dlat / 2.0
        
                mask = (lat >= latmin) * (lat < latmax) * (lon >= lonmin) * (lon < lonmax)
        
                if np.sum(mask) == 0:
                    sic_out[jt, jy, jx] = 0.0
                else:
                    sic_out[jt, jy, jx] = np.mean(sic[mask])
    else:
        print("File " + root + files[jt] + " not found, skipping")
    
# Create file
fileout = "../data/" + myyear + "/netcdf/barreira_001_concentration.nc"
f = Dataset(fileout, mode = "w")

# Create dimensions
time = f.createDimension('time', nt)
y =  f.createDimension('y', ny)
x =  f.createDimension('x', nx)
# Create variables
times = f.createVariable('time', np.int32, ('time',))
times.units = "days since 2018 " + myyear[:4]  + "-" + "12" + "-" + "01"

h = f.createVariable("siconc", np.float32, ('time', 'y', 'x'))
h.units = "%"
h[:] = sic_out

lons = f.createVariable("longitude", np.float32, ("y", "x"))
lons[:] = lon_out

lats = f.createVariable("latitude", np.float32, ("y", "x"))
lats[:] = lat_out

cellarea = f.createVariable("areacello", np.float32, ("y", "x"))
cellarea.units = "m2"
cellarea[:] = areacello

# Close
f.close()
