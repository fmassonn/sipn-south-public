#!/usr/bin/python

import csv
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from netCDF4 import Dataset

R = 6378000.0

dlon, dlat = 1.0, 1.0

lat_out = np.arange(-90.0 + dlat / 2.0, 90.0   , dlat)
lon_out = np.arange(-180.0 + dlon / 2.0, 180.0 , dlon)

lon_out, lat_out = np.meshgrid(lon_out, lat_out)
ny, nx = lon_out.shape

areacello = R * np.cos( (lat_out * 2.0 * np.pi / 360.0 )) * (dlon * 2.0 * np.pi / 360.0)  * R * (dlat * 2.0 * np.pi / 360.0)

root = "/storepelican/fmasson/sipn-south-data/2018-2019/"
files = ["dic2018", "ene2019", "feb2019"]
nt = len(files)



for jt in range(nt):

    data = list()
    with open (root + "/other/" + files[jt], "r") as f:
        for row in csv.reader(f):
            data.append(np.array([float(r) for r in row[0].split()]))
    
    data = np.array(data)
    
    lat = data[:, 0]
    lon = data[:, 1] 
    are = data[:, 2]
    sic = data[:, 3]
    
    # Recenter lon
    lon[lon > 180.0] = lon[lon > 180.0] - 360.0
    
    if jt == 0:
        sic_out = np.empty((nt, ny, nx))
        sic_out[:] = np.nan
    
    for jy in np.arange(ny):
        print(str(jy) + " / " + str(ny))
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
    
# Create file
fileout = "../data/2018-2019/netcdf/barreira_001_concentration_monthly.nc"
f = Dataset(fileout, mode = "w")
# Create dimensions
time = f.createDimension('time', nt)
y =  f.createDimension('y', ny)
x =  f.createDimension('x', nx)
# Create variables
times = f.createVariable('time', np.int32, ('time',))
times.units = "months since 2018" + "-" + "12" + "-" + "01"

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
