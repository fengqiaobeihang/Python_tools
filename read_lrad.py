# -*- coding: utf-8 -*-
import  math
import  numpy            as     np
from    numpy            import array
from    random           import random
from    math             import sin, sqrt
import  time        
import  datetime         as     datetime
from    datetime         import timedelta
import  netCDF4          as     nc
from    netCDF4          import Dataset
from    decimal          import *
import  matplotlib.pylab as     plt
nc_path = 'J:/lrad_ITPCAS-CMFD_V0106_B-01_03hr_010deg_197901.nc/'
out_path='J:/1/'
# outpath = 'E:/CGBM/OUTPUT/'
ncfile  = Dataset(nc_path+'lrad_ITPCAS-CMFD_V0106_B-01_03hr_010deg_197901.nc', format='NETCDF4')
lat_nc  = ncfile.dimensions['lat']
lon_nc  = ncfile.dimensions['lon']
lat_nc  = np.array(lat_nc)
lon_nc  = np.array(lon_nc)
# time    = ncfile.variables['time']
time    = 248
short_lrad=ncfile.variables['lrad']
short_lrad=np.array(short_lrad)
print 'short_lrad=',short_lrad
for i in range(time):
    lrad = short_lrad[i,:,:]
    fin_short_lrad=lrad[250,350]
    print 'fin_short_lrad',fin_short_lrad
    write_file=open(out_path +'short_lrad'+'.txt','a')
    write_file.write(str(fin_short_lrad)+'\n')
    # np.savetxt(out_path +'short_lrad'+'.txt',fin_short_lrad)