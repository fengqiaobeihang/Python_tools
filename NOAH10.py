# -*- coding: utf-8 -*-
"""
Created on 2017/5/18
author: shaodonghang
"""
import time
import numpy             as     np
import datetime          as     datetime
import matplotlib.pylab  as     plt
import netCDF4           as     nc
import sys
import math
# import gdal
import ogr
from   gdalconst         import *
from   osr               import SpatialReference
from   osgeo             import gdal
from   osgeo.gdalconst   import *
# import libtiff           as     tif
# from   libtiff           import TIFF     #read tiff data
from   datetime          import timedelta
from   scipy.interpolate import griddata
from   netCDF4           import Dataset
# from   matplotlib.pyplot import savefig  
#---------------------README----------------
# Particle-Size Distribution: silt(SI),sand(SA),clay(CL),Soil Organic Matter(SOM), % weight, g/100g
# Bulk Density (BD),	g/cm3
# Gravel content(GRAV), % volume
############### old region of NC data #########################
#nc_path = 'F:/BaiduYunDownload/'
#outpath = 'F:/UHB_soil_result/'
###############################################################
############### old region of NC data #########################
nc_path = 'I:/zain/'
outpath = 'I:/zain/output/'
ncfile  = Dataset(nc_path+'GLDAS_NOAH10_M.A200301.001.nc', format='NETCDF4')
lat_nc  = ncfile.dimensions['g0_lat_0']
lon_nc  = ncfile.dimensions['g0_lon_1']
lat_nc  = np.array(lat_nc)
lon_nc  = np.array(lon_nc)
#print np.shape(lon_nc)
# time=672
# g0_lat_0=150
# g0_lon_1=360
SWE_nc   = ncfile.variables['SWE_GDS0_SFC'] 
SWE_nc   = np.array(SWE_nc)
#read soil and estimate data
driver  = gdal.GetDriverByName('HFA')
driver.Register()
#gdal.AllRegister() 单独注册某一类型的数据驱动，这样的话可以读也可以写，可以新建数据集
gdal.AllRegister()
fn      ='I:/zain/NOAA.tif'
ds      = gdal.Open(fn, GA_ReadOnly)
# print ds.GetDriver().ShortName
if ds is None:
   print 'Could not open ' + fn
   sys.exit(1)
#read the X and Y direction pixels of raster dataset 0, 0, cols, rows
cols    = ds.RasterXSize
print 'cols=',cols
rows    = ds.RasterYSize
print 'rows=',rows
#read snow reflectance
band1         = ds.GetRasterBand(1)
reflectance1  = band1.ReadAsArray(0, 0, cols, rows)
reflectance1  = np.array(reflectance1)

plt.imshow(reflectance1,origin='lower left')
plt.colorbar()
plt.savefig('D:\\REFLECTANCE1.jpg')  #save the figure
plt.show()
plt.clf()
#read angle data
fn_angle  ='I:/zain/NOAA.tif'
ds_angle  = gdal.Open(fn_angle, GA_ReadOnly)
# plt.imshow(angle1,origin='lower left')
# plt.colorbar()
# plt.savefig('D:\\angle1.jpg')  #save the figure
# plt.show()
# plt.clf()
driver      = ds.GetDriver()
#复制一份数据驱动
outfilename ='I:/zain/output/swe.img'
outDataset  = driver.Create(outfilename,cols,rows,1,GDT_Float32)
if outDataset is None:
    print 'Could not create albedo.tif'
    sys.exit(1)
#create new dataset
outBand=outDataset.GetRasterBand(1)
#write the data
outBand.WriteArray(SWE_nc,0,0)
#flush data to disk, set the NoData value and calculate stats
outBand.FlushCache()
outBand.SetNoDataValue(-99)
# georeference the image and set the projection
outDataset.SetGeoTransform(ds.GetGeoTransform())
outDataset.SetProjection(ds.GetProjection())
#***************************************************************** 
plt.imshow(SWE_nc,origin='lower left')
plt.colorbar()
plt.savefig('D:/SWE.jpg')  #save the figure
plt.show()
plt.clf()
print 'end' 