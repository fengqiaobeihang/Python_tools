# -*- coding: utf-8 -*-
'''
Created on 2015-10-15
@author: Hongyi LI
'''
import time
import numpy as np
import datetime as datetime
from   datetime import timedelta
# import matplotlib.pylab as plt
from   scipy.interpolate import griddata
import netCDF4 as nc
from   netCDF4 import Dataset

# print pyforcing.__doc__
############### old region of NC data #########################
mode     = 1 # 0:windows OS; 1:linux OS
if mode == 0:
    nc_path = 'G:/WRFdata/'
    outpath= 'E:/WRFdata/'
if mode == 1:
    nc_path = 'I:/2014WRF_WindSpeed/'
    outpath = 'I:/2014WRF_WindSpeed/'

# indexlat = np.arange(130)
# indexlon = np.arange(120)
# indexlongrid,indexlatgrid = np.meshgrid(indexlon,indexlat)
# plt.imshow(indexlatgrid,origin='ll')
# plt.colorbar()
# plt.show()
lat_uhh = np.loadtxt('para/lat_upstreamHEIHE.txt',skiprows=6);lat_uhh=np.flipud(lat_uhh)
lon_uhh = np.loadtxt('para/lon_upstreamHEIHE.txt',skiprows=6);lon_uhh=np.flipud(lon_uhh)
newnc   = Dataset(nc_path+'wrfout_heihe_2014-01-01.nc', format='NETCDF4')

wrflat = np.loadtxt('para/wrf_lat.txt')
wrflon = np.loadtxt('para/wrf_lon.txt')
wrfhgt = np.loadtxt('para/wrf_hgt.txt')

south_north = 130
west_east   = 120
newcols     = 226
newrows     = 157
lon_points  = newrows  # number of lat points, noted there is an error to be corrected
lat_points  = newcols  # number of lon points
############### time ###########################
start_time=time.time()
first=datetime.date(2014,1,1)
last =datetime.date(2014,12,31)
leng =(last-first).days
print leng
dates=[first+datetime.timedelta(n) for n in range(leng+1)]
temp_year  = 0
temp_month = 0
outdata = []
#垭口(100.2421,38.0142);祁连(100.2500,38.1833);野牛沟(99.6000,38.4333);肃南(99.6167,38.8333)
lat_yakou = 38.1833
lon_yakou = 100.2500
dis_yakou = (lat_uhh-lat_yakou)**2 +(lon_uhh-lon_yakou)**2
ivt_uhh = np.loadtxt('para/ivt_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)
mask = ivt_uhh
mask = np.where(ivt_uhh>-1,1,0)
print np.sum(mask)
for i_date in dates:
    idate    = np.array([i_date.year,i_date.month,i_date.day])
    realdate = i_date; print realdate
    iyear = realdate.year
    imonth = realdate.month
    iday = realdate.day
    date_delta = realdate - datetime.date(iyear,1,1) + timedelta(days = 1)
    idays = date_delta.days-1   
    print idays
    if temp_year != iyear or temp_month != imonth:
        newnc.close()
        newnc  = Dataset(outpath+'UHB_CGBM'+str(iyear)+'%02d'%(imonth)+'.nc',  format='NETCDF4')
        temp_year = iyear
        temp_month= imonth
        hour      = 0
    
        # old_airt_nc = newnc.variables['forc_t']
        old_us_nc   = newnc.variables['forc_us']
        old_vs_nc   = newnc.variables['forc_vs']
        old_us_nc   = np.array(old_us_nc)
        old_vs_nc   = np.array(old_vs_nc)
        old_airt_nc = (old_us_nc**2+old_vs_nc**2)**0.5
 
    ihour = 0
    for ihour in range(0,24):
        iold_airt_nc = old_airt_nc[hour]
        yakouairt = iold_airt_nc[dis_yakou==np.min(dis_yakou)] 
        outdata.append(yakouairt)
        meant = np.sum(iold_airt_nc*mask)/np.sum(mask)
        outdata.append(meant)
        ihour = ihour+1
        hour  = hour +1
    if i_date  == last:
        outdata = np.array(outdata)
        np.savetxt('2000_2010_qilian.txt',outdata)     
newnc.close() 
outdata = np.array(outdata)
np.savetxt('2000_2010_qilian.txt',outdata)       
print 'end'   
