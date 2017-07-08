# -*- coding: utf-8 -*-
'''
Created on 2016-11-1
author: shaodonghang
'''
#提取WRF数据中大冬树垭口站点的气象数据
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

############### old region of NC data #########################
nc_path = 'H:/UHB_WRF_DATA/Yakou_data/'
out_path= 'H:/UHB_WRF_DATA/Yakou_data/txtdata/200807/'

ncfile  = Dataset(nc_path+'UHB_CGBM200807.nc', format='NETCDF4')
lat_nc  = ncfile.dimensions['lat']
lon_nc  = ncfile.dimensions['lon']
lat_nc  = np.array(lat_nc)
lon_nc  = np.array(lon_nc)
#print np.shape(lon_nc)
time = 744
lat  = 157
lon  = 226
print 'The program is run!'
forc_us_nc       = ncfile.variables['forc_us']
forc_vs_nc       = ncfile.variables['forc_vs']
forc_t_nc        = ncfile.variables['forc_t']
forc_q_nc        = ncfile.variables['forc_q']
forc_rh_nc       = ncfile.variables['forc_rh']
forc_prec_nc     = ncfile.variables['forc_prec']
forc_psrf_nc     = ncfile.variables['forc_psrf']
forc_solarin_nc  = ncfile.variables['forc_solarin']
forc_frl_nc      = ncfile.variables['forc_frl']   

forc_us_nc       = np.array(forc_us_nc)
forc_vs_nc       = np.array(forc_vs_nc)
forc_t_nc        = np.array(forc_t_nc) 
forc_q_nc        = np.array(forc_q_nc)
forc_rh_nc       = np.array(forc_rh_nc)
forc_prec_nc     = np.array(forc_prec_nc)
forc_psrf_nc     = np.array(forc_psrf_nc)
forc_solarin_nc  = np.array(forc_solarin_nc)
forc_frl_nc      = np.array(forc_frl_nc)

oldwindspeed     = (forc_us_nc**2+forc_vs_nc**2)**0.5
############### cycle time ###########################
for i in range(time):
    print i
    
    windspeed         = oldwindspeed    [i,:,:]
    i_forc_t_nc       = forc_t_nc       [i,:,:]
    i_forc_q_nc       = forc_q_nc       [i,:,:]
    i_forc_rh_nc      = forc_rh_nc      [i,:,:]
    i_forc_prec_nc    = forc_prec_nc    [i,:,:]
    i_forc_psrf_nc    = forc_psrf_nc    [i,:,:]
    i_forc_solarin_nc = forc_solarin_nc [i,:,:]
    i_forc_frl_nc     = forc_frl_nc     [i,:,:]

    Yakou_windspeed       = windspeed         [124,144]
    Yakou_forc_t_nc       = i_forc_t_nc       [124,144]
    Yakou_forc_q_nc       = i_forc_q_nc       [124,144]
    Yakou_forc_rh_nc      = i_forc_rh_nc      [124,144]
    Yakou_forc_prec_nc    = i_forc_prec_nc    [124,144]
    Yakou_forc_psrf_nc    = i_forc_psrf_nc    [124,144]
    Yakou_forc_solarin_nc = i_forc_solarin_nc [124,144]
    Yakou_forc_frl_nc     = i_forc_frl_nc     [124,144]

    write_file=open(out_path +'Yakou_windspeed'+'.txt','a')
    write_file.write(str(Yakou_windspeed)+'\n')
    
    write_file=open(out_path +'Yakou_temperature'+'.txt','a')
    write_file.write(str(Yakou_forc_t_nc)+'\n')

    write_file=open(out_path +'Yakou_Precipitation'+'.txt','a')
    write_file.write(str(Yakou_forc_prec_nc)+'\n')

    write_file=open(out_path +'Yakou_forc_q_nc'+'.txt','a')
    write_file.write(str(Yakou_forc_q_nc)+'\n')

    write_file=open(out_path +'Yakou_forc_rh_nc'+'.txt','a')
    write_file.write(str(Yakou_forc_rh_nc)+'\n')

    write_file=open(out_path +'Yakou_forc_psrf_nc'+'.txt','a')
    write_file.write(str(Yakou_forc_psrf_nc)+'\n')

    write_file=open(out_path +'Yakou_forc_solarin_nc'+'.txt','a')
    write_file.write(str(Yakou_forc_solarin_nc)+'\n')

    write_file=open(out_path +'Yakou_forc_frl_nc'+'.txt','a')
    write_file.write(str(Yakou_forc_frl_nc)+'\n')

ncfile.close()    
print 'The program is end!'