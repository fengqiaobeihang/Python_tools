# -*- coding: utf-8 -*-
import time
import numpy as np
import datetime as datetime
from   datetime import timedelta
import matplotlib.pylab as plt
from   scipy.interpolate import griddata
import netCDF4 as nc
from   netCDF4 import Dataset

############### old region of NC data #########################
nc_path = 'H:/UHB_WRF_DATA/wrfout_heihe2013/'
outpath = 'H:/UHB_WRF_DATA/UHB2013/'
# ncfile  = Dataset(nc_path+'wrfout_heihe_2014-01-01.nc', format='NETCDF4')
# lat_nc  = ncfile.variables['LAT']
# lon_nc  = ncfile.variables['LONG']
# lat_nc  = np.array(lat_nc)
# lon_nc  = np.array(lon_nc)
# print np.shape(lon_nc)
south_north = 130
west_east   = 120
newcols     = 226
newrows     = 157

############### time ###########################
start_time=time.time()
first=datetime.date(2013,1,1)
last =datetime.date(2014,1,1)
leng =(last-first).days
print leng
dates=[first+datetime.timedelta(n) for n in range(leng+1)]
temp_year  = 0
temp_month = 0

# plt.imshow(lon_uhh,origin='lower left')
# plt.colorbar()
# plt.show()
# plt.clf()


for idate in dates:
    realdate = idate; print realdate
    iyear = realdate.year
    imonth = realdate.month
    iday = realdate.day
    date_delta = realdate - datetime.date(iyear,1,1) + timedelta(days = 1)
    idays = date_delta.days-1   
    print idays
    
    if temp_year != iyear or temp_month != imonth:
   
        # ncfile.close()
        ncfile = Dataset(nc_path+'UHB'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        temp_year = iyear
        temp_month= imonth
        hour      = 0
        
        old_glw_nc  = ncfile.variables['GLW']
        old_psfc_nc = ncfile.variables['PSFC'] 
        old_swd_nc  = ncfile.variables['SWDOWN']
        old_airt_nc = ncfile.variables['T2']
        old_U10m_nc = ncfile.variables['U10']
        old_V10m_nc = ncfile.variables['V10']
        old_Q2_nc   = ncfile.variables['Q2']
        old_prec_nc = ncfile.variables['PREC'] 

    hourly_glw_nc   = np.array([])
    hourly_psfc_nc  = np.array([])
    hourly_swd_nc   = np.array([])
    hourly_airt_nc  = np.array([])
    hourly_wind_nc  = np.array([])
    hourly_rh_nc    = np.array([])
    hourly_prec_nc  = np.array([]) 
    for ihour in range(0,24):
        iold_glw_nc  = old_glw_nc [hour]
        iold_psfc_nc = old_psfc_nc[hour]
        iold_swd_nc  = old_swd_nc [hour]
        iold_airt_nc = old_airt_nc[hour]
        iold_U10m_nc = old_U10m_nc[hour]
        iold_V10m_nc = old_V10m_nc[hour]
        iold_Q2_nc   = old_Q2_nc  [hour]
        iold_prec_nc = old_prec_nc[hour]

        
        iold_glw_nc  =  np.array(iold_glw_nc )
        iold_psfc_nc =  np.array(iold_psfc_nc)
        iold_swd_nc  =  np.array(iold_swd_nc )
        iold_airt_nc =  np.array(iold_airt_nc)
        iold_U10m_nc =  np.array(iold_U10m_nc)
        iold_V10m_nc =  np.array(iold_V10m_nc)
        iold_Q2_nc   =  np.array(iold_Q2_nc  )
        iold_prec_nc =  np.array(iold_prec_nc)
        
        iold_wind    =  (iold_U10m_nc**2+iold_V10m_nc**2)**0.5
        iold_airt_nc =  iold_airt_nc-273.15
        
        a_base = 7.5+iold_airt_nc*0.
        a_base[iold_airt_nc<=0.] = 9.5
        b_base = 237.3 +iold_airt_nc*0.
        b_base[iold_airt_nc<=0.]=265.5
        
        rh = (iold_Q2_nc * iold_psfc_nc/0.622)/(6.11 * 10.**((a_base * iold_airt_nc) \
        /(iold_airt_nc + b_base))) * 100.
        
        hourly_glw_nc  = np.append(hourly_glw_nc , np.flipud(iold_glw_nc ))
        hourly_psfc_nc = np.append(hourly_psfc_nc, np.flipud(iold_psfc_nc))
        hourly_swd_nc  = np.append(hourly_swd_nc , np.flipud(iold_swd_nc ))
        hourly_airt_nc = np.append(hourly_airt_nc, np.flipud(iold_airt_nc))
        hourly_wind_nc = np.append(hourly_wind_nc, np.flipud(iold_wind))
        hourly_rh_nc   = np.append(hourly_rh_nc  , np.flipud(rh  ))
        hourly_prec_nc = np.append(hourly_prec_nc, np.flipud(iold_prec_nc))
       
        # hourly_glw_nc  = np.append(hourly_glw_nc , iold_glw_nc )
        # hourly_psfc_nc = np.append(hourly_psfc_nc, iold_psfc_nc)
        # hourly_swd_nc  = np.append(hourly_swd_nc , iold_swd_nc )
        # hourly_airt_nc = np.append(hourly_airt_nc, iold_airt_nc)
        # hourly_U10m_nc = np.append(hourly_U10m_nc, iold_U10m_nc)
        # hourly_V10m_nc = np.append(hourly_V10m_nc, iold_V10m_nc)
        # hourly_Q2_nc   = np.append(hourly_Q2_nc  , iold_Q2_nc  )
        # hourly_prec_nc = np.append(hourly_prec_nc, iold_prec_nc)
        
        hour  = hour +1  
        
    #save to TXT format also.
    hourly_glw_nc  =  np.reshape(hourly_glw_nc ,(157*24,226))
    hourly_psfc_nc =  np.reshape(hourly_psfc_nc,(157*24,226))
    hourly_swd_nc  =  np.reshape(hourly_swd_nc ,(157*24,226))
    hourly_airt_nc =  np.reshape(hourly_airt_nc,(157*24,226))
    hourly_wind_nc =  np.reshape(hourly_wind_nc,(157*24,226))
    hourly_rh_nc   =  np.reshape(hourly_rh_nc  ,(157*24,226))
    hourly_prec_nc =  np.reshape(hourly_prec_nc,(157*24,226))
    
    np.savetxt(outpath +'glw' +str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_glw_nc )
    np.savetxt(outpath +'psfc'+str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_psfc_nc)
    np.savetxt(outpath +'swd' +str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_swd_nc )
    np.savetxt(outpath +'airt'+str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_airt_nc)
    np.savetxt(outpath +'wind'+str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_wind_nc)
    # np.savetxt(outpath +'V10m'+str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_V10m_nc)
    np.savetxt(outpath +'RH'  +str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_rh_nc  )
    np.savetxt(outpath +'prec'+str(iyear) +'%02d'%(imonth)+'%03d'%(idays)+'.txt',hourly_prec_nc)
        

    #END FOR loop
#END FOR loop
newnc.close()        
print 'end'   




    
# x,y = np.meshgrid(lon_scf,lat_scf)
# points = np.array([x.flatten(),y.flatten()]).T



# oldscf = np.asfortranarray(np.array(oldscf),dtype='f')     
# objectarray = oldscf
# oldscf_day = objectarray.flatten()

# newscf = griddata(points, oldscf_day, (xi, yi),method='nearest')  #'nearest','linear','cubic'
   # newscf = np.array(newscf)
    # plt.imshow(np.flipud(newscf))
    # plt.colorbar()
    # plt.show()
    # plt.clf()
    # np.savetxt(newarea_path+str(iyear)+'%03d'%(idays+1)+'.txt',newscf)   




# ncfile.close()
            # newnc.close() 
# ncfile = Dataset('../data/babaohe_data_'+str(iyear)+'_1.nc', format='NETCDF4')
# scffile = Dataset('../data/modis_scf/'+str(iyear)+'_1.nc', format='NETCDF4')

# AIRT = ncfile.variables['temperature_K']    
# WSP = ncfile.variables['wind_speed']
# DSR = ncfile.variables['downward_shortwave_radiation']
# DLR = ncfile.variables['downward_longwave_radiation']
# PREC = ncfile.variables['precipitation_mm']
# PRES= ncfile.variables['atmosphere_pressure']
# RH = ncfile.variables['relative_humidity']