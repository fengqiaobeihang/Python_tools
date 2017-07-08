resample_forcing_ICEv2

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
import forcing_class
forc = forcing_class.forcing_calc()
# print pyforcing.__doc__
############### old region of NC data #########################
mode     = 1 # 0:windows OS; 1:linux OS
if mode == 0:
    nc_path = 'G:/WRFdata/'
    outpath= 'E:/WRFdata/'
if mode == 1:
    nc_path = '../../data/WRFdatav1/wrfout_heihe1/'
    outpath= '../../data/hybird/'

ncfile  = Dataset(nc_path+'wrfout_heihe_2000-02-01.nc', format='NETCDF4')
newnc   = Dataset(nc_path+'wrfout_heihe_2000-03-01.nc', format='NETCDF4')
lat_nc  = ncfile.variables['LAT']
lon_nc  = ncfile.variables['LONG']
lat_nc  = np.array(lat_nc)
lon_nc  = np.array(lon_nc)
print np.shape(lon_nc)
south_north = 130
west_east   = 120
newcols     = 226
newrows     = 157
lon_points  = newrows  # number of lat points, noted there is an error to be corrected
lat_points  = newcols  # number of lon points
############### time ###########################
start_time=time.time()
first=datetime.date(2000,1,1)
last =datetime.date(2016,1,1)
leng =(last-first).days
print leng
dates=[first+datetime.timedelta(n) for n in range(leng+1)]
temp_year  = 0
temp_month = 0

lat_uhh = np.loadtxt('data/lat_upstreamHEIHE.txt',skiprows=6);lat_uhh=np.flipud(lat_uhh)
lon_uhh = np.loadtxt('data/lon_upstreamHEIHE.txt',skiprows=6);lon_uhh=np.flipud(lon_uhh)
ivt_uhh = np.loadtxt('data/ivt_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)
# plt.imshow(lon_uhh,origin='lower left')
# plt.colorbar()
# plt.show()
# plt.clf()

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
        ncfile.close()
        newnc.close()
        ncfile = Dataset(nc_path+'wrfout_heihe_'+str(iyear)+'-'+'%02d'%(imonth)+'-01.nc', format='NETCDF4')
        newnc  = Dataset(outpath+'UHB_CGBM'+str(iyear)+'%02d'%(imonth)+'_v4.nc', 'w', format='NETCDF4')
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
        old_prec_nc = ncfile.variables['RAINC'] 
        
        # =======creat new NC file for resampled data========
        newnc.description = 'Upstream of HEIHE Basin'
        # dimensions    
        newnc.createDimension('time', None)   
        newnc.createDimension('lat', newrows)
        newnc.createDimension('lon', newcols)
        # variables                                            
        # glw_nc  = newnc.createVariable('GLW'   , 'f4', ('time', 'lat', 'lon', )) 
        # psfc_nc = newnc.createVariable('PSFC'  , 'f4', ('time', 'lat', 'lon', )) 
        # swd_nc  = newnc.createVariable('SWDOWN', 'f4', ('time', 'lat', 'lon', )) 
        # airt_nc = newnc.createVariable('T2'    , 'f4', ('time', 'lat', 'lon', ))
        # U10m_nc = newnc.createVariable('U10'   , 'f4', ('time', 'lat', 'lon', ))
        # V10m_nc = newnc.createVariable('V10'   , 'f4', ('time', 'lat', 'lon', ))
        # Q2_nc   = newnc.createVariable('Q2'    , 'f4', ('time', 'lat', 'lon', ))
        # prec_nc = newnc.createVariable('PREC'  , 'f4', ('time', 'lat', 'lon', ))
        
        nc_forc_us   = newnc.createVariable('forc_us'      , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_vs   = newnc.createVariable('forc_vs'      , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_t    = newnc.createVariable('forc_t'       , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_q    = newnc.createVariable('forc_q'       , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_rh   = newnc.createVariable('forc_rh'      , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_prec = newnc.createVariable('forc_prec'    , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_psrf = newnc.createVariable('forc_psrf'    , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_swd  = newnc.createVariable('forc_solarin' , 'f4', ('time', 'lat', 'lon', )) 
        nc_forc_frl  = newnc.createVariable('forc_frl'     , 'f4', ('time', 'lat', 'lon', )) 

    ihour = 0
    for ihour in range(0,24):
        iold_glw_nc  = old_glw_nc [hour]
        iold_psfc_nc = old_psfc_nc[hour]
        iold_swd_nc  = old_swd_nc [hour]
        iold_airt_nc = old_airt_nc[hour]
        iold_U10m_nc = old_U10m_nc[hour]
        iold_V10m_nc = old_V10m_nc[hour]
        iold_Q2_nc   = old_Q2_nc  [hour]
        iold_prec_nc = old_prec_nc[hour]

        iold_prec_nc = np.where(iold_prec_nc>=0.,iold_prec_nc,0.)
        iold_prec_nc = np.where(iold_prec_nc<100.,iold_prec_nc,0.)
        iold_prec_nc = iold_prec_nc/3600.

        points = np.array([lat_nc.flatten(),lon_nc.flatten()]).T
        # Griddata methods: 'nearest','linear','cubic'
        inew_glw_nc  = griddata(points, iold_glw_nc .flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_psfc_nc = griddata(points, iold_psfc_nc.flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_swd_nc  = griddata(points, iold_swd_nc .flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_airt_nc = griddata(points, iold_airt_nc.flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_U10m_nc = griddata(points, iold_U10m_nc.flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_V10m_nc = griddata(points, iold_V10m_nc.flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_Q2_nc   = griddata(points, iold_Q2_nc  .flatten(), (lat_uhh,lon_uhh ),method='linear')  
        inew_prec_nc = griddata(points, iold_prec_nc.flatten(), (lat_uhh,lon_uhh ),method='linear')  

        # tranfer to CLM format
        forc_us      = inew_U10m_nc*1.
        forc_vs      = inew_V10m_nc*1.
        forc_frl     = inew_glw_nc*1.
        forc_t       = inew_airt_nc*1.
        forc_q       = inew_Q2_nc*1.
        forc_psrf    = inew_psfc_nc*1.
        forc_prec    = inew_prec_nc*1.
        forc_solarin = inew_swd_nc*1.

        forc_t = np.where(forc_t<=0.,270.,forc_t)
        forc_t = np.where(forc_t<=180.,180.,forc_t)
        forc_t = np.where(forc_t>=326.,326.,forc_t)

        forc_us     [np.isnan(forc_us     )] = 2.
        forc_vs     [np.isnan(forc_vs     )] = 2.
        forc_t      [np.isnan(forc_t      )] = 273.
        forc_q      [np.isnan(forc_q      )] = 0.001
        forc_prec   [np.isnan(forc_prec   )] = 0.
        forc_psrf   [np.isnan(forc_psrf   )] = 101325.
        forc_solarin[np.isnan(forc_solarin)] = 200. 
        forc_frl    [np.isnan(forc_frl    )] = 100.

        forc_us     [np.abs(forc_us)>100.] = 2.
        forc_vs     [np.abs(forc_vs)>100.] = 2.
        forc_us     [np.abs(forc_us)<=0.]  = 0.1
        forc_vs     [np.abs(forc_vs)<=0.]  = 0.1        
        forc_prec   [forc_prec<0.]         = 0.
        forc_prec   [forc_prec>100.]       = 1.
        forc_psrf   [forc_psrf<1000.]      = 101325.
        forc_psrf   [forc_psrf>120000.]    = 101325.
        forc_solarin[forc_solarin<0.]      = 0. 
        forc_solarin[forc_solarin>1500.]   = 1000. 
        forc_frl    [forc_frl<10.]         = 50.
        forc_frl    [forc_frl>1000.]       = 500.
        forc_q      [forc_q<=0.]           = 0.001

        forc_rh, forc_q = forc.etorh(forc_q,forc_t,forc_psrf)

        nc_forc_us    [hour,:,:] = forc_us     
        nc_forc_vs    [hour,:,:] = forc_vs     
        nc_forc_t     [hour,:,:] = forc_t      
        nc_forc_q     [hour,:,:] = forc_q      
        nc_forc_rh    [hour,:,:] = forc_rh     
        nc_forc_prec  [hour,:,:] = forc_prec    
        nc_forc_psrf  [hour,:,:] = forc_psrf   
        nc_forc_swd   [hour,:,:] = forc_solarin  
        nc_forc_frl   [hour,:,:] = forc_frl    

        ihour = ihour+1
        hour  = hour +1
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

