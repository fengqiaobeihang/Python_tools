# -*- coding: utf-8 -*-
import time
import numpy as np
import datetime as datetime
from   datetime import timedelta
from   scipy.interpolate import griddata
import netCDF4 as nc
from   netCDF4 import Dataset

############### old region of NC data #########################
nc_path = 'H:/UHB_WRF_DATA/wrfout_heihe2014/'
outpath = 'H:/UHB_WRF_DATA/UHB2014/'

ncfile  = Dataset(nc_path+'wrfout_heihe_2014-01-01.nc', format='NETCDF4')
# presfile= Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_197901.nc', format='NETCDF4')
# sradfile= Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_197901.nc', format='NETCDF4')
# tempfile= Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_197901.nc', format='NETCDF4')
# windfile= Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_197901.nc', format='NETCDF4')
# shumfile= Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_197901.nc', format='NETCDF4')
# precfile= Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_197901.nc', format='NETCDF4')
newnc   = Dataset(nc_path+'wrfout_heihe_2014-01-01.nc', format='NETCDF4')
lat_nc  = ncfile.variables['LAT']
lon_nc  = ncfile.variables['LONG']
lat_nc  = np.array(lat_nc)
lon_nc  = np.array(lon_nc)
lon_nc,lat_nc = np.meshgrid(lon_nc,lat_nc)
print np.shape(lon_nc)
south_north = 130
west_east   = 120
newcols     = 226
newrows     = 157
lon_points  = newrows  # number of lat points, noted there is an error to be corrected
lat_points  = newcols  # number of lon points
############### time ###########################
start_time=time.time()
first=datetime.date(2014,1,1)
last =datetime.date(2015,1,1)
leng =(last-first).days
print leng
dates=[first+datetime.timedelta(n) for n in range(leng+1)]
temp_year  = 0
temp_month = 0

lat_uhh = np.loadtxt('H:/UHB_WRF_DATA/data/lat_upstreamHEIHE.txt',skiprows=6);lat_uhh=np.flipud(lat_uhh)
lon_uhh = np.loadtxt('H:/UHB_WRF_DATA/data/lon_upstreamHEIHE.txt',skiprows=6);lon_uhh=np.flipud(lon_uhh)
ivt_uhh = np.loadtxt('H:/UHB_WRF_DATA/data/ivt_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)
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
        # presfile.close()
        # sradfile.close()
        # tempfile.close()
        # windfile.close()
        # shumfile.close()
        # precfile.close()
        newnc.close()
        temp_year = iyear
        temp_month= imonth
        hour      = 0

        newnc  = Dataset(outpath+'UHB'+str(iyear)+'%02d'%(imonth)+'.nc', 'w', format='NETCDF4')    
        
        ncfile = Dataset(nc_path+'wrfout_heihe_'+str(iyear)+'-'+'%02d'%(imonth)+'-01'+'.nc', format='NETCDF4')
        old_glw_nc    = ncfile.variables['GLW']
        old_prec_nc   = ncfile.variables['GLW']
        old_psfc_nc   = ncfile.variables['PSFC']
        old_q2_nc     = ncfile.variables['Q2']
        old_swdown_nc = ncfile.variables['SWDOWN']
        old_t2_nc     = ncfile.variables['T2']
        old_u10_nc    = ncfile.variables['U10']
        old_v10_nc    = ncfile.variables['V10']
        # presfile = Dataset(nc_path+'pres/pres_ITPCAS-CMFD_V0106_B-01_'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        # old_psfc_nc = presfile.variables['pres'] 
        # sradfile = Dataset(nc_path+'srad/srad_ITPCAS-CMFD_V0106_B-01_'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        # old_swd_nc  = sradfile.variables['srad']
        # tempfile = Dataset(nc_path+'temp/temp_ITPCAS-CMFD_V0106_B-01_'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        # old_airt_nc = tempfile.variables['temp']
        # windfile = Dataset(nc_path+'wind/wind_ITPCAS-CMFD_V0106_B-01_'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        # old_wind_nc = windfile.variables['wind']
        # shumfile = Dataset(nc_path+'shum/shum_ITPCAS-CMFD_V0106_B-01_'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        # old_sh_nc   = shumfile.variables['shum']
        # precfile = Dataset(nc_path+'prec/prec_ITPCAS-CMFD_V0106_B-01_'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        # old_prec_nc = precfile.variables['prec'] 
        
        # =======creat new NC file for resampled data========
        newnc.description = 'Upstream of HEIHE Basin'
        # dimensions    
        newnc.createDimension('time', None)   
        newnc.createDimension('lat', newrows)
        newnc.createDimension('lon', newcols)
        # variables                                            
        glw_nc  = newnc.createVariable('GLW'   , 'f4', ('time', 'lat', 'lon', )) 
        psfc_nc = newnc.createVariable('PSFC'  , 'f4', ('time', 'lat', 'lon', )) 
        swd_nc  = newnc.createVariable('SWDOWN', 'f4', ('time', 'lat', 'lon', )) 
        T2_nc   = newnc.createVariable('T2'    , 'f4', ('time', 'lat', 'lon', ))
        U10m_nc = newnc.createVariable('U10'   , 'f4', ('time', 'lat', 'lon', ))
        V10m_nc = newnc.createVariable('V10'   , 'f4', ('time', 'lat', 'lon', ))
        Q2_nc   = newnc.createVariable('Q2'    , 'f4', ('time', 'lat', 'lon', ))
        prec_nc = newnc.createVariable('PREC'  , 'f4', ('time', 'lat', 'lon', ))
        
        # nc_forc_us   = newnc.createVariable('forc_us'      , 'f4', ('time', 'lat', 'lon', )) 
        # nc_forc_t    = newnc.createVariable('forc_t'       , 'f4', ('time', 'lat', 'lon', )) 
        # nc_forc_rh   = newnc.createVariable('forc_rh'      , 'f4', ('time', 'lat', 'lon', )) 
        # nc_forc_prec = newnc.createVariable('forc_prec'    , 'f4', ('time', 'lat', 'lon', )) 
        # nc_forc_psrf = newnc.createVariable('forc_psrf'    , 'f4', ('time', 'lat', 'lon', )) 
        # nc_forc_swd  = newnc.createVariable('forc_solarin' , 'f4', ('time', 'lat', 'lon', )) 
        # nc_forc_frl  = newnc.createVariable('forc_frl'     , 'f4', ('time', 'lat', 'lon', )) 

    ihour = 0
    for ihour in range(0,8):
        print ihour

        iold_glw_nc    = old_glw_nc [hour]
        iold_prec_nc   = old_prec_nc[hour]
        iold_psfc_nc   = old_psfc_nc[hour]
        iold_q2_nc     = old_q2_nc [hour]
        iold_swdown_nc = old_swdown_nc[hour]
        iold_t2_nc     = old_t2_nc[hour]
        iold_u10_nc    = old_u10_nc  [hour]
        iold_v10_nc    = old_v10_nc[hour]

        iold_prec_nc = np.where(iold_prec_nc>=0.,iold_prec_nc,0.)
        iold_prec_nc = np.where(iold_prec_nc<100.,iold_prec_nc,0.)
        iold_prec_nc = iold_prec_nc/3600.

        points = np.array([lat_nc.flatten(),lon_nc.flatten()]).T
        # Griddata methods: 'nearest','linear','cubic'
        inew_glw_nc    = griddata(points, iold_glw_nc.flatten(),    (lat_uhh,lon_uhh ), method='nearest')  
        inew_prec_nc   = griddata(points, iold_prec_nc.flatten(),   (lat_uhh,lon_uhh ), method='nearest')  
        inew_psfc_nc   = griddata(points, iold_psfc_nc.flatten(),   (lat_uhh,lon_uhh ), method='nearest')  
        inew_q2_nc     = griddata(points, iold_q2_nc.flatten(),     (lat_uhh,lon_uhh ), method='nearest')  
        inew_swdown_nc = griddata(points, iold_swdown_nc.flatten(), (lat_uhh,lon_uhh ), method='nearest')  
        inew_t2_nc     = griddata(points, iold_t2_nc.flatten(),     (lat_uhh,lon_uhh ), method='nearest')  
        inew_u10_nc    = griddata(points, iold_u10_nc.flatten(),    (lat_uhh,lon_uhh ), method='nearest')
        inew_v10_nc    = griddata(points, iold_v10_nc.flatten(),    (lat_uhh,lon_uhh ), method='nearest')   

        # tranfer to CLM format
        GLW    = inew_glw_nc*1.
        PREC   = inew_prec_nc*1.
        PSFC   = inew_psfc_nc*1.
        Q2     = inew_q2_nc*1.
        SWDOWN = inew_swdown_nc*1.
        T2     = inew_t2_nc*1.
        U10    = inew_u10_nc*1.
        V10    = inew_v10_nc*1.

        T2 = np.where(T2<=0.,270.,T2)
        T2 = np.where(T2<=180.,180.,T2)
        T2 = np.where(T2>=326.,326.,T2)

        U10    [np.isnan(U10     )]  = 2.
        V10    [np.isnan(V10     )]  = 2.
        T2     [np.isnan(T2      )]  = 273.
        Q2     [np.isnan(Q2      )]  = 0.001
        PREC   [np.isnan(PREC    )]  = 0.
        PSFC   [np.isnan(PSFC    )]  = 101325.
        SWDOWN [np.isnan(SWDOWN  )]  = 200. 
        GLW    [np.isnan(GLW     )]  = 100.

        U10     [np.abs(U10)>100.] = 2.
        U10     [np.abs(U10)<=0.]  = 0.1
        PREC    [PREC<0.]          = 0.
        PREC    [PREC>100.]        = 1.
        PSFC    [PSFC<1000.]       = 101325.
        PSFC    [PSFC>120000.]     = 101325.
        SWDOWN  [SWDOWN<0.]        = 0. 
        SWDOWN  [SWDOWN>1500.]     = 1000. 
        GLW     [GLW<10.]          = 50.
        GLW     [GLW>1000.]        = 500.
        Q2      [Q2<=0.]           = 0.001


        glw_nc   [hour,:,:] = GLW     
        psfc_nc  [hour,:,:] = PSFC
        swd_nc   [hour,:,:] = SWDOWN
        T2_nc    [hour,:,:] = T2
        U10m_nc  [hour,:,:] = U10
        V10m_nc  [hour,:,:] = V10
        Q2_nc    [hour,:,:] = Q2
        prec_nc  [hour,:,:] = PREC

        ihour = ihour+1
        hour  = hour +1
newnc.close()        
print 'end'   

 