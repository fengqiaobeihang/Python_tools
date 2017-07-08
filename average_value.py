'''
Created on 2016-8-1
author: shaodonghang
'''
import time
import numpy            as  np
import datetime         as  datetime
from   datetime             import timedelta
import matplotlib.pylab as  plt
from   scipy.interpolate    import griddata
import netCDF4          as  nc
from   netCDF4              import Dataset

############### old region of NC data #########################
nc_path = 'E:/CGBM/'
outpath = 'E:/CGBM/OUTPUT/'
ncfile  = Dataset(nc_path+'UHB_CGBM200801.nc', format='NETCDF4')
lat_nc  = ncfile.dimensions['lat']
lon_nc  = ncfile.dimensions['lon']
lat_nc  = np.array(lat_nc)
lon_nc  = np.array(lon_nc)
#print np.shape(lon_nc)
time=672
lat=157
lon=226

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
windspeed        = (forc_us_nc**2+forc_vs_nc**2)**0.5
temp             = forc_t_nc-273.15
print 'temp=',temp
sum_ws      = [0]
sum_t       = [0]
sum_q       = [0]
sum_rh      = [0]
sum_prec    = [0]
sum_psrf    = [0]
sum_solarin = [0]
sum_frl     = [0]
print 'forc_us_nc=',forc_us_nc
############### cycle time ###########################
for i in range(time):
    print i
    sum_ws      = sum_ws      + windspeed[i]
    sum_t       = sum_t       + temp[i]
    sum_q       = sum_q       + forc_q_nc[i]
    sum_rh      = sum_rh      + forc_rh_nc[i]
    sum_prec    = sum_prec    + forc_prec_nc[i]
    sum_psrf    = sum_psrf    + forc_psrf_nc[i]
    sum_solarin = sum_solarin + forc_solarin_nc[i]
    sum_frl     = sum_frl     + forc_frl_nc[i]
    
    print 'sum_ws=',sum_ws
    print 'sum_t=',sum_t
    print 'sum_q=',sum_q
    print 'sum_rh=',sum_rh
    print 'sum_prec=',sum_prec
    print 'sum_psrf=',sum_psrf
    print 'sum_solarin=',sum_solarin
    print 'sum_frl=',sum_frl
    
    forc_ws_nc_mean      = sum_ws/672
    forc_t_nc_mean       = sum_t/672
    forc_q_nc_mean       = sum_q/672
    forc_rh_nc_mean      = sum_rh/672
    forc_prec_nc_mean    = sum_prec/672
    forc_psrf_nc_mean    = sum_psrf/672
    forc_solarin_nc_mean = sum_solarin/672
    forc_frl_nc_mean     = sum_frl/672
    # average_us       = np.mean(forc_us_nc)
    # average_vs       = np.mean(forc_vs_nc)
    # average_t        = np.mean(forc_t_nc)
    # average_q        = np.mean(forc_q_nc)
    # average_rh       = np.mean(forc_rh_nc)
    # average_prec     = np.mean(forc_prec_nc)
    # average_psrf     = np.mean(forc_psrf_nc)
    # average_solarin  = np.mean(forc_solarin_nc)
    # average_frl      = np.mean(forc_frl_nc)
# write_file=open(outpath +'Yakou_windspeed'+'.txt','a')
# # write_file.write(str(forc_us_nc_mean)+'\n')
# write_file.write(str(forc_us_nc_mean))
np.savetxt(outpath +'windspeed'+'.txt',forc_ws_nc_mean)
np.savetxt(outpath +'T'+'.txt',forc_t_nc_mean)
np.savetxt(outpath +'rh'+'.txt',forc_rh_nc_mean)
np.savetxt(outpath +'prec'+'.txt',forc_prec_nc_mean)
# print 'forc_us_nc_mean=',forc_us_nc_mean
# plt.imshow(forc_us_nc_mean,origin='lower left')
# plt.title('forc_us')
# plt.colorbar()
# plt.savefig('D:\\forc_us_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_vs_nc_mean=',forc_vs_nc_mean
# plt.imshow(forc_vs_nc_mean,origin='lower left')
# plt.title('forc_vs')
# plt.colorbar()
# plt.savefig('D:\\forc_vs_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_t_nc_mean=',forc_t_nc_mean
# plt.imshow(forc_t_nc_mean,origin='lower left')
# plt.title('forc_t')
# plt.colorbar()
# plt.savefig('D:\\forc_t_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_q_nc_mean=',forc_q_nc_mean
# plt.imshow(forc_q_nc_mean,origin='lower left')
# plt.title('forc_q')
# plt.colorbar()
# plt.savefig('D:\\forc_q_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_rh_nc_mean=',forc_rh_nc_mean
# plt.imshow(forc_rh_nc_mean,origin='lower left')
# plt.title('forc_rh')
# plt.colorbar()
# plt.savefig('D:\\forc_rh_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_prec_nc_mean=',forc_prec_nc_mean
# plt.imshow(forc_prec_nc_mean,origin='lower left')
# plt.title('forc_prec')
# plt.colorbar()
# plt.savefig('D:\\forc_prec_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_psrf_nc_mean=',forc_psrf_nc_mean
# plt.imshow(forc_psrf_nc_mean,origin='lower left')
# plt.title('forc_psrf')
# plt.colorbar()
# plt.savefig('D:\\forc_psrf_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_solarin_nc_mean=',forc_solarin_nc_mean
# plt.imshow(forc_solarin_nc_mean,origin='lower left')
# plt.title('forc_solarin')
# plt.colorbar()
# plt.savefig('D:\\forc_solarin_mean.jpg')  #save the figure
# plt.show()
# plt.clf()

# print 'forc_frl_nc_mean=',forc_frl_nc_mean
# plt.imshow(forc_frl_nc_mean,origin='lower left')
# plt.title('forc_frl')
# plt.colorbar()
# plt.savefig('D:\\forc_frl_mean.jpg')  #save the figure
# plt.show()
# plt.clf()
print 'program to finish running!'
