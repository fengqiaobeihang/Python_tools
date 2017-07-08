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
import  matplotlib.gridspec as gridspec

# imagepath='../image/'
# nc_path = 'I:/snow_data/'
outpath ='H:/UHB_WRF_DATA/Yakou_data/OUTPUT/result/'
elev1    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201401windspeed.txt')
elev2    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201402windspeed.txt')
elev3    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201403windspeed.txt')
elev4    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201404windspeed.txt')
elev5    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201405windspeed.txt')
elev6    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201406windspeed.txt')
elev7    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201407windspeed.txt')
elev8    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201408windspeed.txt')
elev9    = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201409windspeed.txt')
elev10   = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201410windspeed.txt')
elev11   = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201411windspeed.txt')
elev12   = np.loadtxt('H:/UHB_WRF_DATA/Yakou_data/OUTPUT/windspeed/201412windspeed.txt')
average=(elev1+elev2+elev3+elev4+elev5+elev6+elev7+elev8+elev9+elev10+elev11+elev12)/12
np.savetxt(outpath +'windspeed'+'.txt',average)