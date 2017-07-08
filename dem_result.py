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

imagepath='../image/'
nc_path = 'I:/snow_data/'
outpath ='I:/snow_data/output/'
nyear = 12.
elev    = np.loadtxt('C:/Users/sdh/Desktop/WRF_verify/para/dem_uhh.txt')
lat_uhh = np.loadtxt('C:/Users/sdh/Desktop/WRF_verify/para/lat_upstreamHEIHE.txt',skiprows=6);lat_uhh=np.flipud(lat_uhh)
lon_uhh = np.loadtxt('C:/Users/sdh/Desktop/WRF_verify/para/lon_upstreamHEIHE.txt',skiprows=6);lon_uhh=np.flipud(lon_uhh)
ivt_uhh = np.loadtxt('C:/Users/sdh/Desktop/WRF_verify/para/ivt_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)
ks1     = np.loadtxt('C:/Users/sdh/Desktop/WRF_verify/para/3900-4700.asc',skiprows=6); ks1=np.flipud(ks1)
mask = np.where(ks1>0,1,0)
ncresult= Dataset(nc_path+'2000_result.nc', format='NETCDF4')
xmin = np.min(lon_uhh[lon_uhh>0])
xmax = np.max(lon_uhh[lon_uhh>0])
ymin = np.min(lat_uhh[lon_uhh>0])
ymax = np.max(lat_uhh[lon_uhh>0])

yearlist = np.arange(2001,2016)
hrank= np.arange(1670,4772,100)
hrank= np.arange(2250,4472,100)
hrank= np.append(hrank,4800)
# print hrank
elevlist = []
for ihrank in hrank[0:-1]:
    elevlist.append(str(ihrank)+'m-'+str(ihrank+100)+'m')
elevlist[-1]='4450m-4800m'
# print elevlist
# print np.sum(maskelv),np.min(elev),np.max(elev)
# ymin=37.8;ymax=39.2

color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                  '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
                  '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
                  '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
color_sequence = ['dimgray','skyblue','coral','peru','red','khaki',
                    'y','olive','green','yellowgreen',
                    'springgreen','navy','c','dodgerblue','violet','magenta']
font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 14,
        }
avgsevap = lat_uhh*0.
avgsmelt = lat_uhh*0.
avgsfall = lat_uhh*0.
avgswe   = lat_uhh*0.

for iyear in range(2001,2016):
    print iyear
    ncresult.close()
    ncresult  = Dataset(nc_path+str(iyear)+'_result.nc', format='NETCDF4')
    sevap = ncresult.variables['annual_sevap']
    smelt = ncresult.variables['annual_smelt']
    sfall = ncresult.variables['annual_snow']
    swe   = ncresult.variables['annual_swe']
    airt  = ncresult.variables['annual_temp']
    avgsevap = avgsevap + sevap[0,:,:]/nyear
    avgsmelt = avgsmelt + smelt[0,:,:]/nyear
    avgsfall = avgsfall + sfall[0,:,:]/nyear
    avgswe   = avgswe   + swe[0,:,:]/nyear

    tmp_melt = []
    tmp_airt = []
    tmp_fall = []
    tmp_swe  = []
    # print ih,hrank[np.argmin((hrank-ih)**2)+1]
    maskelv = mask*(elev>=1669)*(elev<2300)
    # tmp_airt.append(np.sum(smelt[0,:,:]*maskelv)/np.sum(sfall[0,:,:]*maskelv))
    tmp_melt.append(np.sum(smelt[0,:,:]*maskelv)/np.sum(maskelv))
    tmp_airt.append(np.sum(airt[0,:,:]*maskelv)/np.sum(maskelv))
    tmp_fall.append(np.sum(sfall[0,:,:]*maskelv)/np.sum(maskelv))
    tmp_swe.append(np.sum(swe[0,:,:]*maskelv)/np.sum(maskelv))
    # print 'ih=',ih,np.sum(smelt[0,:,:]*maskelv),np.sum(maskelv)
    tmp_airt = np.array(tmp_airt).reshape(1,np.size(tmp_airt))
    tmp_fall = np.array(tmp_fall).reshape(1,np.size(tmp_fall))
    tmp_melt = np.array(tmp_melt).reshape(1,np.size(tmp_melt))
    tmp_swe  = np.array(tmp_swe).reshape(1,np.size(tmp_swe))
    # print tmp_airt
    if iyear == 2001:
        ann_airt = tmp_airt
        ann_melt = tmp_melt
        ann_fall = tmp_fall
        ann_swe  = tmp_swe
    else:
        ann_airt = np.append(ann_airt,tmp_airt,axis=0)
        ann_melt = np.append(ann_melt,tmp_melt,axis=0)
        ann_fall = np.append(ann_fall,tmp_fall,axis=0)
        ann_swe  = np.append(ann_swe,tmp_swe,axis=0)
print 'ann_airt=',ann_airt
print 'ann_melt=',ann_melt
print 'ann_fall=',ann_fall
print 'ann_swe=' ,ann_swe