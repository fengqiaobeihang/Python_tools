'''
Created on 2015-10-15
Modified on 2016-7-5
@author: Hongyi LI
'''
import time
import numpy             as     np
import datetime          as     datetime
from   datetime          import timedelta
from   scipy.interpolate import griddata
import netCDF4           as     nc
from   netCDF4           import Dataset
from   osgeo             import gdal
from   osgeo.gdalconst   import *
import matplotlib.pylab  as     plt
# import pandas as pd
driver  = gdal.GetDriverByName('HFA')
driver.Register()
gdal.AllRegister()
############### old region of NC data #########################
# datapath = '/home/ecohydro/data/Cloud free snow products_China/amsrmts_txt/'
datapath = 'E:/Python_gdal/'
outpath = 'result/'
parapath= 'para/'
newnc   = Dataset(parapath+'2000_result.nc', format='NETCDF4')

lat_uhh = np.loadtxt('para/lat_upstreamHEIHE.txt',skiprows=6);lat_uhh=np.flipud(lat_uhh)
lon_uhh = np.loadtxt('para/lon_upstreamHEIHE.txt',skiprows=6);lon_uhh=np.flipud(lon_uhh)
ivt_uhh = np.loadtxt('para/ivt_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)


lat_nc  = lat_uhh
lon_nc  = lon_uhh

dims_lat,dims_lon = np.shape(lat_uhh)  

#-----------------------------------
iyear = 2003
fn      ='E:\\Python_gdal\\MODIS_Dysno_Cloudfree_20130509.tif'
#filename = datapath+str(iyear)+'/'+'MODIS_Dysno_Cloudfree_20130509.tif'
dataset = gdal.Open(fn, GA_ReadOnly)
band = dataset.GetRasterBand(1)

print band
print 'Driver: ', dataset.GetDriver().ShortName,'/', \
      dataset.GetDriver().LongName
print 'Size is ',dataset.RasterXSize,'x',dataset.RasterYSize, \
      'x',dataset.RasterCount
print 'Projection is ',dataset.GetProjection()

geotransform = dataset.GetGeoTransform()
if not geotransform is None:
    print 'Origin = (',geotransform[0], ',',geotransform[3],')'
    print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
#-----------------------------------

# print np.shape(lon_nc)

# old region infromation
# -----------------------------------
# ncols         15777
# nrows         7913
# xllcorner     59.875
# yllcorner     14.875
# cellsize      0.0050865689142938

ncols = dataset.RasterXSize
nrows = dataset.RasterYSize
resolution= geotransform[1]
xllcorner = geotransform[0]
yllcorner = geotransform[3]-resolution*nrows
bottom=yllcorner+resolution/2;top = bottom+resolution*(nrows-1)
left=xllcorner+resolution/2 ; right=left+resolution*(ncols-1)
# NC dimensions and its boundary coordination
xmin = left;   xmax = right
ymin = bottom; ymax = top


old_xi = np.arange(xmin,xmax,resolution)
old_yi = np.arange(ymin,ymax,resolution)
print np.shape(old_xi),np.shape(old_yi)
# determining the subset region covered DOI region
doi_ximin = np.abs(old_xi - np.min(lon_uhh)).argmin() 
doi_yimin = np.abs(old_yi - np.min(lat_uhh)).argmin()
doi_ximax = np.abs(old_xi - np.max(lon_uhh)).argmin()
doi_yimax = np.abs(old_yi - np.max(lat_uhh)).argmin()
if old_yi[doi_yimin] > np.min(lat_uhh): doi_yimin = doi_yimin - 1
if old_xi[doi_ximax] < np.max(lon_uhh): doi_ximax = doi_ximax + 1
if old_xi[doi_ximin] > np.min(lon_uhh): doi_ximin = doi_ximin - 1
if old_yi[doi_yimax] < np.max(lat_uhh): doi_yimax = doi_yimax + 1
old_xi        = old_xi[doi_ximin:doi_ximax+1]   
old_yi        = old_yi[doi_yimin:doi_yimax+1]   
old_xi,old_yi = np.meshgrid(old_xi, old_yi)     
old_points    = np.array([old_xi.flatten(), old_yi.flatten()]).T
doicols = np.arange(doi_ximin,doi_ximax+1)
data = band.ReadAsArray(0, 0, ncols, nrows)
data = np.array(data)
#data = np.flipud(data)
print 'data',data.shape,doi_ximin,doi_ximax+1,doi_yimin,doi_yimax+1
plt.imshow(data[doi_yimin:doi_yimax+1,doi_ximin:doi_ximax+1],origin='ll')
plt.show()
############### time ###########################
start_time=time.time()
first=datetime.date(2005,1,1)
last =datetime.date(2013,1,1)
leng =(last-first).days
print leng
dates=[first+datetime.timedelta(n) for n in range(leng+1)]
temp_year  = 0
temp_month = 0

for idate in dates:
    realdate = idate; print realdate
    iyear = realdate.year
    imonth = realdate.month
    iday = realdate.day
    date_delta = realdate - datetime.date(iyear,1,1) + timedelta(days = 1)
    idays = date_delta.days-1   
    print idays
    if temp_year != iyear or temp_month != imonth:
        newnc.close()
        newnc  = Dataset(outpath+'UHBSCA'+str(iyear)+'%02d'%(imonth)+'.nc', 'w', format='NETCDF4')
        temp_year = iyear
        temp_month= imonth
       
        # =======creat new NC file for resampled data========
        newnc.description = 'SCA in Upstream of HEIHE Basin'
        # dimensions    
        newnc.createDimension('time', None)   
        newnc.createDimension('lat', dims_lat)
        newnc.createDimension('lon', dims_lon)
        # variables                                            
        sca_nc  = newnc.createVariable('SCA'   , 'f4', ('time', 'lat', 'lon', )) 
        #END IF block

    scafile = np.loadtxt(datapath+'am'+str(iyear)+'%03d'%(idays+1)+'.asc',skiprows=6,usecols=doicols)
    scafile = np.flipud(scafile)
    # scafile = pd.read_csv(datapath+'am'+str(iyear)+'%03d'%(59)+'.asc',skiprows=6,index_col=False);scafile=np.array(scafile);scafile=np.flipud(scafile)
    scadoi = scafile[doi_yimin:doi_yimax+1,:]
    # print scafile
    # print np.shape(old_points),np.shape(scadoi)
    scadoi = scadoi.flatten()
    scadoi = np.array(scadoi)
    # print scadoi
    # points = np.array([lat_nc.flatten(),lon_nc.flatten()]).T
    #Griddata methods: 'nearest','linear','cubic'
    inew_sca_nc  = griddata(old_points, scadoi, (lon_uhh,lat_uhh),method='nearest')          
    sca_nc  [iday-1,:,:] = inew_sca_nc 
#END FOR loop
newnc.close()        
print 'end'   