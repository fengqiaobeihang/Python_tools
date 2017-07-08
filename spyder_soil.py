# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 10:58:01 2016

@author: lhy
"""
import time
import numpy             as     np
import datetime          as     datetime
import matplotlib.pylab  as     plt
import netCDF4           as     nc
import sys
# import gdal
import ogr
from   gdalconst         import *
from   osr               import SpatialReference
from   osgeo             import gdal, gdalconst
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
nc_path = 'F:/BaiduYunDownload/'
outpath = 'F:/UHB_soil_result/'
###############################################################
#read soil and estimate data
driver  = gdal.GetDriverByName('HFA')
driver.Register()
fn      ='F:/BaiduYunDownload/estimate_sm_1km_babao/SM_1km_2013213.tif'
ds      = gdal.Open(fn, GA_ReadOnly)
if ds is None:
   print 'Could not open ' + fn
   sys.exit(1)
cols         = ds.RasterXSize
rows         = ds.RasterYSize
band         = ds.GetRasterBand(1)
data         = band.ReadAsArray(0, 0, cols, rows)
print 'data',np.shape(data)
noDataValue  = band.GetNoDataValue()
projection   = ds.GetProjection()
geotransform = ds.GetGeoTransform()
print 'cols',cols
print 'rows',rows
# cols         =np.array(cols)
# rows         =np.array(rows)
# xList=[]
# yList=[]
x0=geotransform[0]  #originX
x1=geotransform[1]  #pixelWidth
x2=geotransform[2]
x3=geotransform[3]  #originY
x4=geotransform[4]
x5=geotransform[5]  #pixelHeight
print 'x0',x0
print 'x1',x1
print 'x2',x2
print 'x3',x3
print 'x4',x4
print 'x5',x5
ncols      =144
nrows      =84
xllcorner  =100.041666667
yllcorner  =38.3958333333
Xresolution=0.0083333333
Yresolution=-0.0083333333
left       =xllcorner; right=left+Xresolution*ncols
bottom     =yllcorner; top  =bottom+Yresolution*nrows
# dimensions and its boundary coordination
xmin = left+0.5*Xresolution;   xmax = right-0.5*Xresolution
ymin = bottom+0.5*Yresolution; ymax = top-0.5*Yresolution
xi    = np.arange(xmin,xmax,Xresolution)
yi    = np.arange(ymin,ymax,Yresolution)
print 'yi',np.shape(xi)
print 'xi',np.shape(yi)

xi,yi = np.meshgrid(xi, yi)
point =np.array([yi.flatten(),xi.flatten()]).T
print 'point',point
# for col in range(0,cols):
    # xTempList=[]
    # yTempList=[]
    # for row in range(0,rows):
    #     x    =geotransform[0]+col*geotransform[1]+row*geotransform[2]
    #     y    =geotransform[3]+col*geotransform[5]
    #     x    =np.array(x)
    #     y    =np.array(y)
        # x,y  =np.meshgrid(x,y)
        # xTempList.append(x)
        # yTempList.append(y)
        # point=np.array([y.flatten(),x.flatten()]).T
    # xList.append(xTempList)
    # yList.append(yTempList)
# xList,yList=np.meshgrid(xList,yList)
# point=np.array([xList.flatten(),yList.flatten()]).T
# x,y  =np.meshgrid(x,y)
# point=np.array([y.flatten(),x.flatten()]).T
# print 'xList',xList
# print 'yList',yList
# lat_nc        = ncfile.variables['lat']
# lon_nc        = ncfile.variables['lon']
# lat_nc        = np.array(lat_nc[2000:3000])
# lon_nc        = np.array(lon_nc[3000:4000])
# lon_nc,lat_nc = np.meshgrid(lon_nc,lat_nc)
# points        = np.array([lat_nc.flatten(),lon_nc.flatten()]).T
#for i in range(0, rows, yBSize): 
#if i + yBSize < rows: 
#        numRows = yBSize 
#else: 
#        numRows = rows – i 
#    for j in range(0, cols, xBSize): 
#        if j + xBSize < cols: 
#            numCols = xBSize 
#        else: 
#            numCols = colsnumCols = cols – j
#        data = band.ReadAsArray(j, i, numCols, numRows)
###############################################################
newcols = 226
newrows = 157

lat_uhh = np.loadtxt('F:/Python/data/ivt_upstreamHEIHE.txt',skiprows=6);lat_uhh=np.flipud(lat_uhh)
lon_uhh = np.loadtxt('F:/Python/data/lat_upstreamHEIHE.txt',skiprows=6);lon_uhh=np.flipud(lon_uhh)
ivt_uhh = np.loadtxt('F:/Python/data/lon_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)
plt.imshow(data,origin='lower left')
plt.colorbar()
plt.savefig('D:\\ivt_uhh.png')  #save the figure
plt.show()
plt.clf()

ncfile     = Dataset(nc_path+'CL.nc' , format='NETCDF4')
safile     = Dataset(nc_path+'SA.nc' , format='NETCDF4')
sifile     = Dataset(nc_path+'SI.nc' , format='NETCDF4')
somfile    = Dataset(nc_path+'SOM.nc', format='NETCDF4')
bdfile     = Dataset(nc_path+'BD.nc' , format='NETCDF4')
grav1file  = Dataset(nc_path+'GRAV1.nc' , format='NETCDF4')
grav2file  = Dataset(nc_path+'GRAV2.nc' , format='NETCDF4')
oc1file    = Dataset(nc_path+'OC1.nc' , format='NETCDF4')
oc2file    = Dataset(nc_path+'OC2.nc' , format='NETCDF4')
kschfile   = Dataset(nc_path+'K_SCH.nc' , format='NETCDF4')
ksvgfile   = Dataset(nc_path+'K_SVG.nc' , format='NETCDF4')
psisfile   = Dataset(nc_path+'PSI_S.nc' , format='NETCDF4')
alphafile  = Dataset(nc_path+'ALPHA.nc' , format='NETCDF4')
lfile      = Dataset(nc_path+'L.nc' , format='NETCDF4')
lambdafile = Dataset(nc_path+'LAMBDA.nc' , format='NETCDF4')
nfile      = Dataset(nc_path+'N.nc' , format='NETCDF4')
th33file   = Dataset(nc_path+'TH33.nc' , format='NETCDF4')
thrfile    = Dataset(nc_path+'THR.nc' , format='NETCDF4')
thschfile  = Dataset(nc_path+'THSCH.nc' , format='NETCDF4')
thsgmfile  = Dataset(nc_path+'THSGM.nc' , format='NETCDF4')

#read soil and estimate data
# smfile     = TIFF(nc_path+'SM_1km_2013213.tif' , format='libtiff')
# soilfile   = TIFF(nc_path+'20080123.tif' , format='libtiff')

lat_nc        = ncfile.variables['lat']
lon_nc        = ncfile.variables['lon']
lat_nc        = np.array(lat_nc[2000:3000])
lon_nc        = np.array(lon_nc[3000:4000])
lon_nc,lat_nc = np.meshgrid(lon_nc,lat_nc)
points        = np.array([lat_nc.flatten(),lon_nc.flatten()]).T
print 'lat_nc',lat_nc
print 'lon_nc',lon_nc
print 'points',points
#the range is different from the other dataset
lat_grav          = grav1file.variables['lat']
lon_grav          = grav1file.variables['lon']
lat_grav          = np.array(lat_grav[5000:6000])
lat_grav          = np.flipud(lat_grav)
lon_grav          = np.array(lon_grav[33000:34000])
lon_grav,lat_grav = np.meshgrid(lon_grav,lat_grav)
points_grav       = np.array([lat_grav.flatten(),lon_grav.flatten()]).T

old_CL            = ncfile.variables['CL']
old_SA            = safile.variables['SA']
old_SI            = sifile.variables['SI']
old_BD            = bdfile.variables['BD']
old_SOM           = somfile.variables['SOM']
grav1             = grav1file.variables['GRAV']
grav2             = grav2file.variables['GRAV']
oc1               = oc1file.variables['OC']
oc2               = oc2file.variables['OC']
ksch              = kschfile.variables['K_SCH']
ksvg              = ksvgfile.variables['K_SVG']
psis              = psisfile.variables['PSI_S']
depth             = safile.variables['depth']
thsch =thschfile .variables['THSCH'] 
l     =lfile.variables['L']  
alpha =alphafile .variables['ALPHA']    
lambd =lambdafile.variables['LAMBDA']
n     =nfile     .variables['N']     
th33  =th33file  .variables['TH33']  
thr   =thrfile   .variables['THR']   
thsgm =thsgmfile .variables['THSGM'] 

# =======creat new NC file for resampled data========
newnc  = Dataset(nc_path+'SOIL_VEG_UHB.nc', 'w', format='NETCDF4')
newnc.description = 'Upstream of HEIHE Basin'
# dimensions    
newnc.createDimension('depth', None)   
newnc.createDimension('LAT', newrows)
newnc.createDimension('LON', newcols)
# variables   
                                         
new_CL    = newnc.createVariable('CL'   , 'f4', ('depth', 'LAT', 'LON', )) 
new_SA    = newnc.createVariable('SA'   , 'f4', ('depth', 'LAT', 'LON', )) 
new_SI    = newnc.createVariable('SI'   , 'f4', ('depth', 'LAT', 'LON', )) 
new_BD    = newnc.createVariable('BD'   , 'f4', ('depth', 'LAT', 'LON', ))
new_SOM   = newnc.createVariable('SOM'  , 'f4', ('depth', 'LAT', 'LON', )) 
new_grav  = newnc.createVariable('GRAV' , 'f4', ('depth', 'LAT', 'LON', )) 
new_oc    = newnc.createVariable('OC'   , 'f4', ('depth', 'LAT', 'LON', )) 
new_ivt   = newnc.createVariable('ivt'  , 'f4', ( 'LAT', 'LON', ))
new_ksch  = newnc.createVariable('KSCH' , 'f4', ('depth', 'LAT', 'LON', )) 
new_ksvg  = newnc.createVariable('KSVG' , 'f4', ('depth', 'LAT', 'LON', )) 
new_psis  = newnc.createVariable('PSIS' , 'f4', ('depth', 'LAT', 'LON', )) 
new_alpha = newnc.createVariable('ALPHA', 'f4', ('depth', 'LAT', 'LON', ))
new_l     = newnc.createVariable('L', 'f4', ('depth', 'LAT', 'LON', ))
new_lambd = newnc.createVariable('LAMBDA', 'f4', ('depth', 'LAT', 'LON', ))
new_n     = newnc.createVariable('N', 'f4', ('depth', 'LAT', 'LON', ))
new_th33  = newnc.createVariable('TH33', 'f4', ('depth', 'LAT', 'LON', ))
new_thr   = newnc.createVariable('THR', 'f4', ('depth', 'LAT', 'LON', ))
new_thsch = newnc.createVariable('THSCH', 'f4', ('depth', 'LAT', 'LON', ))
new_thsgm = newnc.createVariable('THSGM', 'f4', ('depth', 'LAT', 'LON', ))

new_SM    = newnc.createVariable('SM', 'f4', ('depth', 'LAT', 'LON', ))#estimate_sm_1km_babao

varnames = ['alpha','l' ,'lambd','n','th33','thr','thsch','thsgm']

for i in range(7):
    print i

    iold_CL_nc    = old_CL [i,2000:3000,3000:4000]
    iold_SA_nc    = old_SA [i,2000:3000,3000:4000]
    iold_SI_nc    = old_SI [i,2000:3000,3000:4000]
    iold_BD_nc    = old_BD [i,2000:3000,3000:4000]
    iold_SOM_nc   = old_SOM[i,2000:3000,3000:4000]
    iold_ksch_nc  = ksch[i,2000:3000,3000:4000]
    iold_ksvg_nc  = ksvg[i,2000:3000,3000:4000]
    iold_psis_nc  = psis[i,2000:3000,3000:4000]
    
    iold_alpha_nc =alpha[i,2000:3000,3000:4000]
    iold_l_nc     =l    [i,2000:3000,3000:4000]
    iold_lambd_nc =lambd[i,2000:3000,3000:4000]
    iold_n_nc     =n    [i,2000:3000,3000:4000]
    iold_th33_nc  =th33 [i,2000:3000,3000:4000]
    iold_thr_nc   =thr  [i,2000:3000,3000:4000]
    iold_thsch_nc =thsch[i,2000:3000,3000:4000]
    iold_thsgm_nc =thsgm[i,2000:3000,3000:4000]

    #iold_SM_nc    = data[i,cols,rows]
    #Griddata methods: 'nearest','linear','cubic'
    tmp_lat = lat_nc.flatten()
    tmp_lon = lon_nc.flatten()    
    for varname in varnames:
        # tmp_iold_var = str('tmp_iold_'+varname)
        iold_var_nc = eval('iold_'+varname+'_nc')
        new_var = eval('new_'+varname)
        tmp_iold_var = iold_var_nc.flatten()
        mask_iold_var=tmp_iold_var[tmp_iold_var>-9900]
        mask_lat = tmp_lat[tmp_iold_var>-9900]
        mask_lon = tmp_lon[tmp_iold_var>-9900]
        mask_points  = np.array([mask_lat,mask_lon]).T
        mask_var_nc = griddata(mask_points, mask_iold_var, (lat_uhh,lon_uhh ),method='nearest') 
        if i==6: 
            inew_var_nc  = griddata(points, tmp_iold_var, (lat_uhh,lon_uhh ),method='nearest') 
            inew_var_nc  = np.array(inew_var_nc) 
            inew_var_nc = np.where(inew_var_nc<-9000,np.mean(mask_iold_var),inew_var_nc)
        else:    
            inew_var_nc = mask_var_nc
        new_var[i,:,:] = inew_var_nc

    inew_psis_nc  = griddata(points, iold_psis_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_psis_nc  = np.array(inew_psis_nc) 
    tmp_iold_psis = iold_psis_nc.flatten()
    mask_iold_psis=tmp_iold_psis[tmp_iold_psis>-9900]
    tmp_lat = lat_nc.flatten()
    tmp_lon = lon_nc.flatten()
    mask_lat = tmp_lat[tmp_iold_psis>-9900]
    mask_lon = tmp_lon[tmp_iold_psis>-9900]
    mask_points  = np.array([mask_lat,mask_lon]).T
    mask_psis_nc = griddata(mask_points, mask_iold_psis, (lat_uhh,lon_uhh ),method='nearest') 
    if i==6: 
        inew_psis_nc = np.where(inew_psis_nc<-9000,np.mean(mask_iold_psis),inew_psis_nc)
    else:
        inew_psis_nc = np.where(inew_psis_nc<-9000,mask_psis_nc,inew_psis_nc)
    new_psis  [i,:,:] = inew_psis_nc

    inew_ksch_nc  = griddata(points, iold_ksch_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_ksch_nc  = np.array(inew_ksch_nc)  
    tmp_iold_ksch = iold_ksch_nc.flatten()
    mask_iold_ksch=tmp_iold_ksch[tmp_iold_ksch>-9900]
    tmp_lat = lat_nc.flatten()
    tmp_lon = lon_nc.flatten()
    mask_lat = tmp_lat[tmp_iold_ksch>-9900]
    mask_lon = tmp_lon[tmp_iold_ksch>-9900]
    mask_points  = np.array([mask_lat,mask_lon]).T
    mask_ksch_nc = griddata(mask_points, mask_iold_ksch, (lat_uhh,lon_uhh ),method='nearest')
    if i==6: 
        inew_ksch_nc = np.where(inew_ksch_nc<-9000,np.mean(mask_iold_ksch),inew_ksch_nc)
    else:
        inew_ksch_nc = np.where(inew_ksch_nc<-9000,mask_ksch_nc,inew_ksch_nc)
    new_ksch  [i,:,:] = inew_ksch_nc



    inew_ksvg_nc  = griddata(points, iold_ksvg_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_ksvg_nc  = np.array(inew_ksvg_nc)  
    tmp_iold_ksvg = iold_ksvg_nc.flatten()
    mask_iold_ksvg=tmp_iold_ksvg[tmp_iold_ksvg>-9900]
    tmp_lat = lat_nc.flatten()
    tmp_lon = lon_nc.flatten()
    mask_lat = tmp_lat[tmp_iold_ksvg>-9900]
    mask_lon = tmp_lon[tmp_iold_ksvg>-9900]
    mask_points  = np.array([mask_lat,mask_lon]).T
    mask_ksvg_nc = griddata(mask_points, mask_iold_ksvg, (lat_uhh,lon_uhh ),method='nearest') 
    if i==6: 
        inew_ksvg_nc = np.where(inew_ksvg_nc<-9000,np.mean(mask_iold_ksvg),inew_ksvg_nc)
    else:    
        inew_ksvg_nc = np.where(inew_ksvg_nc<-9000,mask_ksvg_nc,inew_ksvg_nc)
    new_ksvg  [i,:,:] = inew_ksvg_nc
     
    print np.shape(points),np.shape(iold_CL_nc .flatten())

    inew_CL_nc  = griddata(points, iold_CL_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_CL_nc  = np.array(inew_CL_nc)  
    inew_CL_nc[inew_CL_nc<0]   = 30.
    inew_CL_nc[inew_CL_nc>100.]= 30.
    new_CL  [i,:,:] = inew_CL_nc
    # print np.shape(points),np.shape(iold_SA_nc .flatten())
    inew_SA_nc  = griddata(points, iold_SA_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_SA_nc  = np.array(inew_SA_nc)  
    inew_SA_nc[inew_SA_nc<0]   = 30.
    inew_SA_nc[inew_SA_nc>100.]= 30.
    new_SA  [i,:,:] = inew_SA_nc

    inew_SI_nc  = griddata(points, iold_SI_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_SI_nc  = np.array(inew_SI_nc)  
    inew_SI_nc[inew_SI_nc<0]   = 30.
    inew_SI_nc[inew_SI_nc>100.]= 30.
    new_SI  [i,:,:] = inew_SI_nc
    
    inew_SOM_nc  = griddata(points, iold_SOM_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_SOM_nc  = np.array(inew_SOM_nc)  
    inew_SOM_nc[inew_SOM_nc<0]   = 30.
    inew_SOM_nc[inew_SOM_nc>100.]= 30.
    new_SOM  [i,:,:] = inew_SOM_nc

    inew_BD_nc  = griddata(points, iold_BD_nc .flatten(), (lat_uhh,lon_uhh ),method='nearest') 
    inew_BD_nc  = np.array(inew_BD_nc)  
    inew_BD_nc[inew_BD_nc<0]   = 1.
    inew_BD_nc[inew_BD_nc>100.]= 1.
    new_BD  [i,:,:] = inew_BD_nc    

    print np.shape(point),np.shape(data.flatten())

    inew_SM_nc  = griddata(point, data.flatten(), (lat_uhh,lon_uhh ),method='linear') 
    inew_SM_nc  = np.array(inew_SM_nc)  
    # inew_SM_nc[inew_SM_nc<0]   = 0.1
    # inew_SM_nc[inew_SM_nc>1.]  = 0.5
    new_SM  [i,:,:] = inew_SM_nc

    if i<= 3:
        iold_grav1  = grav1[i,5000:6000,33000:34000]
        iold_grav1  = np.flipud(iold_grav1)    
        inew_grav   = griddata(points_grav, iold_grav1.flatten(), (lat_uhh,lon_uhh ),method='nearest') 
        inew_grav   = np.array(inew_grav)  
        # inew_BD_nc[inew_BD_nc<0]   = 1.
        # inew_BD_nc[inew_BD_nc>100.]= 1.
        new_grav  [i,:,:] = inew_grav    

        iold_oc1  = oc1[i,5000:6000,33000:34000]*0.01
        iold_oc1  = np.flipud(iold_oc1)    
        inew_oc   = griddata(points_grav, iold_oc1.flatten(), (lat_uhh,lon_uhh ),method='nearest') 
        inew_oc   = np.array(inew_oc)  
        # inew_BD_nc[inew_BD_nc<0]   = 1.
        # inew_BD_nc[inew_BD_nc>100.]= 1.
        new_oc  [i,:,:] = inew_oc    
    else:
        iold_grav2  = grav2[i-4,5000:6000,33000:34000]
        iold_grav2  = np.flipud(iold_grav2)
        inew_grav   = griddata(points_grav, iold_grav2.flatten(), (lat_uhh,lon_uhh ),method='nearest') 
        inew_grav   = np.array(inew_grav)  
        # inew_BD_nc[inew_BD_nc<0]   = 1.
        # inew_BD_nc[inew_BD_nc>100.]= 1.
        new_grav  [i,:,:] = inew_grav   
        
        iold_oc2  = oc2[i-4,5000:6000,33000:34000]*0.01
        iold_oc2  = np.flipud(iold_oc2)    
        inew_oc   = griddata(points_grav, iold_oc2.flatten(), (lat_uhh,lon_uhh ),method='nearest') 
        inew_oc   = np.array(inew_oc)  
        # inew_BD_nc[inew_BD_nc<0]   = 1.
        # inew_BD_nc[inew_BD_nc>100.]= 1.
        new_oc  [i,:,:] = inew_oc    
    # plt.imshow(inew_grav,origin='ll')
    # plt.show()
    #save to TXT format also.
    np.savetxt(outpath +'CL_UHB'  +str(i)+'.txt',inew_CL_nc )
    np.savetxt(outpath +'SI_UHB'  +str(i)+'.txt',inew_SI_nc )
    np.savetxt(outpath +'SA_UHB'  +str(i)+'.txt',inew_SA_nc )
    np.savetxt(outpath +'BD_UHB'  +str(i)+'.txt',inew_BD_nc )
    np.savetxt(outpath +'SOM_UHB' +str(i)+'.txt',inew_SOM_nc)
    np.savetxt(outpath +'GRAV_UHB'+str(i)+'.txt',inew_grav  )
    np.savetxt(outpath +'OC_UHB'  +str(i)+'.txt',inew_oc    )
#END FOR loop
#---------transfer vegetation classification to USGS standard
# 1	寒温性针叶林	Needleleaf  forest -> 14
# 2	灌丛	        Shrub              ->  8
# 3	农田	        Farmland           ->  5
# 4	草原	        Alpine steppe      ->  7
# 5	高寒草甸	    Alpine meadow      ->  7
# 6	高山苔原	    Alpine Tundra      -> 23
# 7	冰川	        Glacier            -> 24
# 8	城市与村镇	City and town      ->  1
# 9	水体	        Water body         -> 16
# 10祼地	        Desert             -> 23
#USGS_CLASSIFICATION
#---------------------------
# GLCC USGS Land Use/Land Cover System Legend 
# 0  Ocean
# 1  Urban and Built-Up Land
# 2  Dryland Cropland and Pasture
# 3  Irrigated Cropland and Pasture
# 4  Mixed Dryland/Irrigated Cropland and Pasture
# 5  Cropland/Grassland Mosaic
# 6  Cropland/Woodland Mosaic
# 7  Grassland
# 8  Shrubland
# 9  Mixed Shrubland/Grassland
#10  Savanna
#11  Deciduous Broadleaf Forest 
#12  Deciduous Needleleaf Forest 
#13  Evergreen Broadleaf Forest
#14  Evergreen Needleleaf Forest
#15  Mixed Forest
#16  Inland Water
#17  Herbaceous Wetland
#18  Wooded Wetland
#19  Barren or Sparsely Vegetated
#20  Herbaceous Tundra
#21  Wooded Tundra
#22  Mixed Tundra
#23  Bare Ground Tundra
#24  Snow or Ice
# ivt = np.array(ivt_uhh)
# ivt[ivt == 1] = 14+1000
# ivt[ivt == 2] =  8+1000
# ivt[ivt == 3] =  5+1000
# ivt[ivt == 4] =  7+1000
# ivt[ivt == 5] =  7+1000
# ivt[ivt == 6] = 23+1000
# ivt[ivt == 7] = 24+1000
# ivt[ivt == 8] =  1+1000
# ivt[ivt == 9] = 16+1000
# ivt[ivt == 10]= 23+1000
# ivt[ivt == 8] =  1+1000
# ivt           = ivt-1000
# ivt[ivt<0.]   = -1
# new_ivt[:,:]  = np.int8(ivt)
# np.savetxt(outpath +'ivt_UHB.txt',np.int8(ivt) )

newnc    .close() 
ncfile   .close()       
safile   .close() 
sifile   .close() 
bdfile   .close() 
somfile  .close()
grav1file.close() 
grav2file.close() 
oc1file  .close() 
oc2file  .close() 
alphafile .close()
lfile     .close()
lambdafile.close()
nfile     .close()
th33file  .close()
thrfile   .close()
thschfile .close()
thsgmfile .close()
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