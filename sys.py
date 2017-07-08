import sys
import gdal
from gdalconst import *
from osr import SpatialReference

fn = r'F:\BaiduYunDownload\estimate_sm_1km_babao\SM_1km_2013213.tif'
print fn
ds = gdal.Open(fn,GA_ReadOnly)
if ds is None:
    print 'cannot open ',fn
    sys.exit(1)
print 'size is ',ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount
print 'projection is:',ds.GetProjection()
geotransform = ds.GetGeoTransform()
if not geotransform is None:
    print 'origin(x,y) is:',geotransform[0],',',geotransform[3]
    print 'Pixel size is:',geotransform[1],',',geotransform[5]
band = ds.GetRasterBand(1)
print 'Band Type=',gdal.GetDataTypeName(band.DataType)
data=band.ReadAsArray()
print 'data is ',data
min = band.GetMinimum()
max = band.GetMaximum()
print 'Nodatavalue is:',band.GetNoDataValue()
if min is None or max is None:
    (min,max)=band.ComputeRasterMinMax(1)
print 'Min=%.3f,Max=%.3f' % (min,max)
if band.GetOverviewCount() > 0:
    print 'Band has', band.GetOverviewCount(), 'overviews.'
if not band.GetRasterColorTable() is None:
    print 'Band has a color table with ',band.GetRasterColorTable().GetCount(),\
          'entries'