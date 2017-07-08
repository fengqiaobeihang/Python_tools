import arcpy
from arcpy import env
from arcpy.sa import *
import os
import os.path
import sys

arcpy.env.workspace = "H:\\2014MOD09GA006\\MOD09GA_angle\\angle"
rootdir = "H:\\2014MOD09GA006\\MOD09GA_angle\\angle"
for dirpath,filename,filenames in os.walk(rootdir):
    for filename in filenames:
        if os.path.splitext(filename)[1] == '.img':
            filepath = os.path.join(dirpath,filename)
            inASCII = filepath
            outname = filename.replace('.img','')
            print "ok"
            print filepath
            print outname
            outRaster= "H:\\2014MOD09GA006\\angle_fin\\"+outname+".tif"
            arcpy.Resample_management(filepath, outRaster, "463.312717", "NEAREST")

