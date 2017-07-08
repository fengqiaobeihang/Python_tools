# -*- coding: utf-8 -*-
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
import  os

data   = np.loadtxt('C:/Users/sdh/Desktop/data1.txt') 
line   = np.shape(data)
database=np.array([[row[i] for i in range(0,19) if i !=0] for row in data])
print 'database=',database# print 'data=',data
# print 'line=',line
# data=array(data)
# print 'data=',data
# for i in range(10):
# 	for j in range(19):
# 		database=data[1,j]
#         print 'database=',database