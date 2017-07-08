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

mode     = 1 # 0:windows OS; 1:linux OS
run_pbsm = 0
run_gbhm = 0
run_plot = 0
plot_img = 0
record   = 1
builtin  = 1
builtmod = 1
if mode == 0:
    nc_path = 'H:/UHB_WRF_DATA/WRF/'
    out_path= 'H:/UHB_WRF_DATA/Result/'
    # plt.ion()
if mode == 1:
    modis_path             = 'H:/UHB_WRF_DATA/WRF/'
    out_path_Temperature   = 'H:/UHB_WRF_DATA/Result/Temperature/'
    out_path_Windspeed     = 'H:/UHB_WRF_DATA/Result/Windspeed/'
    out_path_Precipitation = 'H:/UHB_WRF_DATA/Result/Precipitation/'

modisdata= Dataset(modis_path+'UHB200001.nc', format='NETCDF4')
newcols  = 226
newrows  = 157

lon_points  = newrows  # number of lat points, noted there is an error to be corrected
lat_points  = newcols  # number of lon points

############### time ###########################
start_time = time.time()
first      = datetime.date(2000,1,1)
last       = datetime.date(2014,1,1)
leng       = (last-first). days
print leng 
dates      = [first+datetime.timedelta(n) for n in range(leng+1)]
temp_year  =  0
temp_month =  0

ivt_uhh    = np.loadtxt('H:/SCA_DATA/data/ivt_upstreamHEIHE.txt',skiprows=6);ivt_uhh=np.flipud(ivt_uhh)
mask       = ivt_uhh
mask       = np.where(ivt_uhh>-1,1,0)

# startyear = first.year
# endyear   = last.year
# startmonth= first.month
# startday  = first.day 
# endmonth  = last.month
# endday    = last.day

asca                 = np.zeros((157,226) )  
msca                 = np.zeros((157,226) )  
testrnf              = np.zeros((157,226) )  
testp                = np.zeros((157,226) ) 
annual_windspeed     = np.zeros((157,226) )
annual_Precipitation = np.zeros((157,226) )
# dsca    = [1,2,2]
sum_rssca1=[0]
for i_date in dates:
    realdate   = i_date
    iyear      = realdate.year
    imonth     = realdate.month
    iday       = realdate.day
    idate      = np.array([i_date.year,i_date.month,i_date.day])
    date_delta = realdate - datetime.date(iyear,1,1) + timedelta(days = 1)
    idays      = date_delta.days-1   
    idc        = date_delta.days  #contineous day in a year (1-366)
    ihc0       = (idc-1)*24
    ihc        = ihc0

    if temp_year != iyear:
        # dsca = np.array(dsca)
        # print 'year:',iyear-1, np.mean(dsca),'var:',np.std(dsca),'max:',np.max(dsca)
        # np.savetxt(str(iyear)+'sca.txt',asca)
        ########################Statistics snow cover days
        # print 'asca=',asca
        mean_asca=asca/365
        # print 'mean_asca',mean_asca
        plt.imshow(mean_asca,origin='lower left')
        plt.title(str(iyear-1)+' '+'Temperature')
        plt.colorbar()
        figure_path='H:/UHB_WRF_DATA/figure/T2/'
        plt.savefig(figure_path+str(iyear-1)+' '+'Temperature.jpg')  #save the figure
        plt.show()
        plt.clf()
        ##################################################
        #风速画图
        mean_windspeed=annual_windspeed/365
        plt.imshow(mean_windspeed,origin='lower left')
        plt.title(str(iyear-1)+' '+'Windspeed')
        plt.colorbar()
        figure_path='H:/UHB_WRF_DATA/figure/Wind speed/'
        plt.savefig(figure_path+str(iyear-1)+' '+'Windspeed.jpg')  #save the figure
        plt.show()
        plt.clf()
        ##################################################
         #降水量画图
        mean_Precipitation=annual_Precipitation/365
        plt.imshow(mean_Precipitation,origin='lower left')
        plt.title(str(iyear-1)+' '+'Precipitation')
        plt.colorbar()
        figure_path='H:/UHB_WRF_DATA/figure/PREC/'
        plt.savefig(figure_path+str(iyear-1)+' '+'Precipitation.jpg')  #save the figure
        plt.show()
        plt.clf()
        ##################################################       
        asca                 = mask*0.
        annual_windspeed     = mask*0.
        annual_Precipitation = mask*0.
        # dsca = []
    if temp_year != iyear or temp_month != imonth:        
        if i_date>=datetime.date(2000,1,1) and i_date<datetime.date(2014,1,1):
            modisdata.close()
            modisdata  = Dataset(modis_path+'UHB'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
            #read data
            PREC_dataset = modisdata.variables['PREC']#针对降水量数据需要累加
            T2_dataset   = modisdata.variables['T2']
            U10_dataset  = modisdata.variables['U10']
            V10_dataset  = modisdata.variables['V10']
            # windspeed    = (U10_dataset**2+V10_dataset**2)**0.5
            # windspeed=pow((pow(U10_dataset,2)+pow(U10_dataset,2)),0.05)
        # if i_date>=datetime.date(2004,4,1) and i_date<datetime.date(2015,1,1):
        #     modisdata.close()
        #     modisdata  = Dataset(modis_path+'UHB'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
        #     modisdataset = modisdata.variables['SCA']          
        temp_year  = iyear
        temp_month = imonth  
    #日均温度
    Temperature        = T2_dataset[iday-1,:,:]
    Yakou_Temperature  = Temperature[124,144]#取出第89行150的数据
    # print 'Yakou_Temperature=',Yakou_Temperature
    write_file=open(out_path_Temperature +'Yakou_Temperature'+str(iyear)+'.txt','a')
    write_file.write(str(Yakou_Temperature)+'\n')
    # np.savetxt(str(iyear)+'sca.txt',asca)
    # np.savetxt(out_path +'Yakou_Temperature'+str(iyear)+'.txt',str(Yakou_Temperature)+'\t')
    # rssca  = np.where(rssca1==200,1.,0.)#取出所有值为200的数据
    # rssca1  = rssca1*mask
    #统计年温度
    asca   = asca + Temperature
    # dsca.append(np.sum(rssca1))
    #日均风速
    U10              = U10_dataset[iday-1,:,:]
    V10              = V10_dataset[iday-1,:,:]
    windspeed        = (U10**2+V10**2)**0.5
    Yakou_windspeed  = windspeed[124,144]
    write_file=open(out_path_Windspeed +'Yakou_windspeed'+str(iyear)+'.txt','a')
    write_file.write(str(Yakou_windspeed)+'\n')
    # print 'Yakou_windspeed=',Yakou_windspeed
    # np.savetxt(out_path +'Yakou_windspeed.txt',Yakou_windspeed)
    annual_windspeed = annual_windspeed + windspeed
    #日均降水量
    Precipitation       = PREC_dataset[iday-1,:,:]
    Yakou_Precipitation = Precipitation[124,144]
    write_file=open(out_path_Precipitation +'Yakou_Precipitation'+str(iyear)+'.txt','a')
    write_file.write(str(Yakou_Precipitation)+'\n')
    # print 'Yakou_Precipitation',Yakou_Precipitation
    # np.savetxt(out_path +'Yakou_Precipitation.txt',Yakou_Precipitation)
    annual_Precipitation= annual_Precipitation + Precipitation
#################figure######################################
# plt.imshow(asca,origin='lower left')
# plt.title('asca')
# plt.colorbar()
# plt.savefig('D:\\asca.jpg')  #save the figure
# plt.show()
# plt.clf()
print 'The program is end!'
modisdata.close()