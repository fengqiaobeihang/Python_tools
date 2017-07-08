import  numpy              as     np
from    numpy              import array
from    random             import random
from    math               import sin, sqrt
import  time          
import  datetime           as     datetime
from    datetime           import timedelta
import  netCDF4            as     nc
from    netCDF4            import Dataset
from    decimal            import *
import  matplotlib.pylab   as     plt
from    matplotlib.pyplot  import figure, show, cm

mode     = 1 # 0:windows OS; 1:linux OS
run_pbsm = 0
run_gbhm = 0
run_plot = 0
plot_img = 0
record   = 1
builtin  = 1
builtmod = 1
if mode == 0:
    nc_path = 'H:/SCA_DATA/resample_SCA/result/'
    out_path= 'E:/result_UHB3/'
    import matplotlib.pylab as plt
    # plt.ion()
if mode == 1:
    modis_path  = 'H:/SCA_DATA/resample_SCA/result/'

modisdata= Dataset(modis_path+'UHBSCA200801.nc', format='NETCDF4')
newcols  = 226
newrows  = 157

lon_points  = newrows  # number of lat points, noted there is an error to be corrected
lat_points  = newcols  # number of lon points

############### time ###########################
start_time = time.time()
first      = datetime.date(2000,3,1)
last       = datetime.date(2015,1,1)
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

asca    = np.zeros((157,226) )  
msca    = np.zeros((157,226) )  
testrnf = np.zeros((157,226) )  
testp   =  np.zeros((157,226) ) 
dsca    = [1,2,2]
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
        dsca = np.array(dsca)
        print 'year:',iyear-1, np.mean(dsca),'var:',np.std(dsca),'max:',np.max(dsca)
        np.savetxt(str(iyear)+'sca.txt',asca)
        ########################Statistics snow cover days
        # import matplotlib as mpl
        # cmap = mpl.cm.gray_r
        # norm = mpl.colors.Normalize(vmin=0)
        # plt.imshow(pdata,cmap=cmap)
        ###########################################
        # from matplotlib.pyplot import figure, show, savefig
        # from matplotlib import cm, colors
        # from numpy import ma
        # fig   = figure()
        # ax    = fig.add_subplot(111)
        # # ax.set_axis_bgcolor("#bdb76b")
        # ax.set_axis_bgcolor('white')
        # cs=ax.pcolormesh(asca, shading='gouraud')
        # cs=ax.imshow(asca,origin='lower left',cmap=cm.jet)
        # fig.colorbar(cs)
        # show()
        ###########################################
        plt.imshow(asca,origin='lower left')
        plt.title(str(iyear)+' '+'Snow Cover Days')
        plt.colorbar(orientation='horizontal')
        plt.clim(vmin=0,vmax=365)
        figure_path='H:/SCA_DATA/resample_SCA/snow cover days/'
        plt.savefig(figure_path+str(iyear)+' '+'Snow Cover Days.jpg')  #save the figure
        plt.show()
        plt.clf()
        #############################################
        asca = mask*0.
        dsca = []
    if temp_year != iyear or temp_month != imonth:        
        if i_date>=datetime.date(2000,3,1) and i_date<datetime.date(2004,4,1):
            modisdata.close()
            modisdata  = Dataset(modis_path+'UHBSCA'+str(iyear)+'%02d'%(imonth)+'.nc', format='NETCDF4')
            modisdataset = modisdata.variables['SCA']
        if i_date>=datetime.date(2004,4,1) and i_date<datetime.date(2015,1,1):
            modisdata.close()
            modisdata  = Dataset(modis_path+'UHBSCA'+str(iyear)+'%02d'%(imonth)+'_H.nc', format='NETCDF4')
            modisdataset = modisdata.variables['SCA']          
        temp_year = iyear
        temp_month = imonth  

    rssca1 = modisdataset[iday-1,:,:]
    rssca  = np.where(rssca1==200,1.,0.)
    rssca  = rssca*mask
    asca   = asca + rssca
    dsca.append(np.sum(rssca))
#################figure######################################
# plt.imshow(asca,origin='lower left')
# plt.title('asca')
# plt.colorbar()
# plt.savefig('D:\\asca.jpg')  #save the figure
# plt.show()
# plt.clf()
print 'The program is end!'
modisdata.close()