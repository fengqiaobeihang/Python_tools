'''
@author: hongyi
'''
import numpy as np
import datetime
from datetime import timedelta
import matplotlib.pylab as plt
from netCDF4 import Dataset
from netCDF4 import MFDataset
import pandas as pd

class stationinfor:
    def __init__(self):
        self.dempath = '../data/dem.txt'        
        self.filepath = '../data/stationInfor.dat'
        dem = np.flipud(np.genfromtxt(self.dempath,skiprows = 6)) 
        self.dem = dem
        stationInfor = np.genfromtxt(self.filepath)
        stationInfor =np.array(stationInfor)
#         print stationInfor
        self.stationsID = np.array(stationInfor[:,0], dtype='int')
        self.NO = len(open(self.filepath,'rU').readlines())
        self.station_lat = stationInfor[:,1]
        self.station_lon = stationInfor[:,2]
        self.station_alt = stationInfor[:,3]

    def readdata(self,istation):
        station_name = self.stationsID[istation]
        data_temp = np.genfromtxt('stationData/'+str(station_name) + '.txt',dtype='int')
        # np.savetxt('new.txt',data_temp)
        self.station_no = data_temp[:,0]
#         print station_no
        self.year = data_temp[:,1]
        self.month = data_temp[:,2]
        self.day = data_temp[:,3]
        snow_depth = data_temp[:,4]
        air_t_s = data_temp[:,5]
        windspeed_s = data_temp[:,6]
        rh_s = data_temp[:,7]
        press_s = data_temp[:,8]
        preci_s = data_temp[:,9]
        sw_d = data_temp[:,10]
        sw_u = data_temp[:,11]
        lw_d = data_temp[:,12]
        lw_u = data_temp[:,13]
        
        return snow_depth, air_t_s, windspeed_s, rh_s, press_s, preci_s, sw_d, sw_u, lw_d, lw_u


class forcing_input(object):
    '''
    classdocs
    '''
    def __init__(self):
        self.sbc = 5.6704e-8            #Stefan-Boltzman constant
    
    def down_longwave(self,air_t,ev,slope):
        # ev, vapor pressure (mb)
        # air_t,air temperature,(K)
        
        # The parameterization is from Jordan(1991, SNTHERM model)
        #Calculate the emissivity of air
        # First the Idso2 formula 
        slope=slope/180.*np.pi
        emiss_air = (0.70 + 5.95e-5 * ev * np.exp(1500./air_t))
        #Then the Wachtmann correction
        emiss_air = -0.792 + (3.161 * emiss_air) - (1.573 * emiss_air**2.)
        lw_down = 0.98*(self.sbc * emiss_air * air_t**4.)*(np.cos(slope/2.)**2.)
        
        return lw_down

    def down_shortwave(self,lon,lat,Nyear,Nday,Nhour,slope,aspect,st):
        #taut,total transmission coefficient
        #taub,direct shortwave radiation coefficient
        #taud,diffuse radiation transimission coefficient
        #E0,radius vector
        #Zenith
        #I0,solar constant
        #dang,Daily angle
        #hang,Hourly angle
        #tang,ture sun angle
        #decli,solar declination
        #selev,solar elevation
        #sazi,solar azimuth
        #lat,latitute
        #lon,longitude
        #slope,slope
        #aspect,aspect
        #s0,potential solar radiation on horizontal surface
        #ss,potential solar radiation on slope
        #st,total solar radiation(beam+diffuse)at Dadongshu station
        #sb,solar direct(beam) radiation
        #sd,solar diffuse radiation
        #kstar,sb+sd

        lon=lon/180.*np.pi
        lat=lat/180.*np.pi
        slope=slope/180.*np.pi
        aspect=(aspect)/180.*np.pi
        
        I0=1367.
        N0 = 79.6764+0.2422*(Nyear-1985) - np.int((Nyear-1985)/4.)

        dang=2.*np.pi*(Nday-N0)/365.2422
        #correct for ture sun angle
        tang=Nhour-4.*(120.-lon*180./np.pi)/60.+(0.000043+0.002061*np.cos(dang)-0.032040*np.sin(dang)\
                                      -0.014974*np.cos(2.*dang)-0.040685*np.sin(2.*dang))*229.2/60.
        hang=((tang-12.)*15./180.)*np.pi
        
        E0=1.000110+0.034221*np.cos(dang)+0.001280*np.sin(dang)+0.000719*np.cos(2.*dang) \
            +0.000077*np.sin(2.*dang)
        decli=0.006918-0.399912*np.cos(dang)+0.070257*np.sin(dang)-0.006758*np.cos(2*dang) \
            +0.000907*np.sin(2*dang) - 0.002697*np.cos(3*dang) + 0.00148*np.sin(3*dang)
        
        cosz = np.sin(decli)*np.sin(lat) + np.cos(decli)*np.cos(lat)*np.cos(hang)
#         sazi=np.arccos((np.sin(selev)*np.sin(lat)-np.sin(decli))/np.cos(selev)/np.cos(lat))
        
#         U=np.sin(lat)*np.cos(slope)-np.cos(lat)*np.cos(slope)*np.cos(aspect)
#         V=np.sin(lat)*np.sin(slope)*np.cos(aspect)+np.cos(lat)*np.cos(slope)
#         W=np.sin(slope)*np.sin(aspect)
        U = np.arcsin(np.sin(slope)*np.cos(aspect)*np.cos(lat) + np.cos(slope)*np.sin(lat))
        V = hang + np.arctan(np.sin(aspect)*np.sin(slope)/ \
                             (np.cos(slope)*np.cos(lat) - np.cos(aspect)*np.sin(slope)*np.sin(lat)))
#         
        s0 = I0*E0*cosz
#         ss = I0*E0*(U*np.sin(decli)+V*np.cos(decli)*np.cos(hang)+W*np.cos(decli)*np.sin(hang))
        ss = I0*E0*(np.sin(U)*np.sin(decli) + np.cos(U)*np.cos(decli)*np.cos(V)) #(Dewalle,2008)
#         plt.imshow(ss);plt.colorbar();plt.show();plt.clf();plt.close()
#         print ss
        ss = np.where(stationinfor().dem<0., 0., ss)
        #Compute sloar input in Grid

        taut = np.where(s0 <=0., 1., st/s0)
        taud = taut*(1.-np.exp(0.6*(1-0.76/taut)/(0.76-0.4)))
        taud = np.where(taud>1.,1.,taud)
        taud = np.where(taud<0.,0.,taud)

#           dirgrid = ss*(1-taud)

        difgrid = s0*taud*((np.cos(0.5*slope))**2)+(1.-((np.cos(0.5*slope))**2))*0.3*s0*taut    # self-shade effect,(Dewalle,2008)
        dirgrid = ss*taut - ss*taud
        dirgrid = np.where(s0 <= 0., 0., dirgrid)
        dirgrid = np.where(dirgrid <= 0., 0., dirgrid)
        kstar = dirgrid/np.cos(slope) + difgrid  #convert the result to horizontal surface
        kstar = np.where(kstar > I0*E0, I0*E0, kstar)
        kstar = np.where(stationinfor().dem<0., 0., kstar)
        kstar = np.where(kstar <= 0., 0., kstar)

        return kstar

        
    def RH(self, temp_base, rh_base, press, airtemp):
#     def RH(self, temp_base, rh_base, press, airtemp):
#         PRO CALC_rh_wind, datum, temp_base, rh_base, press, temperature, $
#         relhum, ulr, wind, wind_base
        #Extrapolate Relative Humidity
        #airtemp -- air temperature (K)
#         a = 7.5
#         if temp_base <= 273.15:
#             a = 9.5
#         b = 237.3
#         if temp_base <= 273.15:
#             b = 265.5
# 
#         tempC = temp_base - 273.15
# 
#         es_base = 6.11 * 10.**((a * tempC)/(tempC + b)) #mb
#         ev = (rh_base/100 * es_base)
# 
#         q = ((0.622 * ev)/(press - (0.378 * ev)))
#         
#         rh = (q * press/0.622)/(6.11 * 10.**((a * (airtemp - 273.15)) \
#         /((airtemp - 273.15) + b))) * 100.
#         
        a_base = 7.5
        if temp_base <= 273.15:
            a_base = 9.5
        b_base = 237.3
        if temp_base <= 273.15:
            b_base = 265.5

        tempC = temp_base - 273.15

        es_base = 6.11 * 10.**((a_base * tempC)/(tempC + b_base)) #mb
        ev_base = (rh_base/100. * es_base)

        q = ((0.622 * ev_base)/(press - (0.378 * ev_base)))
        
        rh = (q * press/0.622)/(6.11 * 10.**((a_base * (airtemp - 273.15)) \
        /((airtemp - 273.15) + b_base))) * 100.
        
#         a = 0.*rh.copy()
#         b = 0.*rh.copy()
        a = np.where(airtemp<=273.15,9.5, 7.5)
        b = np.where(airtemp<=273.15,265.5, 237.3)
        es = 6.11 * 10.**((a * (airtemp - 273.15))/((airtemp - 273.15) + b))
        ev = (rh/100. * es)

        rh = np.where(stationinfor().dem<0., 0., rh)
        ev = np.where(stationinfor().dem<0., 0., ev)

        return rh, ev

    def wind(self, wind_base1, wind_base2):
        datum = stationinfor().dem - stationinfor().station_alt[0]
        windlapse = datum * (wind_base2-wind_base1)/(stationinfor().station_alt[1] - stationinfor().station_alt[0])
        windspd = windlapse + wind_base1
        windspd = np.where(windspd < wind_base2, wind_base2, windspd)
        windspd = np.where(stationinfor().dem<0., wind_base2, windspd)
        return windspd

    def airt(self,airt_s):
        #airt_s -- air temperature at station (K)
        airtemp= -0.65*(stationinfor().dem - stationinfor().station_alt[0])/100. + airt_s
        airtemp = np.where(stationinfor().dem<0., 273.15, airtemp)   
        return airtemp 

    def press(self,press_s):
        #mb   
        pressure = (stationinfor().dem -stationinfor().station_alt[0])/9. + press_s
        pressure = np.where(stationinfor().dem<0., 0., pressure)   
        return pressure
