import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

yakoutemp = np.loadtxt('outdata.txt')
yakoutempv2 = np.loadtxt('outdatav9.txt')
yakoutempv0 = np.loadtxt('outdatayk.txt')
plt.plot(yakoutemp)
plt.show()
meteodata_index=[
"TIMESTAMP",
"RainSnowacc",
"RainSnow",
"Snow_Depth",
"Soil_T_0_Avg",
"Soil_T_4_Avg",
"Soil_T_10_Avg",
"Soil_T_20_Avg",
"Soil_T_40_Avg",
"Soil_T_80_Avg",
"Soil_T_120_Avg",
"Soil_T_160_Avg",
"Soil_M_4_Avg",
"Soil_M_10_Avg",
"Soil_M_20_Avg",
"Soil_M_40_Avg",
"Soil_M_80_Avg",
"Soil_M_120_Avg",
"Soil_M_160_Avg",
"shf_Avg_1",
"shf_Avg_2",
"PTemp_Avg",
"short_up_Avg",
"short_dn_Avg",
"long_up_corr_Avg",
"long_dn_corr_Avg",
"Rn_Avg",
"albedo_Avg",
"cnr4_T_C_Avg",
"TargTemp_Avg",
"WindDirect_deg",
"WindSpeed_Avg",
"AirTemp_C_Avg",
"RH_Avg",
"Barometer_KPa_Avg",
"Rain_mm_Avg",
"PTemp",
"batt_volt"
]
meteodata = pd.read_csv('meteo2014.csv', index_col=0,names =meteodata_index,parse_dates=True,skiprows=3)
# read seperate data
daily_meteo = meteodata.resample('H')
batt = daily_meteo['AirTemp_C_Avg']
plt.plot(batt,label='obs')
plt.plot(yakoutemp-273.14,label='wrf')
# plt.plot(yakoutempv2-273.14,label='yk')
# plt.plot(yakoutemp[0:8000]-273.14,batt[0:8000],'o')
# plt.plot([-40,20],[-40,20])
# plt.plot([-40,0],[0,0])
# plt.plot([0,0],[-40,0])

# plt.legend()
# batt.plot()
plt.show()
print 'end'