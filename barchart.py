# -*- coding: utf-8 -*-
#画多条折线图
import numpy             as     np
import matplotlib.pyplot as     plt
import xlrd
from   matplotlib.ticker import MultipleLocator, FormatStrFormatter
from   matplotlib.pylab  import *

data  = xlrd.open_workbook('C:/Users/sdh/Desktop/FSC.xlsx')#open file
sheet = data.sheets()[0] #Select table (which tables are shown in brackets)
# Generate fake data
x     = sheet.col_values(0)#[0:365]#表示取第0-第365行
y2001 = sheet.col_values(1)#[0:365] #col = sheet.col_values(0)##获取第一列的数据
y2002 = sheet.col_values(2)#[0:365] #row = sheet.row_values(0)##获取第一行的数据
y2003 = sheet.col_values(3)#[0:365]
y2004 = sheet.col_values(4)#[0:365]
y2005 = sheet.col_values(5)#[0:365]
y2006 = sheet.col_values(6)#[0:365]
y2007 = sheet.col_values(7)#[0:365]
y2008 = sheet.col_values(8)#[0:365]
y2009 = sheet.col_values(9)#[0:365]
y2010 = sheet.col_values(10)#[0:365]
y2011 = sheet.col_values(11)#[0:365]#
# print 'y2001',y2001
fig   = plt.figure(figsize=(20,10))#设置图形界面的大小
ax    = fig.add_subplot(111)
#画折线图
ax.plot(x,y2001,label='2001')
ax.plot(x,y2002,label='2002')
ax.plot(x,y2003,label='2003')
ax.plot(x,y2004,label='2004')
ax.plot(x,y2005,label='2005')
ax.plot(x,y2006,label='2006')
ax.plot(x,y2007,label='2007')
ax.plot(x,y2008,label='2008')
ax.plot(x,y2009,label='2009')
ax.plot(x,y2010,label='2010')
ax.plot(x,y2011,label='2011')
###############################
#修改X轴的字体大小
for xtick in ax.xaxis.get_major_ticks():
	xtick.label1.set_fontsize(16)
#修改Y轴的字体大小
for ytick in ax.yaxis.get_major_ticks():
	ytick.label1.set_fontsize(16)
plt.legend(loc='upper right')#图例
plt.grid(True,axis='y')#设置网格线
# ax.xaxis.grid(True, which='major') #x坐标轴的网格使用主刻度
# ax.yaxis.grid(True, which='minor') #y坐标轴的网格使用次刻度
#设置背景网格线的颜色，样式，尺寸和透明度
# plt.grid(color='#95a5a6',linestyle='--', linewidth=1,axis='y',alpha=0.4)
#显示次刻度标签的位置,没有标签文本
# ax.xaxis.set_minor_locator(xminorLocator)
# ax.yaxis.set_minor_locator(yminorLocator)
# ax.xaxis.grid(True, which='major') #x坐标轴的网格使用主刻度
# ax.yaxis.grid(True, which='minor') #y坐标轴的网格使用次刻度
ax.set_xlim(0,365)#X轴的数值范围
ax.set_ylim(0,100)#Y轴的数值范围
#设置x轴的数据
ax.set_xticks([15, 45, 75, 115, 145,175,205,235,265,295,325,355])
ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']) 
#设置y轴的数据
ax.set_yticks([0,10,20,30,40,50,60,70,80,90,100])
ax.set_yticklabels(['0','10','20','30','40','50','60','70','80','90','100'])  
#字体大小及类型设置
font = matplotlib.font_manager.FontProperties(family='times new roman',size=20)  
ax.set_xlabel('Date',  fontproperties=font)#X轴标签
ax.set_ylabel('FSC(%)',fontproperties=font)#Y轴标签
ax.set_title('The annual change in the FSC',fontproperties=font)#图形标题
plt.show()
# N     = 5
# menMeans = x #(20, 35, 30, 35, 27)
# menStd   = (2, 3, 4, 1, 2)

# ind   = np.arange(N)  # the x locations for the groups
# width = 0.35       # the width of the bars
# fig, ax = plt.subplots()
# rects1  = ax.bar(ind, menMeans, width, color='r', yerr=menStd)

# womenMeans = y #(25, 32, 34, 20, 25)
# womenStd   = (3, 5, 2, 3, 3)
# rects2     = ax.bar(ind + width, womenMeans, width, color='y', yerr=womenStd)

# # add some text for labels, title and axes ticks
# ax.set_xlabel('Time')
# ax.set_ylabel('Windspeed')
# ax.set_title('Changes of Wind Velocity Figure')
# ax.set_xticks(ind + width)
# ax.set_xticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))
