# -*- coding: utf-8 -*-
#画多条折线图
import numpy             as     np
import matplotlib.pyplot as     plt
import xlrd
from   matplotlib.ticker import MultipleLocator, FormatStrFormatter
from   matplotlib.pylab  import *

data  = xlrd.open_workbook('C:/Users/sdh/Desktop/RH2.xlsx')#open file
sheet = data.sheets()[0] #Select table (which tables are shown in brackets)
# Generate fake data
x     = sheet.col_values(0)#[0:365]#表示取第0-第365行
y     = sheet.col_values(1)#[0:365] #col = sheet.col_values(0)##获取第一列的数据
fig   = plt.figure(figsize=(15,5))#设置图形界面的大小
ax    = fig.add_subplot(111)
#画折线图
ax.plot(x,y,label='RH2')
###############################
#修改X轴的字体大小
for xtick in ax.xaxis.get_major_ticks():
	xtick.label1.set_fontsize(16)
#修改Y轴的字体大小
for ytick in ax.yaxis.get_major_ticks():
	ytick.label1.set_fontsize(16)
# plt.legend(loc='upper right')#图例
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
ax.set_xlim(0,14800)#X轴的数值范围
# ax.set_ylim(0,100)#Y轴的数值范围
# #设置x轴的数据
# ax.set_xticks([15, 45, 75, 115, 145,175,205,235,265,295,325,355])
# ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']) 
# #设置y轴的数据
ax.set_yticks([0,10,20,30,40,50,60,70,80,90,100])
# ax.set_yticklabels(['0','10','20','30','40','50','60','70','80','90','100'])  
#字体大小及类型设置
font = matplotlib.font_manager.FontProperties(family='times new roman',size=20)  
ax.set_xlabel('Date',  fontproperties=font)#X轴标签
ax.set_ylabel('RH2(%)',fontproperties=font)#Y轴标签
ax.set_title('2-m relative humidity',fontproperties=font)#图形标题
plt.show()

