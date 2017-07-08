import numpy               as np
import xlrd
import matplotlib
import matplotlib.pyplot   as plt
from   matplotlib.colors   import *
from   scipy.stats         import gaussian_kde
data  =xlrd.open_workbook('C:\\Users\\lhy\\Desktop\\1.xlsx')#open file
sheet =data.sheets()[0] #Select table (which tables are shown in brackets)
# Generate fake data
x     = sheet.col_values(1) #col = sheet.col_values(0)##Get columns:(which columns of data are shown in brackets) 
y     = sheet.col_values(2) #row = sheet.row_values(0)##Line access:(which the data are shown in brackets) 
# Calculate the point density
xy    = np.vstack([x,y])
z     = gaussian_kde(xy)(xy)
cm    = plt.cm.get_cmap('PuRd')
fig   =plt.scatter(x, y, c=z, s=20, edgecolor='',cmap=cm)
plt.colorbar(fig)   
ax    =plt.gca()
ax.set_xlabel('Time')
ax.set_ylabel('windspeed')
plt.title('')
plt.show()
