from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import xlrd

data   =xlrd.open_workbook('C:\\Users\\lhy\\Desktop\\1.xlsx')#open file
sheet  =data.sheets()[0] #Select table (which tables are shown in brackets)
x      = sheet.col_values(1) #col = sheet.col_values(0)##Get columns:(which columns of data are shown in brackets) 
y      = sheet.col_values(2) #row = sheet.row_values(0)##Line access:(which the data are shown in brackets) 
fig    = plt.figure()
ax     = fig.add_subplot(111, projection='3d')
for c, z in zip(['r', 'g', 'b', 'y'], [30, 20, 10, 0]):
    xs = np.arange(20)
    ys = np.random.rand(20)

    # You can provide either a single color or an array. To demonstrate this,
    # the first bar of each set will be colored cyan.
    cs = [c] * len(xs)
    cs[0] = 'c'
    ax.bar(xs, ys, zs=z, zdir='y', color=cs, alpha=0.8)

ax.set_xlabel('Time')
ax.set_ylabel('Windspeed')
ax.set_zlabel('Value')

plt.show()