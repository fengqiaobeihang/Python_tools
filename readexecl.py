import numpy as np
import xlrd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import *
import xlrd
# file='C:\Users\lhy\Desktop\1.xlsx'
wb=xlrd.open_workbook('C:\\Users\\lhy\\Desktop\\1.xlsx')
ws=wb.sheet_by_name('Sheet1')
dataset=[]
for r in xrange(ws.nrows):
	col=[]
	for c in range(ws.ncols):
		col.append(ws.cell(r,c).value)
	dataset.append(col)

from pprint import pprint
pprint(dataset)