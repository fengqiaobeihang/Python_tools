# -*- coding: utf-8 -*-
#shaodonghang
#2017-1-3
import locale
import codecs
import sys
import re
import math

import os
import os.path
import sys
from os.path import walk
import linecache

type=sys.getfilesystemencoding()

folderpath1="D:/abc/data"
# folderpath2="C:/Users/dy/Desktop/data/data2/"
folderpath3="D:/abc/rst"
for dirpath,filename,filenames in os.walk(folderpath1):
    for filename in filenames:
        if os.path.splitext(filename)[1]=='.TXT': #获取指定格式的文件,后缀名一定要对应大小写
            DEMfilepath   = os.path.join(dirpath,filename)
            # FSFilename    = filename.replace('u','v')
            RSFilename    = filename.replace('SURF_CLI_CHN_MUL_DAY-WIN-11002','wind-speed')
            # FSfilepath    = os.path.join(folderpath2,FSFilename)#合并绝对路径
            RESULTfilepath= os.path.join(folderpath3,RSFilename)

            savefile=codecs.open(RESULTfilepath,"w","utf-8")
            
            fDEMPointer= open(DEMfilepath,"r")
            # fFSPointer = open(FSfilepath,"r")
            nocls=linecache.getlines(DEMfilepath)
            rows =len(nocls)#文件的总行数
            print 'rows=',rows
            # print 'nocls=',nocls
            data=fDEMPointer.readlines()
            for line in data:
                v=line.strip(' ').split()
                if v[0]=='52643':
                    savefile.writelines(v[7]+"\n")
            aa=['dfadf','dfasf\n','fsaf\n','dfasdf']
            savefile.writelines(aa)
            savefile.close()


            # while True:
            #     ResultLine=''
            #     DEMLine=(fDEMPointer.readline()).decode('utf-8').encode(type)
            #     # FSLine=(fFSPointer.readline()).decode('utf-8').encode(type)
            #     if DEMLine:
            #         # print DEMLine[39:41]
            #         DEMLine=DEMLine.strip(' ')#删除文件中的指定字符
            #         # FSLine=FSLine.strip()
            #         DEMValues=DEMLine.split()#通过指定分隔符对字符串进行切片,默认为以空格切片,空格数无论多少都可以
            #         # FSValues=FSLine.split(' ')
            #         # print DEMValues[0]
            #         # print DEMValues[7]
            #         # print DEMValues
            #         lines=len(DEMValues)
            #         # print 'lines=',lines
            #     # for i in range(len(DEMValues)):
            #     #     print i
            #     # try:
            #         # tempPlus=[]

            #         if DEMValues[0]=='52643':
            #             # print DEMValues[i]
            #             # DEMValue=float(DEMValues[m])
            #             # FSValue=float(FSValues[m])
            #             tempPlus=DEMValues[7]
                        
            #             # tempws  =tempPlus/10
            #             # Result  =tempws
            #             #ResultLine.append(tempPlus)
            #             print(tempPlus)
            #             savefile.write(tempPlus+"\n")
            #             # ResultLine=ResultLine+str(tempPlus)+" "
            #         else:
            #             # ResultLine=ResultLine+str(tempPlus)
            #             pass
            #     # except:
            #     #     pass
            #     else:
            #         break
                    

            # fDEMPointer.close()
            # # fFSPointer.close()
            # savefile.close()
fp=open("D:/abc/rst/wind-speed-200001.TXT")
print(len(fp.readlines()))
print "done"