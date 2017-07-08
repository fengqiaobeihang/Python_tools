import os,shutil
path = "C:/Users/sdh/Desktop/txtfile"
f=open(path + "/result.txt","a")
for r,d,fi in os.walk(path):
     for files in fi:
         if files.endswith(".txt"):                         
              g=open(os.path.join(r,files))
              shutil.copyfileobj(g,f)
              g.close()
f.close()