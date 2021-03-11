import os
count=-1
maxp = 3837

def loadData(n):
  fname = './splitedData/'+str(n)+'.txt'
  if not os.path.exists(fname):
    return []
  with open(fname,'r') as file:
    datas = file.readlines()
    return datas
    
def mvToDump(n):
  fname = './splitedData/'+str(n)+'.txt'
  nname = './rareCases/'+str(n)+'.txt'
  if not os.path.exists(fname):
    return
  os.rename(fname,nname)

for fn in range(maxp):
  count+=1
  datas = loadData(count)
  if datas:
    print(str(count)+" : "+str(len(datas)))
  if datas and len(datas)>0 and len(datas) < 140:
    mvToDump(count)
  
  

