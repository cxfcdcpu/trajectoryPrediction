import os
count=1
maxp = 3837
for fn in range(maxp):
  count+=1
  if count%10000 ==0:
    print(count)
  fname = "./splitedData/"+str(fn)+".txt"
  f = open(fname,"r")
  if len(f.read()) == 0 :
    f.close()
    os.remove(fname)
  f.close()

