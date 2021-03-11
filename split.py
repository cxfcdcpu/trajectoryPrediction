import numpy as np
import pandas as pd
import geopandas as gpd
import json
from shapely.geometry import LineString, Point
from multiprocessing import  Pool


    
    
df = pd.read_csv('./goodData.csv')
vl = json.loads('[' + df.POLYLINE.str.cat(sep=',') + ']')
maxp = 3837
ft = []
for fn in range(maxp):
  f = open("./splitedData/"+str(fn)+".txt","a")
  ft.append(f)


count=1
for v in vl:
  count+=1
  if count%10000 ==0:
    print(count)
  pnum = len(v)
  ft[pnum].write(str(v))
  ft[pnum].write('\n')
  
for f in ft:
  f.close()
  
print('all done')
