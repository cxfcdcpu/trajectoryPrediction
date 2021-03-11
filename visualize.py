import numpy as np
import pandas as pd
import geopandas as gpd
import json
from shapely.geometry import LineString, Point
from multiprocessing import  Pool
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import chan as ch
import rotatingCalipers as rc
import math
sns.set(style="darkgrid")


def loadData(n):
  fname = './splitedData/'+str(n)+'.txt'
  with open(fname,'r') as file:
    datas = file.readlines()
    return datas
    
def createLineList(m,datas):
  ts = datas[m-1].replace('[','').replace(']','').strip()
  trajectory = ts.split(',')
  return trajectory

def plotHull(trajectory):
  plts = []
  x=[]
  y=[]
  count=0
  for coor in trajectory:
    count+=1
    if count%2 == 1:
      x.append(float(coor))
    else:
      y.append(float(coor)) 
  count=0
  for v in x:
    plts.append([v,y[count]])
    count+=1
  convexHullPts = ch.convex_hull(plts)
  
  finalRect = rc.findMinAreaRect(convexHullPts)
  print(getRectCenter(finalRect))
  print(getLengthWidth(finalRect))
  print(finalRect.angle)
  plotRect(finalRect)
  
def plotRect(finalRect):
  plt.plot(finalRect.x,finalRect.y)

def getRectCenter(finalRect):
  cx=0
  cy = 0
  for i in range(4):
    cx+=finalRect.x[i]
    cy+=finalRect.y[i]
  return [cx/4, cy/4]
  
def getLengthWidth(finalRect):
  rl = 0
  tb = 0
  bp = (finalRect.x[0], finalRect.y[0])
  rp = (finalRect.x[1], finalRect.y[1])
  tb = euDis(bp,rp)
  
  lp = (finalRect.x[3], finalRect.y[3])
  rl = euDis(bp, lp)
  if rl > tb:
    return [rl, tb , finalRect.angle]
  else:
    return [tb, rl, finalRect.angle-math.pi/2]  
  
  
def euDis(p1, p2):
  return math.sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]))
  
def plotTr(trajectory):
  x,y = getTr(trajectory)
  plt.plot(x,y)
  plt.show()  
  
def getTr(trajectory):
  x=[]
  y=[]
  count=0
  for coor in trajectory:
    count+=1
    if count%2 == 1:
      x.append(float(coor))
    else:
      y.append(float(coor)) 
  return x,y  

def getMiniRect(trajectory):
  plts = []
  x=[]
  y=[]
  count=0
  for coor in trajectory:
    count+=1
    if count%2 == 1:
      x.append(float(coor))
    else:
      y.append(float(coor)) 
  count=0
  for v in x:
    plts.append([v,y[count]])
    count+=1
  convexHullPts = ch.convex_hull(plts)
  
  return rc.findMinAreaRect(convexHullPts)
  



def test1():
  try:
      mode=int(input('which data you want to import:'))
  except ValueError:
      print('Not a number')    
  datas = loadData(mode)

  print('There are total '+str(len(datas))+' trajectory')
  try:
      mode=int(input('Which trajectory you want to plot?'))
  except ValueError:
      print('Not a number')
  trajectory = createLineList(mode,datas)
  plotHull(trajectory)
  plotTr(trajectory)
  
def test2():
  try:
      mode=int(input('which data you want to import:'))
  except ValueError:
      print('Not a number')    
  datas = loadData(mode)
  maxRect = [0,0]
  for mode in range(100):  
    #print(mode) 
    trajectory = createLineList(mode,datas)
    try:
      finalRect = getMiniRect(trajectory)
    except:
      print(mode)
    else:  
      curRect = getLengthWidth(finalRect)
      maxRect[0] = max(maxRect[0],curRect[0])
      maxRect[1] = max(maxRect[1], curRect[1])
      print([curRect[0]*20000,curRect[1]*20000])
      
  print( maxRect)
  print([maxRect[0]*20000,maxRect[1]*20000])

def test3():
  try:
      mode=int(input('which data you want to import:'))
  except ValueError:
      print('Not a number')    
  datas = loadData(mode)

  print('There are total '+str(len(datas))+' trajectory')
  try:
      mode1=int(input('Which trajectory you want to plot?'))
  except ValueError:
      print('Not a number')
  trajectory = createLineList(mode1,datas)
  try:
    finalRect = getMiniRect(trajectory)
  except:
    print(mode)
  else:  
    add = plotRightRect(finalRect)
    x,y = plotRightPlot(trajectory, finalRect, add)
    stroke = getStroke(x,y)
    writeStroke(mode,mode1,stroke)
    plt.show()

def test4():
  for mode in range(20,190):   
    datas = loadData(mode)
    
    for mode1 in range(20,150):
      print("processing stroke "+ str(mode)+" traj "+str(mode1))
      trajectory = createLineList(mode1,datas)
      try:
        finalRect = getMiniRect(trajectory)
      except:
        print(mode)
      else:  
        add = findRightRect(finalRect)
        x,y = findRightPlot(trajectory, finalRect, add)
        stroke = getStroke(x,y)
        writeStroke(mode,mode1,stroke)
           
    
def plotRightRect(finalRect):
  curRect = getLengthWidth(finalRect)
  curCenter = getRectCenter(finalRect)
  rightX, rightY, add = getRight(finalRect, curRect,curCenter)
  plt.plot(rightX,rightY)
  return add
  
def findRightRect(finalRect):
  curRect = getLengthWidth(finalRect)
  curCenter = getRectCenter(finalRect)
  rightX, rightY, add = getRight(finalRect, curRect,curCenter)
  
  return add  

def getRight(finalRect,curRect,curCenter):
  shiftedPt = []
  rotatedPt = []
  for i in range(5):
    shiftedPt.append(shiftPt([finalRect.x[i],finalRect.y[i]],curCenter))
  maxX = 0
  maxY = 0
  for i in range(5):
    rotatedPt.append(rotatePt(shiftedPt[i], -curRect[2]))
    maxX = max(maxX, rotatedPt[i][0])
    maxY = max(maxY, rotatedPt[i][1])
  add = 0
  if maxX < maxY:
    add = math.pi/2
    for i in range(5):
      rotatedPt[i]= rotatePt(rotatedPt[i], add)
  rx, ry =getPlotReadyXY(rotatedPt)
  return rx,ry,add
  
  
def shiftPt(pt, curCenter):
  return [pt[0]-curCenter[0],pt[1]-curCenter[1]]
  
def rotatePt(pt, slop):
  c,s = np.cos(slop), np.sin(slop)
  j = np.matrix([[c,-s],[s,c]])
  m = np.dot(j,pt)
  return m.T[0].tolist()[0]+m.T[1].tolist()[0]

def getPlotReadyXY(pts):
  x = []
  y = []
  for pt in pts:
    x.append(20000*pt[0]+1500)
    y.append(20000*pt[1]+600)
  return x, y

def plotRightPlot(trajectory,finalRect, add):
  x,y = getTr(trajectory)
  curRect = getLengthWidth(finalRect)
  curCenter = getRectCenter(finalRect) 
  
  print(curRect)
  print(finalRect.angle)
  rightX,rightY = getRightTr(x,y, curRect,curCenter,add)
  
  
  plt.plot(rightX,rightY)
  return rightX,rightY


def findRightPlot(trajectory,finalRect, add):
  x,y = getTr(trajectory)
  curRect = getLengthWidth(finalRect)
  curCenter = getRectCenter(finalRect) 
  
  print(curRect)
  print(finalRect.angle)
  rightX,rightY = getRightTr(x,y, curRect,curCenter,add)

  return rightX,rightY

  
def getRightTr(x,y,curRect,curCenter, add):
  shiftedPt = []
  rotatedPt = []
  for i in range(len(x)):
    shiftedPt.append(shiftPt([x[i],y[i]],curCenter))
  for i in range(len(shiftedPt)):
    rotatedPt.append(rotatePt(shiftedPt[i], -curRect[2]+add))
  return getPlotReadyXY(rotatedPt)  


def getStroke(x,y):
  res = []
  for i in range(len(x)):
    res.append([x[i],y[i]])
  return res
  
def writeStroke(sID, iID, stroke):
  strokeFileName = "./strokes/"+str(sID)+"_"+str(iID)
  f = open(strokeFileName,"w")
  for p in stroke:
    f.write(str(round(p[0]))+" "+ str(round(p[1]))+"\n")
  f.close()
  

  
try:
    mode=int(input('Which test you want to execute?'))
except ValueError:
    print('Not a number')
else:
  if mode ==1:
    test1()
  if mode ==2:
    test2()
  if mode == 3:
    test3()
  if mode == 4:
    test4()
