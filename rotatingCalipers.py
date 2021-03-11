import copy
import math
import numpy as np
import numpy.linalg as la
class _point:
  
  def __init__(self,  coor=[0,0], index=0):
    self.index = index
    self. coor = coor
    
  def shift(self, originalPoint):
    newX = self.coor[0]- originalPoint.coor[0]
    newY = self.coor[1] - originalPoint.coor[1]
    
    #print ("after Shift: "+str([newX, newY]))
    
    return _point([newX,newY], self.index)
    
  def counterClockRotate(self, slop):
    c,s = np.cos(slop), np.sin(slop)
    j = np.matrix([[c,-s],[s,c]])
    m = np.dot(j,self.coor)
    
    #print(j)
    #print(m)
    
    #print ("after rotated: "+str(m.T[0].tolist()[0]+m.T[1].tolist()[0]))
    
    return _point(m.T[0].tolist()[0]+m.T[1].tolist()[0], self.index)
    
  def getPointAngle(self):
    cosang = np.dot(self.coor, [1,0])
    sinang = la.norm(np.cross(self.coor, [1,0]))
    return np.arctan2(sinang, cosang)
    
  def printPt(self):
  
    print("Point index: "+str(self.index)+" : "+str(self.coor))
    
class _line:

  def __init__(self, slop=0, point=_point()):
    self.slop = slop
    self.passPoint = point
    self.abcForm = [0,0,0]
    self.setAbcForm()
    
  def setAbcForm(self):
    self.abcForm[0] = math.tan(self.slop)
    self.abcForm[1] = -1
    self.abcForm[2] = self.passPoint.coor[1] - self.abcForm[0]*self.passPoint.coor[0]
    #print(self.abcForm)
    
  def distanceToLine(self, line):
    '''
    return abs(self.abcForm[0]*line.passPoint.coor[0] + self.abcForm[1]*line.passPoint.coor[1] + self.abcForm[2])\
            /math.sqrt(self.abcForm[0]*self.abcForm[0]+1)
    '''
    return self.distanceToPoint(line.passPoint.coor)
    
  def distanceToPoint(self, pt):
    return abs(self.abcForm[0]*pt[0] + self.abcForm[1]*pt[1] + self.abcForm[2])\
            /math.sqrt(self.abcForm[0]*self.abcForm[0]+1)
            
  def getAngleToPoint(self, pt):
    
    shiftedPt = pt.shift(self.passPoint)
    #print("shiftedPt: ")
    #shiftedPt.printPt()
    rotatedPt = shiftedPt.counterClockRotate(-self.slop)
    #rotatedPt.printPt()
    return rotatedPt.getPointAngle()
    
  def printLine(self):
    print("Line information: (1)slop: "+str(self.slop)+"  (2): Passing points:")
    self.passPoint.printPt()
    print("(3) abc form: y="+str(self.abcForm[0])+"x+ "+str(self.abcForm[2]))
    
  def crossPoint(self,otherLine):
    x = 0
    y = 0
    if abs(self.slop - otherLine.slop)>0.001:  
      dx = self.abcForm[0]- otherLine.abcForm[0]
      dc = otherLine.abcForm[2]-self.abcForm[2]
      x = dc/dx
      y = self.abcForm[0]* x + self.abcForm[2]
    return x,y         
  

class _rectangle:

  
  def __init__(self, supportPoints=[_point(),_point(),_point(),_point()], angle=0, supportIndex = 0):
    self.supportPoints= supportPoints
    self.angle = angle
    self.supportIndex = supportIndex
    self.lines = [_line(self.angle, self.supportPoints[0]), _line(self.angle + math.pi/2, self.supportPoints[1])\
                 ,_line(self.angle + math.pi, self.supportPoints[2]), _line(self.angle + 3*math.pi/2, self.supportPoints[3])]
    self.area = self.calculateArea()
    
    
  def calculateArea(self):
    bottomLine = self.lines[self.supportIndex]
    rightLine = self.lines[(self.supportIndex+1)%4]
    topLine = self.lines[(self.supportIndex+2)%4]
    leftLine = self.lines[(self.supportIndex+3)%4]
    
    #bottomLine.printLine()
    #topLine.printLine()
    #print("bt: "+ str(bottomLine.distanceToLine(topLine)))
    #rightLine.printLine()
    #leftLine.printLine()
    #print("rl: "+ str(rightLine.distanceToLine(leftLine)))
    bx,by = bottomLine.crossPoint(rightLine)
    rx,ry = rightLine.crossPoint(topLine)
    tx,ty = topLine.crossPoint(leftLine)
    lx,ly = leftLine.crossPoint(bottomLine)
    self.x =[bx,rx,tx,lx,bx]
    self.y =[by,ry,ty,ly,by]
    return bottomLine.distanceToLine(topLine)*rightLine.distanceToLine(leftLine)
    
  def printRect(self):
    print("Rectangle supportPoints: ")
    print("support Index: " + str(self.supportIndex))
    for sp in self.supportPoints:
      print(str(sp.index)+ ": " + str(sp.coor))
    print("lines: ")
    for l in self.lines:
      print(str(l.passPoint.index)+ ": "+ str(l.slop))
    print("angle: "+ str(self.angle))
    print("Area: "+ str(self.area))
    



#import a set of points which constitute a convext hull in counter clock direction
#return a _rectangle which has the minimum area cover all the points
def findMinAreaRect(convexHullPts):
  iniRect = initialRectangle(convexHullPts)
  unvisitedEdges = initialEdgeDic(convexHullPts)
  #print("unvisitedEdges:")
  #print(unvisitedEdges)
  minArea = iniRect.area
  preRect = copy.deepcopy(iniRect)
  finalRect = copy.deepcopy(iniRect)
  #print("preRect:")
  #preRect.printRect()
  count=0
  while unvisitedEdges:
    count+=1
    #print(str(count)+" :"+str(unvisitedEdges))
    
    preRect, newEdges = rotateMinAngle(convexHullPts, preRect)
    #print("*************new Edges**************")
    #print(newEdges)
    #print("preRect:")
    #preRect.printRect()
    updateUnvisitedEdges(unvisitedEdges,newEdges)
    
    if minArea > preRect.area:
      finalRect = copy.deepcopy(preRect)
      minArea = preRect.area
      
    if count >1000:
      print("infinite loop")
      break
      
  return finalRect
      


def initialRectangle(convexHullPts):
  count = 0
  bottomPointIndex = 0
  topPointIndex = 0
  leftPointIndex = 0
  rightPointIndex = 0
  minX = 10000
  minY = 10000
  maxX = -10000
  maxY = -10000
  
  for pt in convexHullPts:
    if pt[0] > maxX:
      maxX = pt[0]
      rightPointIndex = count
    if pt[0] < minX:
      minX = pt[0]
      leftPointIndex = count
    if pt[1] > maxY:
      maxY = pt[1]
      topPointIndex = count
    if pt[1] < minY:
      minY = pt[1]
      bottomPointIndex = count
      
    count+=1
  #print(bottomPointIndex)
  #print(topPointIndex)
  #print(leftPointIndex)
  #print(rightPointIndex)
  bottomPoint = _point(convexHullPts[bottomPointIndex], bottomPointIndex)
  topPoint = _point(convexHullPts[topPointIndex], topPointIndex)
  leftPoint = _point(convexHullPts[leftPointIndex], leftPointIndex)
  rightPoint = _point(convexHullPts[rightPointIndex], rightPointIndex)
  
  return _rectangle([bottomPoint,rightPoint,topPoint ,leftPoint] , 0)

  
def initialEdgeDic(convexHullPts):
  edgeDic = set()
  for i in range(len(convexHullPts)):
    edgeDic.add((i,(i+1)%len(convexHullPts)))
  return edgeDic

def updateUnvisitedEdges(unvisitedEdges, newEdges):
  for newEdge in newEdges:
    if newEdge in unvisitedEdges:
      unvisitedEdges.remove(newEdge)

def rotateMinAngle(convexHullPts, preRect):
  rotates =[]
  totalPoints = len(convexHullPts)
  indexList=[]
  for i in range(4):
    indexList.append(preRect.supportPoints[i].index)
    nextIndex = (preRect.supportPoints[i].index+1)%totalPoints
    #print("nextIndex: "+str(nextIndex))
    rotates.append(preRect.lines[i].getAngleToPoint(_point(convexHullPts[nextIndex],nextIndex)))
  
  miniIndex = rotates.index(min(rotates))
  #print(rotates)
  #print ('min rotates:' + str(min(rotates)))
  #print("miniIndex: "+ str(miniIndex))
  #print(indexList)
  preRect.angle = preRect.angle + min(rotates)
  newEdges =[]
  
  for i in range(4):
    if abs(rotates[i] - min(rotates)) < 0.000001:
      oldIndex = indexList[i]
      newPointIndex = (oldIndex+1)%len(convexHullPts)
      newPoint = _point(convexHullPts[newPointIndex], newPointIndex)
      newEdges.append((oldIndex,newPointIndex))
      preRect.supportPoints[i] = newPoint
      preRect.supportIndex = i
      preRect.lines[i] = _line(preRect.lines[i].slop+min(rotates), _point(convexHullPts[newPointIndex],newPointIndex))
    else:
      preRect.lines[i].slop += min(rotates)
      preRect.lines[i].setAbcForm()
  
  preRect.area = preRect.calculateArea()
  
  return preRect, newEdges
 
  
  
  
