#ifndef USERCLASS_H
#define USERCLASS_H

#include<iostream>	//cout
#include<stdio.h>	//printf
#include<string.h>	//strlen
#include<string>	//string
#include<cstring>
#include<sys/socket.h>	//socket
#include<arpa/inet.h>	//inet_addr
#include<netdb.h>	//hostent
#include<stdlib.h>
#include "allConstant.h"
#include <unistd.h>
#include <thread>
#include <chrono> 
#include "quartic.h"
#include <map>
#include <sstream> 
#include <iterator>
#include <unordered_set>
#include <vector>
#include <math.h>
#include <future>
#include <complex>
#include <algorithm>
#include "Point.h"
#include "CycleTrial.h"
#include "HyperTrial.h"
#include "generalFunction.h"
#include "geometryFunction.h"
#include "EllipseTrial.h"
#include <fstream>
#include <unordered_set>

using namespace std;

class User
{
  
  private:
    string id;
    int nodes;
    int anchor;
    int radioRange;
    int epoch;
    int tSize;
    short xList[nodeSize];
    short yList[nodeSize];
    short tx[strokeSize];
    short ty[strokeSize];
    //int hv[nodeSize][anchorSize];
    bool computingLock;
    int totalStroke;
    
  public:
    unordered_set<string> TAS;
    map<string,string> resultMap;
    User():xList(),yList()
    {
        id="none";
        nodes=-1;
        anchor=-1;
        radioRange=-1;
        epoch=-1;
        tSize=0;
        computingLock=0;
        for(int i=0;i<strokeSize;i++)
        {
          tx[i]=-1;
          ty[i]=-1;
        }
        for(int i=0;i<nodeSize;i++)
        {
          xList[i]=rand()%width;
          yList[i]=rand()%height;
        }
    
    }
    User(string i):xList(),yList()
    {
        id=i;
        epoch=-1;
        nodes=-1;
        anchor=-1;
        radioRange=-1;
        tSize=0;
        computingLock=0;
        for(int i=0;i<strokeSize;i++)
        {
          tx[i]=-1;
          ty[i]=-1;
        }
        for(int i=0;i<nodeSize;i++)
        {
          xList[i]=rand()%width;
          yList[i]=rand()%height;
        }
    }
    void setNodes(int n){nodes=n;}
    void setAnchor(int a){anchor=a;}
    void setRange(int r){radioRange=r;}
    void setEpoch(int e){epoch=e;}
    void setX(int i, int v){xList[i]=v;}
    void setY(int i, int v){yList[i]=v;}
    void setTraj(int x, int y, int ind){tx[ind]=x;ty[ind]=y;}
    int getNodes(){return nodes;}
    int getAnchor(){return anchor;}
    int getRange(){return radioRange;}
    int getEpoch(){return epoch;}
    int getX(int i){return xList[i];}
    int getY(int i){return yList[i];}
    void genTAS(int totalStroke);
    void addToTAS(Point p, int width);
    void addToTAS(Point p1, Point p2, int width);
    void printArea();
    vector<twoCycleTrial> findTwoCycleTrial(short hv[][anchorSize]);
    vector<hyperTrial> findHyperTrial(short hv[][anchorSize]);
    vector<hyperTrial> hyperHelper(short hv[][anchorSize], int, int);
    vector<ellipseTrial> findEllipseTrial(short hv[][anchorSize]);
    vector<ellipseTrial>  ellipseHelper(short hv[][anchorSize],int  ah,int  h3);
    void updateStroke(short st, short st_id, short strokeWidth);
    void getHopInfo(short hv[][anchorSize]);
    void loadNetwork(string uID);
    string randomNetwork(short nodeNum, short anchorNum, short radioRange);
    //void updataHopVector();
};


void insertArea(unordered_set<int> & nodesDic, int dic)
{
  
  int shiftX[15] = {-70000,-60000,-50000,-40000,-30000,-20000,-10000,0,10000,20000,30000,40000,50000,60000,70000};
  int shiftY[15] = {-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7};
  for(int i =0;i<15; i++)
  {
    for(int j = 0; j<15; j++)
    {
      nodesDic.insert(dic+shiftX[i]+shiftY[j]);
    }
  }
}

void User:: updateStroke(short st, short st_id,short strokeWidth)
{
  ifstream strokeFile;
  strokeFile.open("../strokes/"+to_string(st)+"_"+to_string(st_id));
  if(strokeFile.good())
  {
    tx[0]=-1;
    ty[0] = strokeWidth;
    int counter = 1;
    do
    {
      strokeFile >> tx[counter];
    }while(strokeFile >> ty[counter++]);
    genTAS(st);
  }
  else
    cout<<"bad file for: "<<"../strokes/"+to_string(st)+"_"+to_string(st_id)<<endl;
    
  strokeFile.close();
//  for (int i=0; i<totalStroke; i++)
//    cout<<tx[i]<<" "<<ty[i]<<endl;
}

string User:: randomNetwork(short nodeNum, short anchorNum, short radioDis)
{

  int randID=rand();
  
  ofstream newNetwork;
  
  newNetwork.open("../results/"+to_string(randID));
  newNetwork<<nodeNum<<" "<<anchorNum<<" "<<radioDis<<endl;
  int count = 0;
  unordered_set<int> nodesDic;
  while(count<nodeNum)
  {
    int randX = rand()%width;
    int randY = rand()%height;
    int dic = 10000*randX+randY;
    unordered_set<int>::iterator got = nodesDic.find(dic);
    if(got == nodesDic.end())
    {
      insertArea(nodesDic,dic);
      count++;
      newNetwork<<randX<<" "<<randY<<endl;
    }

  }

  newNetwork.close();
  return to_string(randID);
}

void User:: loadNetwork(string uID)
{
  id = uID;
  ifstream network;
  network.open("../results/"+id);
  network>>nodes;
  network>>anchor;
  network>>radioRange;
  for (int i =0; i<nodes; i++)
  {
    network>>xList[i];
    network>>yList[i];
  
  }
  network.close();
//  for (int i=0; i<nodes; i++)
//    cout<<xList[i]<<" "<<yList[i]<<endl;

}

void User:: printArea()
{
    cout<<"printing TAS"<<endl;
    int counter=0;
    for(string point: TAS)
    {
        counter++;
        cout<<"("<<point<<");";
    }
    cout<<endl;
    cout<<"total: "<<counter<<" points"<<endl;
    


}


vector<ellipseTrial> User::findEllipseTrial(short hv[][anchorSize])
{

    vector<ellipseTrial> trials;
    
    vector<future<vector<ellipseTrial>>> VF;
    for(int i=1;i<=19;i++)
    {
      for(int j=1;j<=19;j++)
      {
        VF.push_back(async(&User::ellipseHelper,this,hv,j,i));
      }
    }
    
    for(auto& V:VF)
    {
      vector<ellipseTrial> curT = V.get();
      //cout<<"Threads results number: "<<curT.size()<<endl;
      for(vector<ellipseTrial>::iterator it = curT.begin();it != curT.end();++it)
      {
        trials.push_back(*it);
      }
    }
    
    return trials;
}




vector<ellipseTrial> User::ellipseHelper(short hv[][anchorSize],int  ah,int  h3)
{
    vector<ellipseTrial> trials;

    for(int i=0;i<anchor;i++)
    {
      for(int j=0;j<anchor;j++)
      {
        
          float avgHopDis1= radioRange;
          float avgHopDis2= radioRange;
          float dis=sqrt(euDis(xList[i],yList[i],xList[j],yList[j]));
          Point p1={xList[i],yList[i]};
          Point p2={xList[j],yList[j]};
          if(hv[i][j]!=-1&&hv[i][j]!=0)avgHopDis1=(dis + radioRange/3) / hv[i][j];
          for(int z=0;z<anchor;z++)
          {
            //cout<<"Thread is working"<<endl;
            Point p3={xList[z],yList[z]};
            float dis2=sqrt(euDis(xList[z],yList[z],xList[j],yList[j]));
            float dis3=sqrt(euDis(xList[z],yList[z],xList[i],yList[i]));
            if(hv[z][j]!=-1&&hv[z][j]!=0 && hv[z][i]!=-1&&hv[z][i]!=0)
            {
	            avgHopDis2=(dis+ dis2 +dis3+radioRange) / (hv[i][j]+hv[z][j]+hv[z][i]);
	            //cout<<"avgHopDis2: "<<avgHopDis2<<endl;
	            vector<float> yr1;
	            vector<float> yr2;
	            float area1 = findInterEllipseCycle(p1,p2,p3,ah,h3,avgHopDis1,avgHopDis2,yr1);
	            float area2 = findInterEllipseCycle(p1,p2,p3,(ah-1),h3,avgHopDis1,avgHopDis2,yr2);
	            
	            //cout<<"area1 is : "<< area1<<endl;
	            //cout<<"area2 is : "<< area2<<endl;
	        
	            if(area1-area2>1000 && area1-area2 < 640000)
	            {
	              float yroot[8]={0.0};
	              int startInd=0;
	              for(float& y1:yr1)
	              {
	                yroot[startInd++]=y1;
	              }
	              startInd=4;
	              for(float& y2:yr2)
	              {
	                yroot[startInd++]=y2;
	              }
	              ellipseTrial tryEntry={i,j,z,static_cast<float>(xList[i]),static_cast<float>(yList[i]),
	                                static_cast<float>(xList[j]),static_cast<float>(yList[j]),static_cast<float>(xList[z]),
	                                static_cast<float>(yList[z]),static_cast<float>(ah),static_cast<float>(h3),
	                                static_cast<float>(avgHopDis1),static_cast<float>(avgHopDis2),
	                                static_cast<float>(area1-area2),0,0,0};
	              memcpy(&(tryEntry.yroot), &yroot, 32) ;    
	              trials.push_back(tryEntry);
	            }
            }
          }
          
        
      }
    }
    //cout<<"trials size: "<<trials.size()<<endl;
    return trials;
}


vector<hyperTrial> User::findHyperTrial(short hv[][anchorSize])
{

    vector<hyperTrial> trials;
    
    vector<future<vector<hyperTrial>>> VF;
    for(int i=1;i<=19;i++)
    {
      for(int j=1;j<=19;j++)
      {
        VF.push_back(async(&User::hyperHelper,this,hv,j,i));
      }
    }
    
    for(auto& V:VF)
    {
      vector<hyperTrial> curT = V.get();
      //cout<<"Threads results number: "<<curT.size()<<endl;
      for(vector<hyperTrial>::iterator it = curT.begin();it != curT.end();++it)
      {
        trials.push_back(*it);
      }
    }
    
    return trials;
}






void User :: getHopInfo(short hv[][anchorSize])
{
    nodes = nodes < 5000 ? nodes : 5000; 
    nodes = nodes < 0 ? 5000 : nodes;
    vector<vector<int>> neighborsList(nodes);
    
    for(int i=0;i<nodes;i++)
    {
        int curX=xList[i];
        int curY=yList[i];
        
        for(int j=0;j<nodes;j++)
        {
            if(i!=j){
                int neighborX=xList[j];
                int neighborY=yList[j];
                if(euDis(curX,curY,neighborX,neighborY)<radioRange*radioRange)
                {
                    neighborsList[i].push_back(j);
                }
            }
        }
    }
    
    //cout<<radioRange<<endl;
    
    
    for(int i=0;i<nodes;i++)
    {
        for(int j=0;j<anchor;j++)
        {
            hv[i][j]=-1;
        }
    }
    
    for(int i=0;i<anchor;i++)
    {
        unordered_set<int> visited;
        queue<int> nq;
        nq.push(i);
        while(!nq.empty())
        {
            int cur=nq.front();
            nq.pop();
            visited.insert(cur);
            int curHop=hv[cur][i]==-1?0:hv[cur][i];
          if(!neighborsList[cur].empty()){
            for(int neighbor : neighborsList[cur])
            {
                if(visited.find(neighbor)==visited.end())
                {
                    nq.push(neighbor);
                    visited.insert(neighbor); 
                }
                if(hv[neighbor][i]==-1)hv[neighbor][i]=curHop+1;
                else
                {
                    hv[neighbor][i]=hv[neighbor][i]<curHop+1?hv[neighbor][i]:curHop+1;
                }  
            }
          }
      }
    }

}

vector<twoCycleTrial> User::findTwoCycleTrial(short hv[][anchorSize])
{

    vector<twoCycleTrial> trials;

    for(int i=0;i<anchor;i++)
    {
      for(int j=0;j<anchor;j++)
      {
      if(i!=j){
        float avgHopDis= radioRange;
        float dis=sqrt(euDis(xList[i],yList[i],xList[j],yList[j]));
        Point p1={xList[i],yList[i]};
        Point p2={xList[j],yList[j]};
        //cout<<i<<" and "<<j<<" hv: "<<hv[i][j]<<" dis "<<dis<<endl;
        if(hv[i][j]!=-1&&hv[i][j]!=0)avgHopDis=dis / hv[i][j];
      
        //cout<<"avgHopDis: "<<avgHopDis<<endl;
        for(int h1=1;h1<20;h1++)
        {
          for(int h2=1;h2<20;h2++)
          {
            if ( (h1+h2)*avgHopDis>dis)
            {
                float area=findInterTwoCycle(p1,h1,p2,h2,avgHopDis);
                if (area>1000)
                {
                    twoCycleTrial tryEntry={i,j,static_cast<float>(xList[i]),static_cast<float>(yList[i]),
                                        static_cast<float>(xList[j]),static_cast<float>(yList[j]),static_cast<float>(h1),
                                        static_cast<float>(h2),static_cast<float>(avgHopDis),static_cast<float>(area),0,0,0};  
                    trials.push_back(tryEntry);
                }
            }
          
          }
        }
       }
      }
    }
    return trials;
}

vector<hyperTrial> User::hyperHelper(short hv[][anchorSize],int  ah,int  h3)
{
    vector<hyperTrial> trials;

    for(int i=0;i<anchor;i++)
    {
      for(int j=0;j<anchor;j++)
      {
        if(i!=j){
          float avgHopDis1= radioRange;
          float avgHopDis2= radioRange;
          float dis=sqrt(euDis(xList[i],yList[i],xList[j],yList[j]));
          Point p1={xList[i],yList[i]};
          Point p2={xList[j],yList[j]};
          if(hv[i][j]!=-1&&hv[i][j]!=0)avgHopDis1=dis / hv[i][j];
          //avgHopDis2= avgHopDis1;
          if(dis>2*ah*avgHopDis1)
          {
            for(int z=0;z<anchor;z++)
            {
              //cout<<"Thread is working"<<endl;
              Point p3={xList[z],yList[z]};
              float dis2=sqrt(euDis(xList[z],yList[z],xList[j],yList[j]));
              float dis3=sqrt(euDis(xList[z],yList[z],xList[i],yList[i]));
              if(hv[z][j]!=-1&&hv[z][j]!=0 && hv[z][i]!=-1&&hv[z][i]!=0)
              {
                  avgHopDis2=(dis+ dis2 +dis3) / (hv[i][j]+hv[z][j]+hv[z][i]);
                  //cout<<"avgHopDis2: "<<avgHopDis2<<endl;
                  vector<float> yr1;
                  vector<float> yr2;
                  float area1 = findInterHyperCycle(p1,p2,p3,ah,h3,avgHopDis1,avgHopDis2,yr1);
                  float area2 = findInterHyperCycle(p1,p2,p3,(ah-1),h3,avgHopDis1,avgHopDis2,yr2);
                  
                  //cout<<"area1 is : "<< area1<<endl;
                  //cout<<"area2 is : "<< area2<<endl;
                  
                  if(area1-area2>1000 && area1-area2 < 640000){
                    float yroot[8]={0.0};
                    int startInd=0;
                    for(float& y1:yr1)
                    {
                      yroot[startInd++]=y1;
                    }
                    startInd=4;
                    for(float& y2:yr2)
                    {
                      yroot[startInd++]=y2;
                    }
                    hyperTrial tryEntry={i,j,z,static_cast<float>(xList[i]),static_cast<float>(yList[i]),
                                      static_cast<float>(xList[j]),static_cast<float>(yList[j]),static_cast<float>(xList[z]),
                                      static_cast<float>(yList[z]),static_cast<float>(ah),static_cast<float>(h3),
                                      static_cast<float>(avgHopDis1),static_cast<float>(avgHopDis2),
                                      static_cast<float>(area1-area2),0,0,0};
                    memcpy(&(tryEntry.yroot), &yroot, 32) ;    
                    trials.push_back(tryEntry);
                  }
              }
            }
          }
        }
      }
    }
    //cout<<"trials size: "<<trials.size()<<endl;
    return trials;
}

void User:: genTAS(int ts)
{
    TAS.clear();
    totalStroke=ts;
    int width=1;
    for( int i=0;i<ts;i++)
    {
        if(tx[i]<0)
        {
            width=ty[i];
        
        }
        else
        {
            if(tx[i-1]>0)
            {
                Point preP={tx[i-1],ty[i-1]};
                Point curP={tx[i],ty[i]};
                addToTAS(preP,curP,width);
            
            }
            else
            {
                Point curP={tx[i],ty[i]};
                addToTAS(curP,width);
            
            }
        
        }
    
    }

}

void User:: addToTAS(Point p, int width)
{
    int lx= p.m_X-width / 2 < 0? 0 : p.m_X-width / 2;
    int ly= p.m_Y-width / 2 < 0? 0 : p.m_Y-width / 2;
    for( int i=lx; i<p.m_X+width / 2;i++)
    {
        for(int j=ly;j<p.m_Y+width / 2;j++)
        {
            Point curP={i,j};
            if(pDis(curP, p)<width*width / 4)TAS.insert(to_string(curP.m_X)+" "+to_string(curP.m_Y));
        }
    }


}


void User:: addToTAS(Point p1, Point p2,int width)
{
    int lx1= p1.m_X-width / 2<0?0:p1.m_X-width / 2;
    int ly1= p1.m_Y-width / 2<0?0:p1.m_Y-width / 2;
    int lx2= p2.m_X-width / 2<0?0:p2.m_X-width / 2;
    int ly2= p2.m_Y-width / 2<0?0:p2.m_Y-width / 2;
    int rx= p1.m_X+width / 2<p2.m_X+width / 2?p2.m_X+width / 2:p1.m_X+width / 2;
    int ry= p1.m_Y+width / 2<p2.m_Y+width / 2?p2.m_Y+width / 2:p1.m_Y+width / 2;
    int lx=lx1<lx2?lx1:lx2;
    int ly=ly1<ly2?ly1:ly2;
    
    for( int i=lx; i<rx;i++)
    {
        for(int j=ly;j<ry;j++)
        {
            Point curP={i,j};
            if(lDis(curP, p1,p2)<width*width / 4)TAS.insert(to_string(curP.m_X)+" "+to_string(curP.m_Y));
        }
    }
}

#endif
