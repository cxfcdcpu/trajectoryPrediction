#ifndef GENERALFUNC_H
#define GENERALFUNC_H
#include <iostream>
#include <string>
#include <unordered_set>
#include <math.h>
#include <sstream>
#include "Point.h"
#include "CycleTrial.h"
#include "HyperTrial.h"
#include "EllipseTrial.h"
#include <vector>


void updateTAS(twoCycleTrial ct, unordered_set<string> & st);
void updateTAS(ellipseTrial ct, unordered_set<string> & st);
void updateTAS(hyperTrial ht, unordered_set<string> & st);
string findBestTry(twoCycleTrial * cTri, hyperTrial * hTri, int cs, int hs, unordered_set<string> & TAS);
int myStoi(string data);
int pDis(Point p1, Point p2);
int lDis(Point p1, Point p2, Point p3);
string findBestTry2(ellipseTrial * eTri, hyperTrial * hTri, int cs, int hs, unordered_set<string> & TAS);








void updateTAS(twoCycleTrial ct, unordered_set<string> & TAS)
{
  vector<string> temp;
  for(string st : TAS)
  {
    stringstream ss(st);
    float i,j;
    ss>>i;
    ss>>j;
    //<<i<<endl;
    float di1=ct.c1X-i;
    float dj1=ct.c1Y-j;
    float di2=ct.c2X-i;
    float dj2=ct.c2Y-j;
    float rr1 = ct.h1*ct.d*ct.h1*ct.d;
    float rr2 = ct.h2*ct.d*ct.h2*ct.d;
    float rr3 = (ct.h1-1)*ct.d*(ct.h1-1)*ct.d;
    
    if (di1*di1+dj1*dj1<=rr1 && di1*di1+dj1*dj1>rr3 && di2*di2+dj2*dj2<=rr2)
    {
      //cout<<"erasing cycle tas"<<endl;
      temp.push_back(st);
    }
  }
  for(string st : temp)
  {
      TAS.erase(st);
  }
}
  
void updateTAS(hyperTrial ct, unordered_set<string> & TAS)
{
  vector<string> temp;
  for(string st : TAS)
  {
    stringstream ss(st);
    float i,j;
    ss>>i;
    ss>>j;
    //cout<<j<<endl;
    float tt = sqrt((i-ct.c1X)*(i-ct.c1X)+(j-ct.c1Y)*(j-ct.c1Y))-sqrt((i-ct.c2X)*(i-ct.c2X)+(j-ct.c2Y)*(j-ct.c2Y));
    float di=ct.c3X-i;
    float dj=ct.c3Y-j;
    float rr = ct.h3*ct.avgD2*ct.h3*ct.avgD2;
    float a = ct.ah * ct.avgD1;
    float a2 = (ct.ah-1.0)*ct.avgD1;
    
    if(di * di + dj * dj <= rr && tt <= 2*a && tt >= 2* a2)
    {
      //cout<<"erasing hyper tas"<<endl;
      temp.push_back(st);
    }
  }
  for(string st : temp)
  {
    TAS.erase(st);
  }
}


void updateTAS(ellipseTrial ct, unordered_set<string> & TAS)
{
  vector<string> temp;
  for(string st : TAS)
  {
    stringstream ss(st);
    float i,j;
    ss>>i;
    ss>>j;
    //cout<<j<<endl;
    float tt = sqrt((i-ct.c1X)*(i-ct.c1X)+(j-ct.c1Y)*(j-ct.c1Y))+sqrt((i-ct.c2X)*(i-ct.c2X)+(j-ct.c2Y)*(j-ct.c2Y));
    float di=ct.c3X-i;
    float dj=ct.c3Y-j;
    float rr = ct.h3*ct.avgD2*ct.h3*ct.avgD2;
    float a = ct.ah * ct.avgD1;
    float a2 = (ct.ah-1.0)*ct.avgD1;
    
    if(di * di + dj * dj <= rr && tt <= 2*a && tt >= 2* a2)
    {
      //cout<<"erasing hyper tas"<<endl;
      temp.push_back(st);
    }
  }
  for(string st : temp)
  {
    TAS.erase(st);
  }
  cout<<"done update TAS for ellipse and circle."<<endl;
  
  
}

int pDis(Point p1, Point p2)
{
    return (p1.m_X-p2.m_X)*(p1.m_X-p2.m_X)+(p1.m_Y-p2.m_Y)*(p1.m_Y-p2.m_Y);
}

int lDis(Point p0, Point p1, Point p2)
{
    if(p2.m_X!=p1.m_X && p2.m_Y != p1.m_Y)
        return ((p2.m_X-p1.m_X)*(p1.m_Y-p0.m_Y)-(p1.m_X-p0.m_X)*(p2.m_Y-p1.m_Y))*((p2.m_X-p1.m_X)*(p1.m_Y-p0.m_Y)-(p1.m_X-p0.m_X)*(p2.m_Y-p1.m_Y))/((p2.m_X-p1.m_X)*(p2.m_X-p1.m_X)+(p2.m_Y-p1.m_Y)*(p2.m_Y-p1.m_Y));
    else
        return sqrt(pDis(p0,p1));
}


string findBestTry(twoCycleTrial * cTri, hyperTrial * hTri, int cs, int hs, unordered_set<string> & TAS)
{
    string res="";
    int maxInd=0;
    float max1=0;
    for(int k=0;k<cs;k++)
    {
        if(max1<cTri[k].grAr){
          max1=cTri[k].grAr;
          maxInd=k;
        }
    }
    

    int maxInd2=0;
    float max2=0;
    for(int k=0;k<hs;k++)
    {
        if(max2<hTri[k].grAr){
          max2=hTri[k].grAr;
          maxInd2=k;
        }
    }
    
    if(max1>=max2)
    {
        cout<<cTri[maxInd]<<endl;
        
        updateTAS(cTri[maxInd],TAS);
        
        res+="\n(";
        res+=to_string(cTri[maxInd].c1ID);
        res+="; ";
        res+=to_string(cTri[maxInd].c2ID);
        res+="; ";
        res+=to_string(cTri[maxInd].h1);
        res+="; ";
        res+=to_string(cTri[maxInd].h2);
        res+="; ";
        res+=to_string(cTri[maxInd].d);
        res+=")";
    
    }
    else
    {
        cout<<hTri[maxInd2]<<endl;
        updateTAS(hTri[maxInd2],TAS);
        res+="\n(";
        res+=to_string(hTri[maxInd2].c1ID);
        res+="; ";
        res+=to_string(hTri[maxInd2].c2ID);
        res+="; ";
        res+=to_string(hTri[maxInd2].c3ID);
        res+="; ";
        res+=to_string(hTri[maxInd2].ah);
        res+="; ";
        res+=to_string(hTri[maxInd2].h3);
        res+="; ";
        res+=to_string(hTri[maxInd2].avgD1);
        res+="; ";
        res+=to_string(hTri[maxInd2].avgD2);
        res+=")";
    }
    

    return res;
}


string findBestTry2(ellipseTrial * eTri, hyperTrial * hTri, int cs, int hs, unordered_set<string> & TAS)
{
    string res="";
    int maxInd=0;
    float max1=0;
    for(int k=0;k<cs;k++)
    {
        if(max1<eTri[k].grAr){
          max1=eTri[k].grAr;
          maxInd=k;
        }
    }
    

    int maxInd2=0;
    float max2=0;
    for(int k=0;k<hs;k++)
    {
        if(max2<hTri[k].grAr){
          max2=hTri[k].grAr;
          maxInd2=k;
        }
    }
    
    if(max1>=max2)
    {
        cout<<eTri[maxInd]<<endl;
        
        updateTAS(eTri[maxInd],TAS);
        
        //cout<<"finish update TAS"<<endl;
        
        res+="\n(";
        res+=to_string(eTri[maxInd].c1ID);
        res+="; ";
        res+=to_string(eTri[maxInd].c2ID);
        res+="; ";
        res+=to_string(eTri[maxInd].c3ID);
        res+="; ";
        res+=to_string(-eTri[maxInd].ah);
        res+="; ";
        res+=to_string(eTri[maxInd].h3);
        res+="; ";
        res+=to_string(eTri[maxInd].avgD1);
        res+="; ";
        res+=to_string(eTri[maxInd].avgD2);
        res+=")";
    
    }
    else
    {
        cout<<hTri[maxInd2]<<endl;
        updateTAS(hTri[maxInd2],TAS);
        //cout<<"finish update TAS"<<endl;
        
        res+="\n(";
        res+=to_string(hTri[maxInd2].c1ID);
        res+="; ";
        res+=to_string(hTri[maxInd2].c2ID);
        res+="; ";
        res+=to_string(hTri[maxInd2].c3ID);
        res+="; ";
        res+=to_string(hTri[maxInd2].ah);
        res+="; ";
        res+=to_string(hTri[maxInd2].h3);
        res+="; ";
        res+=to_string(hTri[maxInd2].avgD1);
        res+="; ";
        res+=to_string(hTri[maxInd2].avgD2);
        res+=")";
    }
    

    return res;
}

int myStoi(string data)
{
     stringstream geek(data); 
    
    // The object has the value 12345 and stream 
    // it to the integer x 
    int x = 0; 
    geek >> x;
    //cout<<"transfer string "<<data<<"to integer "<<x<<endl; 
    return x;

}



#endif
