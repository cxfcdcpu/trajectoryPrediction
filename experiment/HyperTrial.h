#ifndef HYPERTRIAL_H
#define HYPERTRIAL_H
#include <iostream>
#include <string>
#include <unordered_set>
#include <math.h>
#include <sstream>
#include "Point.h"
using namespace std;


struct hyperTrial
{
  int c1ID;
  int c2ID;
  int c3ID; 
  float c1X;
  float c1Y;
  float c2X;
  float c2Y;
  float c3X;
  float c3Y;  
  float ah;
  float h3;
  float avgD1;
  float avgD2;
  float tArea;
  float grAr;
  float rate3;
  float acAr;
  float yroot[8];
};


ostream& operator<<(ostream& os,hyperTrial &rhs)
{
  os<<"best HyperBola and cycle:"<<endl;
  os<<"c1ID: "<<rhs.c1ID<<" coor: ("<<rhs.c1X<<", "<<rhs.c1Y<<")"<<endl;
  os<<"c2ID: "<<rhs.c2ID<<" coor: ("<<rhs.c2X<<", "<<rhs.c2Y<<")"<<endl;
  os<<"c3ID: "<<rhs.c3ID<<" coor: ("<<rhs.c3X<<", "<<rhs.c3Y<<")"<<endl;
  os<<"ah and h3 and avgD1 and avgD2: "<<rhs.ah<<"  "<<rhs.h3<<"  "<< rhs.avgD1<<" "<< rhs.avgD2<< endl;
  os<<"grAr and acAr: "<<rhs.grAr<<"  "<<rhs.acAr<<endl;
  os<<"tArea and rate3: "<<rhs.tArea<< "   "<<rhs.rate3<<endl;
  os<<"roots are: ";
  for(int kk=0;kk<8;kk++)
  {
    os<< rhs.yroot[kk]<< " ";
  }
  os<<endl;
  return os;
}

bool sortHyperTrial(hyperTrial ct1, hyperTrial ct2)
{
  return  ct1.grAr * ct1 . rate3 > ct2 . grAr * ct2 . rate3;

}


#endif

