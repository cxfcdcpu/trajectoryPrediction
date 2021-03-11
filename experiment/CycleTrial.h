#ifndef CIRCLETRIAL_H
#define CIRCLETRIAL_H
#include <iostream>
#include <string>
#include <unordered_set>
#include <math.h>
#include <sstream>
#include "Point.h"
using namespace std;


struct twoCycleTrial
{
  int c1ID;
  int c2ID; 
  float c1X;
  float c1Y;
  float c2X;
  float c2Y;
  float h1;
  float h2;
  float d;
  float tArea;
  float grAr;
  float rate3;
  float acAr;
};




ostream& operator<<(ostream& os,twoCycleTrial &rhs)
{
  os<<"best two cycle:"<<endl;
  os<<"c1ID: "<<rhs.c1ID<<" coor: ("<<rhs.c1X<<", "<<rhs.c1Y<<")"<<endl;
  os<<"c2ID: "<<rhs.c2ID<<" coor: ("<<rhs.c2X<<", "<<rhs.c2Y<<")"<<endl;
  os<<"r1 and r2 and d: "<<rhs.h1<<"  "<<rhs.h2<<"  "<< rhs.d<< endl;
  os<<"grAr and acAr: "<<rhs.grAr<<"  "<<rhs.acAr<<endl;
  os<<"tArea and rate3: "<<rhs.tArea<< "   "<<rhs.rate3<<endl;
  return os;
}



bool sortCycleTrial(twoCycleTrial ct1, twoCycleTrial ct2)
{
  return ct1.grAr * ct1 . rate3 > ct2 . grAr * ct2 . rate3;

}







#endif
