#ifndef POINT_H
#define POINT_H
#include <iostream>
#include <string>
#include <unordered_set>
#include <math.h>
#include <sstream>
using namespace std;

struct Point
{
  int m_X;
  int m_Y;

};

ostream & operator<<(ostream& os, Point p){
  os<<"("<<p.m_X<<", "<<p.m_Y<<")"<<endl;
  return os;
}

#endif
