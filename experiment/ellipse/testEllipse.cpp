

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "Point.h"
#include <math.h>
#include "generalFunction.h"
#include "geometryFunction.h"
#include <iostream>
#include <string>
#include <unordered_set>
#include <sstream>
#include <map>
#include <iterator>
#include <future>
#include <complex>
#include <algorithm>
#include "EllipseTrial.h"

using namespace std;
bool checkInEllipse ( Point lf, Point rf, Point c, Point cur, float a , float r)
{

  return sqrt(pDis(lf,cur))+sqrt(pDis(rf,cur)) <= 2*a && pDis(c, cur) <=  r*r;


}

int main()
{
  srand(time(NULL));
  const int ITERATION=100000;
  const int WIDTH = 1000;
  const int HEIGHT = 1000;
  const float avgDis = 80;
  int count =0;
  for(int i=0; i<ITERATION; i++){
    Point lf={rand()%WIDTH+1000, rand()%HEIGHT+1000};
    Point rf={rand()%WIDTH+1000, rand()%HEIGHT+1000};
    Point c={rand()%WIDTH+1000, rand()%HEIGHT+1000};
    
    int ah = rand()%10+1;
    int h3 = rand()%10+1;
    float al = ah*avgDis;
    float rl = h3*avgDis;
    vector<float> yr;
    float area = findInterEllipseCycle(lf,rf,c, ah,h3,avgDis,avgDis, yr);
    float testArea=0;
    for(int w = 0; w< 3000; w++){
      for(int h = 0;h< 3000; h++){
        Point cur = {w, h};
        if ( checkInEllipse(lf,rf,c, cur, al, rl)) testArea+=1;
      
      } 
    }
    count++;
    cout<<count<<" round ";
    if(abs(testArea - area)/testArea >0.1 && abs(testArea - area) >0.05*3.1415926*rl*rl && abs(testArea - area) >0.05*3.1415926*al*al && area >3000){
      cout<<"failed********************************************************************************"<<endl;
      cout<<"testArea : "<<testArea<<"; Math area: "<< area<<endl;//"; circle area: "<< M_PI * rl *rl<<"; dif: "<< M_PI * rl *rl - area << endl;
      cout<<"lf: "<<lf<<"rf: "<<rf<<"c: "<<c<<"a: "<< ah*avgDis<<endl;
      cout<<"r: "<<h3*avgDis<<endl;
      int vsize = yr.size();
      for(int j=0;j <vsize;j++)
	    {
        cout<<yr[j]<<" ";
	    }
	    cout<<endl;
	    return 1;
    }
    else cout<<"success"<<endl;
  }
  

}
