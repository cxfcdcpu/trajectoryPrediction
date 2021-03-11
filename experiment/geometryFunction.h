#ifndef GEOMETRYFUNCTION_H
#define GEOMETRYFUNCTION_H
#include <iostream>
#include <string>
#include <unordered_set>
#include <math.h>
#include <sstream>
#include "Point.h"
#include "CycleTrial.h"
#include "HyperTrial.h"
#include "generalFunction.h"
#include <complex>
#include "quartic.h"
#include <algorithm>
struct Point2
{
  double m_X;
  double m_Y;

};

ostream & operator<<(ostream& os, Point2 p){
  os<<"("<<p.m_X<<", "<<p.m_Y<<")"<<endl;
  return os;
}

Point2 operator -(Point2 &p1, Point2 &p2)
{
  Point2 temp = {p1.m_X-p2.m_X, p1.m_Y-p2.m_Y};
  return temp;

}

int euDis(int x1,int y1,int x2,int y2);
double euDis2(double x1,double y1,double x2,double y2);
double getAngle(Point2 p);
float findInterTwoCycle(Point,int,Point,int,float);
float findInterHyperCycle(Point,Point,Point, int,int,float,float,vector<float>&);
vector<tuple<float,float>> getRoot(float x0,float y0,float r,float a,float b,vector<float>&);
float rootX(float,float,float);
float hyperArea(float a,float b,tuple<float,float> z1,tuple<float,float> z2);
vector<float> realQudricRoot(float b,float c,float d,float e);
float hyperAreaHelper(float a,float b,float m,float z,float y);
vector<Point2> ellipseSpectialCase1(double a,double b,double r);
float findInterEllipseCycle(Point2 h1,Point2 h2,Point2 c1,int ah,int h3,float avgHopDis1,float avgHopDis2,vector<float>& yr );
vector<Point2> getRootEllipse(float cx,float cy,float r,float a, float b, float e);
double getCircleSegArea(Point2 p1,Point2 p2, Point2 orig, double r);
double getEllipseSegArea(Point2 p1,Point2 p2, Point2 orig,double a, double b);
double ellipseSectionArea(Point2 p1,Point2 p2,Point2 orig, double a,double b);
double getPolarAngle(Point2 p, double a,double  b);
double triangleArea(Point2 start, Point2 end, Point2 orig);
double getDAngle(Point2 start, Point2 end, Point2 orig);
double getQuadriArea(Point2 p1,Point2 p2,Point2 p3,Point2 p4);

int euDis(int x1, int y1, int x2, int y2)
{
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

double euDis2(double x1,double y1,double x2,double y2){
  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  
}

double getAngle(Point2 p)
{
  double dis = euDis2(0,0, p.m_X, p.m_Y);  
  if(p.m_Y>=0)return acos(p.m_X/dis);
  else return 2*M_PI-acos(p.m_X/dis);
}

double getDAngle(Point2 start, Point2 end, Point2 orig)
{
  double angle1 = getAngle(start-orig);
  double angle2 = getAngle(end-orig);
  
  //cout<<"angle1: "<<angle1<<" ;angle2: "<<angle2<<endl;
  return angle2-angle1>=0?angle2-angle1: 2*M_PI + angle2-angle1;



}

double triangleArea(Point2 start, Point2 end, Point2 orig)
{
  
  double area = abs(start.m_X*(end.m_Y-orig.m_Y)+end.m_X*(orig.m_Y-start.m_Y)+orig.m_X*(start.m_Y-end.m_Y))/2;
  //cout<<"triangle area: "<<area<<endl;
  return area;
}

double getPolarAngle(Point2 p, double a,double  b)
{
  return p.m_Y>=0?acos(p.m_X/a):(2*M_PI-acos(p.m_X/a));

}

double ellipseSectionArea(Point2 p1,Point2 p2,Point2 orig, double a,double b)
{
  Point2 temp1 = p1 - orig;
  Point2 temp2 = p2 - orig;
  double theta1 = getPolarAngle(temp1, a, b);
  double theta2 = getPolarAngle(temp2, a, b);
  
  double deltaTheta = theta1<theta2?theta2- theta1:2*M_PI+theta2-theta1;
  
  return deltaTheta*a*b/2;

}

double getEllipseSegArea(Point2 p1,Point2 p2, Point2 orig,double a, double b)
{
  double dAngleE = getDAngle(p1, p2, orig);
  double signE = dAngleE>M_PI?1.0:-1.0;
  double triArea = signE*triangleArea(p1,p2,orig);
  double elSecArea = ellipseSectionArea(p1,p2,orig, a, b);
  //cout<<"ellipse triArea: "<<triArea<<" ellipse Section Area: "<<elSecArea<<endl;
  double areaE =  elSecArea + triArea;
  
  return areaE;
}

double getCircleSegArea(Point2 p1,Point2 p2, Point2 orig, double r)
{
  double dAngle = getDAngle(p1, p2, orig);
  double sign =dAngle>M_PI?1.0:-1.0;
  
  double area = dAngle*r*r/2 + sign*triangleArea(p1,p2,orig);
  //cout<<"p1: "<<p1<<"p2: "<<p2<<"orig: "<<orig<<endl;
  //cout<<"dAngle: "<<dAngle<<"; tri area: "<<triangleArea(p1,p2,orig)<<endl;
  return area;
}

double getQuadriArea(Point2 p1,Point2 p2,Point2 p3,Point2 p4)
{
  return abs((p3.m_X-p1.m_X)*(p4.m_Y-p2.m_Y)-(p4.m_X-p2.m_X)*(p3.m_Y-p1.m_Y))/2;


}

vector<Point2> ellipseSpectialCase1(double a,double b,double r){

  vector<Point2> result;
  if(r<a && r>b){
    double y= sqrt((a*a*b*b-r*r*b*b)/(a*a-b*b));
    double x = sqrt((a*a*r*r-a*a*b*b)/(a*a-b*b));
    Point2 p1 = {x,y};
    Point2 p2 = {x,-y};
    Point2 p3 = {-x,y};
    Point2 p4 = {-x,-y};
    result.push_back(p1);
    result.push_back(p2);
    result.push_back(p3);
    result.push_back(p4);
    
  }
  return result;
}

void sortWithAngle(vector<Point2> &vp){
  for(vector<Point2>::iterator it = vp.begin(); it != vp.end(); ++it){
		for(vector<Point2>::iterator it2 = vp.begin(); it2 != vp.end()-1; ++it2){
		
		  //cout<<"p: "<<*it<<" angle: "<<getAngle(*it)<<endl;
		  	 
		  if(getAngle(*it2)>getAngle(*(it2+1))){
		    
		    Point2 temp = *(it2+1);
		    *(it2+1) = *it2;
		    *it2 = temp;
		  }   
		   
		}		
  }

}


float circleOverlapArea(Point p1,int h1,Point p2,int h2,float avgHopDis)
{
    const float pi=3.141592653589793;
    float d=sqrt(pDis(p1,p2));
    float A1=0;
    float r1=h1*avgHopDis;
    float r2=h2*avgHopDis;
    
    if(d>=r1+r2)return 0;
    else if(r2>=r1+d)return pi*(r1*r1);
    else if(r1>r2+d)A1=pi*r2*r2;
    
    if(A1==0)
    {
      if(d > 0 && r2 > 0 )
        A1=r2*r2*acos((d*d+r2*r2-r1*r1) / (2*d*r2))+r1*r1*acos((d*d+r1*r1-r2*r2) / (2*d*r1))-0.5*sqrt((-d+r2+r1)*(d+r2-r1)*(d-r2+r1)*(d+r1+r2));
      else
        A1=10000;
    }
    
    
    return A1;

}




float findInterTwoCycle(Point p1,int h1,Point p2,int h2,float avgHopDis)
{
    const float pi=3.141592653589793;
    float d=sqrt(pDis(p1,p2));
    float A1=0;
    float r1=h1*avgHopDis;
    float r2=h2*avgHopDis;
    float r3=(h1-1)*avgHopDis;
    
    if(d>=r1+r2)return 0;
    else if(r2>=r1+d)return pi*(r1*r1-r3*r3);
    else if(r1>r2+d)A1=pi*r2*r2;
    
    if(A1==0)
    {
      if(d > 0 && r2 > 0 && r3 >0)
        A1=r2*r2*acos((d*d+r2*r2-r1*r1) / (2*d*r2))+r1*r1*acos((d*d+r1*r1-r2*r2) / (2*d*r1))-0.5*sqrt((-d+r2+r1)*(d+r2-r1)*(d-r2+r1)*(d+r1+r2));
      else
        A1=10000;
    }
    
    if (r3==0)return A1;
    float A2=0;
    if(d>=r3+r2)return A1;
    else if(r2>=r3+d)return A1-pi*r3*r3;
    else if(r3>=r2+d)return 0;
    if(d > 0 && r2 > 0 && r3 >0)
        A2=r2*r2*acos((d*d+r2*r2-r3*r3) / (2*d*r2))+r3*r3*acos((d*d+r3*r3-r2*r2) / (2*d*r3))-0.5*sqrt((-d+r2+r3)*(d+r2-r3)*(d-r2+r3)*(d+r3+r2));
    else
        A2=-10000;
    return A1-A2;

}

float findInterHyperCycle(Point h1,Point h2,Point c1,int ah,int h3,float avgHopDis1,float avgHopDis2,vector<float>& yr )
{
  const float pi=3.141592653589793;
  float a = ah*avgHopDis1;
  float r = h3*avgHopDis2;
  float c = sqrt(pDis(h1,h2))/2;
  float b = sqrt( c * c - a * a);
  
  float xs = -(h1.m_X + h2.m_X)/2;
  float ys = -(h1.m_Y + h2.m_Y)/2;
  bool flip = false;
  float ctheta = 1;
  float stheta = 0;
  
  /*
  float stheta = -c*(h2.m_Y+ys)/((h2.m_X+xs)*(h2.m_X+xs)+(h2.m_Y+ys)*(h2.m_Y+ys));
  //not sure if this correct.
  if ( h2.m_Y+ys != 0){
    flip=true;
    ctheta=-(h2.m_X+xs)*stheta/(h2.m_Y+ys);
  }
  */
  //line
  if(h2.m_X+xs==0 && h2.m_Y+ys>0){
    ctheta=0;
    stheta=-1;
  }
  else if(h2.m_X+xs<0 && h2.m_Y+ys==0){
    ctheta=-1;
    stheta=0;
  }
  else if(h2.m_X+xs==0 && h2.m_Y+ys<0){
    ctheta=0;
    stheta=1;
  }
  else if(h2.m_X+xs>0 && h2.m_Y+ys==0){
    ctheta=1;
    stheta=0;
  }
  else {

    ctheta=(h2.m_X+xs)/c;
    stheta=-(h2.m_Y+ys)/c;
  }
              
  float cx=c1.m_X+xs;
  float cy=c1.m_Y+ys;
  float c2x=ctheta*cx-stheta*cy;
  float c2y=stheta*cx+ctheta*cy  ;
  if (abs(a)<5 )
  {
    //cout<<"c2x: "<<c2x<<"; r: "<<r<<endl;
    if(r<abs(c2x))
    {
        yr.push_back(-12345);
        yr.push_back(c2x);
        
        if(c2x>0)return 0;
        else{
         return pi*r*r;
         yr.push_back(pi*r*r);
       }
    }
    else
    {
      //cout<<"0 area: "<<2*acos(c2x/r)<<endl;
      yr.push_back(-100);
      if(flip)yr.push_back(ctheta);
      else yr.push_back(100);
      yr.push_back(c2x);
      float totalArea=acos(c2x/r)*r*r-c2x*sqrt(r*r-c2x*c2x);
      yr.push_back(pi*r*r);
      return totalArea;
    }
  }
  
  vector<tuple<float,float>> rrr = getRoot(c2x, c2y, r, a, b,yr);
  /*
  short count=0;
  cout<<"root numbers: "<<rrr.size()<<endl;
  for( auto& tp : rrr)
  {
     cout<<"root "<<count++<<" is: ("<<get<0>(tp)<<", "<<get<1>(tp)<<")"<<endl;
  
  }
  */
  //handle size=3 later if we have time
  if(rrr.size()>=3)return 400000000;
  
  
  if(rrr.size()<2 && sqrt(pDis(c1,h1))-sqrt(pDis(c1,h2))>2*a)return 0;
  else if(rrr.size()<2 && sqrt(pDis(c1,h1))-sqrt(pDis(c1,h2))<=2*a)return pi*r*r;
  else 
  {
    float area=0;
    
    
        tuple<float,float> z1 = rrr[0];
        tuple<float,float> z2 = rrr[1];
        float i1x = get<0>(z1) - c2x;
        float i1y = get<1>(z1) - c2y;
        float i2x = get<0>(z2) - c2x;
        float i2y = get<1>(z2) - c2y; 
        float sphi = (i1x*i2y-i1y*i2x)/(i1y*i1y+i1x*i1x);
        float cphi = (i2x+i1y*sphi)/i1x*0.999999;
        
        if ((cphi>1 || cphi<-1)&& (1-sphi*sphi>=0))
            cphi=abs(cphi)/cphi*sqrt(1-sphi*sphi);

        float curveArea=hyperArea(a,b,z1,z2);
        if (sphi<0)
            area+=acos(cphi)*r*r/2-abs(sphi)*r*r/2-curveArea;            
        else
            area+=(2*pi-acos(cphi))*r*r/2+abs(sphi)*r*r/2-curveArea;

    return area;
  }
  return 0;
}




float findInterEllipseCycle(Point h1,Point h2,Point c1,int ah,int h3,float avgHopDis1,float avgHopDis2,vector<float>& yr )
{

  

  float a = ah*avgHopDis1;
  float r = h3*avgHopDis2;
  float c = sqrt(pDis(h1,h2))/2;
  
  
  
  
  if( a<c )return 0;
  
  float b = sqrt( a * a - c*c);
  
  float xs = -(h1.m_X + h2.m_X)/2;
  float ys = -(h1.m_Y + h2.m_Y)/2;

  if(c<10){
    Point tempCircleCenter = {(int)-xs,(int) -ys};
    //cout<<"!!!!!!ellipse are too close to have good results, So treat ellipse to be a circle and calculate!!!!!!"<<endl;
    return circleOverlapArea(tempCircleCenter, ah, c1, h3, avgHopDis1);
  }

  float ctheta = 1;
  float stheta = 0;
  

  //line
  if(h2.m_X+xs==0 && h2.m_Y+ys>0){
    ctheta=0;
    stheta=-1;
  }
  else if(h2.m_X+xs<0 && h2.m_Y+ys==0){
    ctheta=-1;
    stheta=0;
  }
  else if(h2.m_X+xs==0 && h2.m_Y+ys<0){
    ctheta=0;
    stheta=1;
  }
  else if(h2.m_X+xs>0 && h2.m_Y+ys==0){
    ctheta=1;
    stheta=0;
  }
  else {

    ctheta=(h2.m_X+xs)/c;
    stheta=-(h2.m_Y+ys)/c;
  }
              
  float cx=c1.m_X+xs;
  float cy=c1.m_Y+ys;
  float c2x=ctheta*cx-stheta*cy;
  float c2y=stheta*cx+ctheta*cy  ;

  
  //cout<<"e: "<< c<< " py: "<<b<<" cx: "<<c2x<<" cy: "<<c2y<<" r: "<<r<<endl;
  
  vector<Point2> result = getRootEllipse(c2x, c2y, r, a, b, c);
  
  
  
  
  
  int vsize = result.size();
  
  if(vsize>0)
  {
    for(int i=0;i <vsize;i++)
	  {
	    float nx = result[i].m_X*ctheta + stheta*result[i].m_Y -xs;
	    float ny = ctheta*result[i].m_Y - stheta*result[i].m_X -ys;
	    yr.push_back(nx);
	    yr.push_back(ny);
	  }
  }
  
  /*
	cout<<"there are "<<vsize<<" of roots"<<endl;
	for(int i=0;i <vsize;i++)
	{
	  cout<<"Point "<<i<<": "<<result[i]<<"; ";
	}
	cout<<endl;
	*/
	
	
	
	double ex =0;
	double ey =0;
	double lx = -c;
	double rx = c;
	double ly =0;
	double ry =0;
	double area =0;
  Point2 newOrig = {c2x, c2y};
  Point2 orig = {0.0,0.0};
	if(vsize==0 || vsize == 1)
	{
	  
	  //case 3,6 seperate
	  if(euDis2(ex,ey,c2x,c2y)>=r && euDis2(lx,ly, c2x,c2y)+euDis2(rx,ry,c2x,c2y)>=2*a){
	    area = 0;
	  
	  }
	  //one contained another.
	  else{
	    area = min(M_PI*a*b, M_PI*r*r);
	    
	  }
	  
	  yr.push_back(-123.45678);
	
	
	}
	else{
	  //case 2,5,9 which has two intersection points.
	  if(vsize ==2 ){
      
      float aa1 = getEllipseSegArea(result[0],result[1],orig,a,b);
      float aa2 = getCircleSegArea(result[0],result[1],newOrig, r);
      float aa3 = getEllipseSegArea(result[1],result[0],orig,a,b);
      float aa4 = getCircleSegArea(result[1],result[0],newOrig, r);
      
      yr.push_back(aa1);
      yr.push_back(aa2);
      yr.push_back(aa3);
      yr.push_back(aa4);

	    area = min(aa1,aa2);
	    area += min(aa3,aa4);
	    //cout<<"area: "<<area<<endl;
	  
	  }
	  //has four intersection points.
	  else if(vsize ==4){
	    double quadArea = getQuadriArea(result[0],result[1],result[2],result[3]);
	    area = quadArea;
	    for(int j=0;j<4;j++)
	    {
	      area+= min(getEllipseSegArea(result[j],result[(j+1)%4],orig,a,b),getCircleSegArea(result[j],result[(j+1)%4],newOrig, r));
	    }
	    //cout<<"quadArea: "<< quadArea<<endl; 
	    //cout<<"area: "<<area<<endl;
	
	  }
	}
  
  return area;
}

vector<Point2> getRootEllipse(float cx,float cy,float r,float a, float b, float e)
{
    

		//cout<<"a: "<<a<<" b: "<<b<<" r: "<<r<<endl;
		
		double m=(b*b-a*a);
		double n = a*a*b*b+b*b*cx*cx+b*b*cy*cy -b*b*r*r;
		//cout<<"m: "<<m<<" n: "<<n<<endl;
		double p1 = -4*b*b*cy/m;
		double p2 = (2*m*n+4*b*b*b*b*cy*cy+4*b*b*a*a*cx*cx)/m/m;
		double p3 = -4*b*b*cy*n/m/m;
		double p4 = (n*n-4*b*b*cx*cx*a*a*b*b)/m/m;
  
    vector<Point2> result;
		if(abs(cx)<0.01 && abs(cy)<0.01){

		   result = ellipseSpectialCase1(a,b,r);
		   sortWithAngle(result);
		   
		   /*
		   for(vector<Point2>::iterator it = result.begin(); it != result.end(); ++it){
		      cout<<*it;		   
		   }
       */
		
		}
		else{
		    complex<double>*  solutions = solve_quartic(p1, p2, p3, p4);
/*		    
		    		std::cout << "x1 = " << (solutions[0].real()>=0. ? " " : "") << solutions[0].real(); if(solutions[0].imag()!=0.0) std::cout << "   +   i * " <<  solutions[0].imag(); std::cout << std::endl;
		std::cout << "x2 = " << (solutions[1].real()>=0. ? " " : "") << solutions[1].real(); if(solutions[1].imag()!=0.0) std::cout << "   -   i * " << -solutions[1].imag(); std::cout << std::endl;
		std::cout << "x3 = " << (solutions[2].real()>=0. ? " " : "") << solutions[2].real(); if(solutions[2].imag()!=0.0) std::cout << "   +   i * " <<  solutions[2].imag(); std::cout << std::endl;
		std::cout << "x4 = " << (solutions[3].real()>=0. ? " " : "") << solutions[3].real(); if(solutions[3].imag()!=0.0) std::cout << "   -   i * " << -solutions[3].imag(); std::cout << std::endl;

		// control / test
		std::cout << std::endl;
		std::cout << polinom_4(solutions[0], p1, p2, p3, p4) << std::endl;
		std::cout << polinom_4(solutions[1], p1, p2, p3, p4) << std::endl;
		std::cout << polinom_4(solutions[2], p1, p2, p3, p4) << std::endl;
		std::cout << polinom_4(solutions[3], p1, p2, p3, p4) << std::endl;
*/		    
		    vector<double> temp;
		    
		    double preY=2147483647;
		    for( int i =0;i<4;i++){
		      double res = solutions[i].real();
		      if(abs(solutions[i].imag())<2 && abs(preY-res)>0.1 && (b*b-res*res)>=-0.1 && (r*r-(res-cy)*(res-cy))>=-0.1){
		        temp.push_back(res);
		        preY = res;
		        //cout<<"actually pushed to result: "<<res<<endl;
		      }		    
		    }
		    
		    for(vector<double>::iterator it = temp.begin(); it != temp.end(); ++it){
		      double y = *it;
		      double res1 = sqrt(abs(a*a*(1-y*y/b/b)));
		      double res2 = -sqrt(abs(a*a*(1-y*y/b/b)));
		      double v1 = abs((res1-cx)*(res1-cx)+(y-cy)*(y-cy)-r*r);
		      double v2 = abs((res2-cx)*(res2-cx)+(y-cy)*(y-cy)-r*r);
		      //cout<<"y: "<<y<<"; res1: "<<res1<<"; res2: "<<res2<<endl;
		      //cout<<"v1: "<<v1<<"; v2: "<<v2<<endl;
		      
		      if(abs(res1)<1){
		          res1 = 0;
		          res2 = 0;
		          if(v1<500||v2<500)
		          {
                Point2 tempPoint = {res1, y};
                result.push_back(tempPoint);
		          }
		      }
		      else {
		          if(v1<256){
		            Point2 tempPoint = {res1, y};
		            result.push_back(tempPoint);
		          }
		          if(v2<256){
		            Point2 tempPoint = {res2, y};
		            result.push_back(tempPoint);
		          }
		      }
		      
		    }
		    sortWithAngle(result);

		    delete[] solutions;
		
		}
		return result;
}


vector<tuple<float,float>> getRoot(float x0,float y0,float r,float a,float b,vector<float>& yr)
{
  float aa=a;
  float bb=b;
  float cc=sqrt(a*a+b*b);
  float alph=-2*y0*b*b/(a*a+b*b);
  float beta=(a*a*b*b+x0*x0*b*b+b*b*y0*y0-r*r*b*b)/(a*a+b*b);
  float gama=2*x0*a*b/(a*a+b*b);
  a=1;
  float e=beta*beta-gama*gama*b*b;
  b=2*alph;
  float c=alph*alph+2*beta-gama*gama;
  float d=2*alph*beta;
  
  vector<float> temp = realQudricRoot(b,c,d,e);
  //cout<<"temp size : "<<temp.size()<<endl;
  
  
  vector<tuple<float,float>> res;
  if(!temp.empty()) std::sort(temp.begin(),temp.end());
  float pre=-55555.5;
  float cnt=0;
  for(vector<float>::iterator it=temp.begin();it!=temp.end();it++)
  {
    float y=*it;
    float x=rootX(y,aa,bb);
    
    //cout<<"root y: "<<y<<"; root x: "<<x<<endl;
    //cout<<"cons1: "<<abs(euDis(x,y,x0,y0)-r*r)<< "; cons2: "<<abs(sqrt(euDis(x,y,-cc,0))-sqrt(euDis(x,y,cc,0))-2*aa)<<endl;
    //abs(euDis(x,y,x0,y0)-r*r)<r*150 &&
    if( abs(euDis(x,y,x0,y0)-r*r)<r*r/16 && abs(sqrt(euDis(x,y,-cc,0))-sqrt(euDis(x,y,cc,0))-2*aa)<16 && abs(pre-y) > 5)
    {
      tuple<float,float> curP=make_tuple(x,y);
      cnt+=1;
      res.push_back(curP);
      yr.push_back(y);
      pre=y;
    }
    else
    {
      yr.push_back(abs(euDis(x,y,x0,y0)-r*r));    
    }
    
  }
  return res;
}

float rootX(float y,float a,float b)
{
  return sqrt(b*b+y*y)*a/b; 
}




float hyperArea(float a,float b,tuple<float,float> z1,tuple<float,float> z2)
{
    float z=get<0>(z1)-get<1>(z1)*(get<0>(z2)-get<0>(z1))/(get<1>(z2)-get<1>(z1));
    float m=(get<0>(z2)-get<0>(z1))/(get<1>(z2)-get<1>(z1));
    float y2=get<1>(z2);
    float y1=get<1>(z1);
    return abs(hyperAreaHelper(a,b,m,z,y1)-hyperAreaHelper(a,b,m,z,y2));
}
 
float hyperAreaHelper(float a,float b,float m,float z,float y)
{
    return a/b*(y*sqrt(b*b+y*y)/2 + b*b*log(y+sqrt(b*b+y*y))/2) - m*y*y/2-z*y;
}




vector<float> realQudricRoot(float a, float b, float c, float d)
{
    vector<float> res;
    complex<double>*  solutions = solve_quartic(a, b, c, d);
    for(int i = 0 ; i < 4; i++)
    {
        if(abs(solutions[i].imag())<0.0001)res.push_back(static_cast<float>(solutions[i].real()));
    }
    delete[] solutions;
    return res;
}

#endif
