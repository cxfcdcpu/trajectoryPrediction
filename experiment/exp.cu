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
#include "safeQueue.h"
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
#include "userClass.h"
#include "EllipseTrial.h"

__global__ void goOver3(int n, ellipseTrial *data, float *area,int m){
    int index = threadIdx.x+blockIdx.x*blockDim.x;
    int stride=blockDim.x*gridDim.x;
    for(int k=index;k<n;k+=stride){
        float x = data[k].c3X;
        float y = data[k].c3Y;
        float r = data[k].h3 * data[k].avgD2;
        float a = data[k].ah * data[k].avgD1;
        
        float h1x = data[k].c1X;
        float h1y = data[k].c1Y;
        float h2x = data[k].c2X;
        float h2y = data[k].c2Y;
        
        float rr = r*r;
        float a2 = (data[k].ah-1.0)*data[k].avgD1;
        float total =0.0;
        
        for(int l = 0; l<m;){
            float i = area[l++];
            float j = area[l++];
            float tt = sqrtf((i-h1x)*(i-h1x)+(j-h1y)*(j-h1y))+sqrtf((i-h2x)*(i-h2x)+(j-h2y)*(j-h2y));
            float di = x-i;
            float dj = y-j;
            if(di * di + dj * dj <= rr && tt <= 2*a && tt >= 2* a2) total+=1.0;
        }

        float rate3 = data[k].rate3;
        data[k].grAr=rate3*total;
        data[k].acAr = total;

    }
}

__global__ void bestEllipse(int n, ellipseTrial *data, float *area,int m){
    int index = threadIdx.x+blockIdx.x*blockDim.x;
    int stride=blockDim.x*gridDim.x;
    for(int k=index;k<n;k+=stride){
        float x = data[k].c3X;
        float y = data[k].c3Y;
        float r = data[k].h3 * data[k].avgD2;
        float a = data[k].ah * data[k].avgD1;
        
        float h1x = data[k].c1X;
        float h1y = data[k].c1Y;
        float h2x = data[k].c2X;
        float h2y = data[k].c2Y;
        
        float rr = r*r;
        float a2 = (data[k].ah-1.0)*data[k].avgD1;
        float total =0.0;
        
        for(int l = 0; l<m;){
            float i = area[l++];
            float j = area[l++];
            float tt = sqrtf((i-h1x)*(i-h1x)+(j-h1y)*(j-h1y))+sqrtf((i-h2x)*(i-h2x)+(j-h2y)*(j-h2y));
            float di = x-i;
            float dj = y-j;
            if(di * di + dj * dj <= rr && tt <= 2*a && tt >= 2* a2) total+=1.0;
        }
        float rate = total / data[k].tArea;
        float rate3 = rate*rate*rate;
        if ( rate3<1.1){
            data[k].grAr=rate3*total;
            data[k].acAr = total;
            data[k].rate3 = rate3;
        }
    }
}




__global__ void goOver2(int n, hyperTrial *data, float *area,int m){
    int index = threadIdx.x+blockIdx.x*blockDim.x;
    int stride=blockDim.x*gridDim.x;
    for(int k=index;k<n;k+=stride){
        float x = data[k].c3X;
        float y = data[k].c3Y;
        float r = data[k].h3 * data[k].avgD2;
        float a = data[k].ah * data[k].avgD1;
        
        float h1x = data[k].c1X;
        float h1y = data[k].c1Y;
        float h2x = data[k].c2X;
        float h2y = data[k].c2Y;
        
        float rr = r*r;
        float a2 = (data[k].ah-1.0)*data[k].avgD1;
        float total =0.0;
        
        for(int l = 0; l<m;){
            float i = area[l++];
            float j = area[l++];
            float tt = sqrtf((i-h1x)*(i-h1x)+(j-h1y)*(j-h1y))-sqrtf((i-h2x)*(i-h2x)+(j-h2y)*(j-h2y));
            float di = x-i;
            float dj = y-j;
            if(di * di + dj * dj <= rr && tt <= 2*a && tt >= 2* a2) total+=1.0;
        }

        float rate3 = data[k].rate3;
        data[k].grAr=rate3*total;
        data[k].acAr = total;

    }
}

__global__ void bestHyper(int n, hyperTrial *data, float *area,int m){
    int index = threadIdx.x+blockIdx.x*blockDim.x;
    int stride=blockDim.x*gridDim.x;
    for(int k=index;k<n;k+=stride){
        float x = data[k].c3X;
        float y = data[k].c3Y;
        float r = data[k].h3 * data[k].avgD2;
        float a = data[k].ah * data[k].avgD1;
        
        float h1x = data[k].c1X;
        float h1y = data[k].c1Y;
        float h2x = data[k].c2X;
        float h2y = data[k].c2Y;
        
        float rr = r*r;
        float a2 = (data[k].ah-1.0)*data[k].avgD1;
        float total =0.0;
        
        for(int l = 0; l<m;){
            float i = area[l++];
            float j = area[l++];
            float tt = sqrtf((i-h1x)*(i-h1x)+(j-h1y)*(j-h1y))-sqrtf((i-h2x)*(i-h2x)+(j-h2y)*(j-h2y));
            float di = x-i;
            float dj = y-j;
            if(di * di + dj * dj <= rr && tt <= 2*a && tt >= 2* a2) total+=1.0;
        }
        float rate = total / (data[k].tArea +0.1);
        float rate3 = rate*rate*rate;
        if ( rate3<1.1){
            data[k].grAr=rate3*total;
            data[k].acAr = total;
            data[k].rate3 = rate3;
        }
    }
}

__global__ void bestTwoCycle(int n, twoCycleTrial *data, float *area,int m){
    int index = threadIdx.x+blockIdx.x*blockDim.x;
    int stride=blockDim.x*gridDim.x;
    for(int k=index;k<n;k+=stride){
        float x1=data[k].c1X;
        float y1=data[k].c1Y;
        float x2=data[k].c2X;
        float y2=data[k].c2Y;
        float r1=data[k].h1*data[k].d;
        float r2=data[k].h2*data[k].d;
        float r3=r1-data[k].d;
        float rr3=r3*r3;
        float rr1=r1*r1;
        float rr2=r2*r2;
        float total=0.0;
        for(int l=0;l<m;){
            float i=area[l++];
            float j=area[l++];
            float di1=x1-i;
            float dj1=y1-j;
            float di2=x2-i;
            float dj2=y2-j;
            if (di1*di1+dj1*dj1<=rr1 && di1*di1+dj1*dj1>rr3 && di2*di2+dj2*dj2<=rr2)
                total+=1.0;
        }
        float rate=total / (data[k].tArea+0.1);
        float rate3=rate*rate*rate;
        if (rate3<1.1){
                
            data[k].grAr=rate3*total ; 
            data[k].acAr=total;
            data[k].rate3=rate3;
        }
    }    

}

__global__ void goOver1(int n, twoCycleTrial *data, float *area,int m){
    int index = threadIdx.x+blockIdx.x*blockDim.x;
    int stride=blockDim.x*gridDim.x;
    for(int k=index;k<n;k+=stride){
        float x1=data[k].c1X;
        float y1=data[k].c1Y;
        float x2=data[k].c2X;
        float y2=data[k].c2Y;
        float r1=data[k].h1*data[k].d;
        float r2=data[k].h2*data[k].d;
        float r3=r1-data[k].d;
        float rr3=r3*r3;
        float rr1=r1*r1;
        float rr2=r2*r2;
        float total=0.0;
        for(int l=0;l<m;){
            float i=area[l++];
            float j=area[l++];
            float di1=x1-i;
            float dj1=y1-j;
            float di2=x2-i;
            float dj2=y2-j;
            if (di1*di1+dj1*dj1<=rr1 && di1*di1+dj1*dj1>rr3 && di2*di2+dj2*dj2<=rr2)
                total+=1.0;
        }

        float rate3=data[k].rate3;

        data[k].grAr=rate3*total ; 
        data[k].acAr=total;

    }    

}

string getRoutingMSG_without_ellipse(User & curUser,vector<twoCycleTrial> & cycleTrials,vector<hyperTrial> & hyperTrials)
{
  string res = "";
  float *TAS,*d_TAS;
  int tasSize=curUser.TAS.size();
  TAS=(float*)malloc(sizeof(float)*tasSize*2);
  cudaMalloc((void**)&d_TAS, sizeof(float) *tasSize*2);
  int counter=0;
  for(string t:curUser.TAS)
  {
      stringstream tt(t);
      float x,y;
      tt>>x;
      tt>>y;
      TAS[counter++]=x;
      TAS[counter++]=y;
  }
  cudaMemcpy(d_TAS, TAS, sizeof(float) *tasSize*2, cudaMemcpyHostToDevice);
  
  //cout<<"I'm OK Here"<<endl;
  counter=0; 
  twoCycleTrial *cTri;
  twoCycleTrial *d_cTri;
  if(!cycleTrials.empty())
  {
      cTri=(twoCycleTrial*)malloc(sizeof(twoCycleTrial)*cycleTrials.size());
      
      //cout<<"I'm OK after here"<<endl;
      for(twoCycleTrial ct: cycleTrials)
      {
          cTri[counter++]=ct;
      }
      cudaMalloc((void**)&d_cTri, sizeof(twoCycleTrial) *cycleTrials.size());
      cudaMemcpy(d_cTri,cTri,sizeof(twoCycleTrial) *cycleTrials.size(),cudaMemcpyHostToDevice);
      cout<<"finish copy totoal trial: "<<cycleTrials.size()<<endl;        
      bestTwoCycle<<<2048,256>>>(cycleTrials.size(),d_cTri,d_TAS,tasSize*2);
     
      //cudaFree(d_cTri);
  } 
  //cout<<"number of hyperTrial: "<<hyperTrials.size()<<endl;
  hyperTrial *hTri;
  hyperTrial *d_hTri;
  int counter2=0;
  if(!hyperTrials.empty())
  {
      hTri=(hyperTrial*)malloc(sizeof(hyperTrial)*hyperTrials.size());

      //cout<<"I'm OK after here"<<endl;
      for(hyperTrial ht: hyperTrials)
      {
          hTri[counter2++]=ht;
      }
      cudaMalloc((void**)&d_hTri, sizeof(hyperTrial) *hyperTrials.size());        
      cudaMemcpy(d_hTri,hTri,sizeof(hyperTrial) *hyperTrials.size(),cudaMemcpyHostToDevice);
      cout<<"finish copy totoal trial: "<<hyperTrials.size()<<endl;
      cout<<"********TAS size is: "<<tasSize<<"******"<<endl;
      bestHyper<<<2048,256>>>(hyperTrials.size(),d_hTri,d_TAS,tasSize*2);
      //cudaDeviceSynchronize();
      
      //cudaFree(d_hTri);
  }
  
  if(cycleTrials.empty() || hyperTrials.empty())return "No result";
  cudaDeviceSynchronize();
  cudaMemcpy(cTri,d_cTri,sizeof(twoCycleTrial) *cycleTrials.size(),cudaMemcpyDeviceToHost); 
  cudaMemcpy(hTri,d_hTri,sizeof(hyperTrial) *hyperTrials.size(),cudaMemcpyDeviceToHost);  
    
  cout<<counter<<" ||||||| "<<counter2<<" ||||||||||  "<<tasSize<<endl;
  res+=findBestTry(cTri,hTri, counter, counter2, curUser.TAS);
  
  sort(cTri, cTri+counter, sortCycleTrial);
  sort(hTri, hTri+counter2, sortHyperTrial);
  counter = 1000000 < counter ? 1000000 : counter;
  counter2 = 5000000 < counter2 ? 5000000 : counter2;
  
  
  //cudaMalloc((void**)&d_cTri, sizeof(twoCycleTrial) *counter);
  cudaMemcpy(d_cTri,cTri,sizeof(twoCycleTrial) *counter, cudaMemcpyHostToDevice);
  //cudaMalloc((void**)&d_hTri, sizeof(hyperTrial) *counter2);  
  cudaMemcpy(d_hTri,hTri,sizeof(hyperTrial) *counter2, cudaMemcpyHostToDevice);

  int newSize = 0;
  int ccc = 0;
  do
  {
    ccc++;
    newSize=curUser.TAS.size();
    cout<<counter<<" ||||||| "<<counter2<<" ||||||||||  "<<newSize<<endl;
    int tasInd=0;
    for(string t:curUser.TAS)
    {
        stringstream tt(t);
        float x,y;
        tt>>x;
        tt>>y;
        TAS[tasInd++]=x;
        TAS[tasInd++]=y;
    }
    //cudaFree(d_TAS);
    //cudaMalloc((void**)&d_TAS, sizeof(float) *newSize*2);
    cudaMemcpy(d_TAS, TAS, sizeof(float) *newSize*2, cudaMemcpyHostToDevice);
    goOver1<<<2048,256>>>(counter, d_cTri, d_TAS, newSize*2);
    goOver2<<<2048,256>>>(counter2, d_hTri, d_TAS, newSize*2);
    cudaDeviceSynchronize();
    
    //free(cTri);
    //free(hTri);
    //cTri=(twoCycleTrial*)malloc(sizeof(twoCycleTrial)*counter);
    //hTri=(hyperTrial*)malloc(sizeof(hyperTrial)*counter2);
    cudaMemcpy(cTri,d_cTri,sizeof(twoCycleTrial) *counter,cudaMemcpyDeviceToHost);  
    cudaMemcpy(hTri,d_hTri,sizeof(hyperTrial) *counter2,cudaMemcpyDeviceToHost);
    res+=findBestTry(cTri,hTri, counter, counter2, curUser.TAS);
  }while(newSize>0.1*tasSize && ccc<30);
  
  cudaFree(d_cTri);
  cudaFree(d_hTri);
  cudaFree(d_TAS);
  free(TAS);
  free(cTri); 
  free(hTri);     
  return res;
  
}

string getRoutingMSG_with_ellipse(User & curUser,vector<ellipseTrial> ellipseTrials,vector<hyperTrial> & hyperTrials)
{
  string res="";
  float *TAS,*d_TAS;
  int tasSize=curUser.TAS.size();
  TAS=(float*)malloc(sizeof(float)*tasSize*2);
  cudaMalloc((void**)&d_TAS, sizeof(float) *tasSize*2);
  int counter=0;
  for(string t:curUser.TAS)
  {
      stringstream tt(t);
      float x,y;
      tt>>x;
      tt>>y;
      TAS[counter++]=x;
      TAS[counter++]=y;
  }
  cudaMemcpy(d_TAS, TAS, sizeof(float) *tasSize*2, cudaMemcpyHostToDevice);
  
  counter=0; 
  ellipseTrial *eTri;
  ellipseTrial *d_eTri;
  if(!ellipseTrials.empty())
  {
      eTri=(ellipseTrial*)malloc(sizeof(ellipseTrial)*ellipseTrials.size());
      
      //cout<<"I'm OK after here"<<endl;
      for(ellipseTrial et: ellipseTrials)
      {
          eTri[counter++]=et;
      }
      cudaMalloc((void**)&d_eTri, sizeof(ellipseTrial) * ellipseTrials.size());
      cudaMemcpy(d_eTri,eTri,sizeof(ellipseTrial) * ellipseTrials.size(),cudaMemcpyHostToDevice);
      cout<<"finish copy totoal trial: "<<ellipseTrials.size()<<endl;        
      bestEllipse<<<2048,256>>>(ellipseTrials.size(),d_eTri,d_TAS,tasSize*2);
     
      
      
  }
  cudaDeviceSynchronize();
  cudaMemcpy(eTri,d_eTri,sizeof(ellipseTrial) *ellipseTrials.size(),cudaMemcpyDeviceToHost); 
  cudaDeviceSynchronize();
  cudaFree(d_eTri);
  
  
  hyperTrial *hTri;
  hyperTrial *d_hTri;
  int counter2=0;
  
  if(!hyperTrials.empty())
  {
      hTri=(hyperTrial*)malloc(sizeof(hyperTrial)*hyperTrials.size());

      //cout<<"I'm OK after here"<<endl;
      for(hyperTrial ht: hyperTrials)
      {
          hTri[counter2++]=ht;
      }
      cudaMalloc((void**)&d_hTri, sizeof(hyperTrial) *hyperTrials.size());        
      cudaMemcpy(d_hTri,hTri,sizeof(hyperTrial) *hyperTrials.size(),cudaMemcpyHostToDevice);
      cout<<"finish copy totoal trial: "<<hyperTrials.size()<<endl;
      cout<<"********TAS size is: "<<tasSize<<"******"<<endl;
      bestHyper<<<2048,256>>>(hyperTrials.size(),d_hTri,d_TAS,tasSize*2);
      //cudaDeviceSynchronize();
      
      //cudaFree(d_hTri);
  }
  
  if(ellipseTrials.empty() || hyperTrials.empty())return "No result";
  cudaDeviceSynchronize();
  cudaMemcpy(hTri,d_hTri,sizeof(hyperTrial) *hyperTrials.size(),cudaMemcpyDeviceToHost);  
  cudaDeviceSynchronize();
  cudaFree(d_hTri);
  
    
  cout<<counter<<" ||||||| "<<counter2<<" ||||||||||  "<<tasSize<<endl;
  res+=findBestTry2(eTri,hTri, counter, counter2, curUser.TAS);
  
  sort(eTri, eTri+counter, sortEllipseTrial);
  sort(hTri, hTri+counter2, sortHyperTrial);
  counter = 5000000 < counter ? 5000000 : counter;
  counter2 = 5000000 < counter2 ? 5000000 : counter2;
  
  
  cudaMalloc((void**)&d_eTri, sizeof(ellipseTrial) *counter);
  cudaMemcpy(d_eTri,eTri,sizeof(ellipseTrial) *counter, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_hTri, sizeof(hyperTrial) *counter2);  
  cudaMemcpy(d_hTri,hTri,sizeof(hyperTrial) *counter2, cudaMemcpyHostToDevice);
  cout<<"copy to device again"<<endl;
  int newSize = curUser.TAS.size();
  cout<<"new TAS SIZE is: "<<newSize<<"; old TAS size is"<<tasSize<<endl;
  cudaDeviceSynchronize();
  while(newSize>0.10*tasSize)
  {
    
    newSize=curUser.TAS.size();
    cout<<counter<<" ||||||| "<<counter2<<" ||||||||||  "<<newSize<<endl;
    int tasInd=0;
    for(string t:curUser.TAS)
    {
        stringstream tt(t);
        float x,y;
        tt>>x;
        tt>>y;
        TAS[tasInd++]=x;
        TAS[tasInd++]=y;
    }
    cudaFree(d_TAS);
    cudaMalloc((void**)&d_TAS, sizeof(float) *newSize*2);
    cudaMemcpy(d_TAS, TAS, sizeof(float) *newSize*2, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    goOver3<<<2048,256>>>(counter, d_eTri, d_TAS, newSize*2);
    cudaDeviceSynchronize();
     
    
    
    
    
    cudaDeviceSynchronize();
    goOver2<<<2048,256>>>(counter2, d_hTri, d_TAS, newSize*2);
    cudaDeviceSynchronize();
    cudaMemcpy(eTri,d_eTri,sizeof(ellipseTrial) *counter,cudaMemcpyDeviceToHost); 
    cudaMemcpy(hTri,d_hTri,sizeof(hyperTrial) *counter2,cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    
    //free(cTri);
    //free(hTri);
    //cTri=(twoCycleTrial*)malloc(sizeof(twoCycleTrial)*counter);
    //hTri=(hyperTrial*)malloc(sizeof(hyperTrial)*counter2);
    
    
    res+=findBestTry2(eTri,hTri, counter, counter2, curUser.TAS);
  }
  cudaDeviceSynchronize();
  cudaFree(d_eTri);
  cudaFree(d_hTri);
  cudaFree(d_TAS);
  free(TAS);
  free(eTri); 
  free(hTri);
  
  //cout<<"sending back result!!!!!!!!!!!!!!!"<<endl;
    
  return res;


}

void writeToFile(string msg, int sn, int id, int mode,string network)
{

  string filename =  "../results/"+network+"_d/"+ to_string(sn)+"_"+to_string(id)+"_";
  if(mode==0)filename+="withoutEllipse";
  else filename+="withEllipse";
  ofstream resFile;
  resFile.open(filename);
  resFile<<msg;  
  resFile.close();
}

int main(int argc , char *argv[])
{  
  srand(time(NULL));
  User curUser;
  int mode;
  while(1){
    cout << "type in the mode: "<<endl;
    cout << "1: create random network"<<endl;
    cout << "2: find route msg for specific random network"<<endl;
    
    
    cin >> mode;
    if(mode ==1){
      int nn, an , rr;
      cout<< "Please type in number of nodes"<<endl;
      cin >> nn;
      cout<< "please type in maximum number of anchor nodes"<<endl;
      cin >> an;
      cout<< "please type in the radio range"<<endl;
      cin >> rr;
      string nid = curUser.randomNetwork(nn,an,rr);
      cout<< "created random network: "<<nid<<endl;
    }
    else if(mode == 2){
      string network;
      cout<< "Please type in the network ID"<<endl;
      cin >> network;
      curUser.loadNetwork(network);
      short hv[nodeSize][anchorSize];
      curUser.getHopInfo(hv);
      
      vector<twoCycleTrial> cycleTrials=curUser.findTwoCycleTrial(hv);
      cout<<"number of cycleTrials: "<<cycleTrials.size()<<endl;
      vector<ellipseTrial> ellipseTrials=curUser.findEllipseTrial(hv);
      cout<<"number of ellipseTrials: "<<ellipseTrials.size()<<endl;
      
      vector<hyperTrial> hyperTrials=curUser.findHyperTrial(hv);
      cout<<"number of hyperTrial: "<<hyperTrials.size()<<endl; 
      for(int i=180; i>=20; i--)
      {
        for(int j=0; j<20; j++)
        {
          curUser.updateStroke(i,j,50);
          cout<<"curUser TAS: "<<curUser.TAS.size()<<endl;
          string res = getRoutingMSG_without_ellipse(curUser,cycleTrials,hyperTrials);
          writeToFile(res,i,j,0,network);
          cout<<i<<"_"<<j<<" result without ellipse: "<<res<<endl;
          
          //string res2 = getRoutingMSG_with_ellipse(curUser,ellipseTrials,hyperTrials);
          //writeToFile(res2,i,j,1,network);
          //cout<<i<<"_"<<j<<" result with ellipse: "<< res2<<endl;
        }
      }
    
    }
    else{
      cout<<" please type 1 or 2"<<endl;
    
    }
  }


	return 0;
}



