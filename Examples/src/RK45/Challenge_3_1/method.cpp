#include <vector>
#include "method.hpp"

void RK45(std::vector<std::vector<double>> & a, std::vector<double> & b1, std::vector<double> & b2, std::vector<double> & c)
{
a=std::vector<std::vector<double>> (6,std::vector<double> (6,0.));
b1=std::vector<double> (6,0.);
b2=std::vector<double> (6,0.);
c=std::vector<double> (6,0.);
  
  c[0]=0.;
  c[1]=1./4.;
  c[2]=9./32.+3./32.;
  c[3]=1932./2197.-7200./2197.+7296./2197.;
  c[4]=439./216.- 8.+3680./513.-845./4104.;
  c[5]=-8./27.+2.-3544./2565.+1859./4104.-11./40.;
  b1[0]=16./135.;
  b1[1]=0.;
  b1[2]=6656./12825.;
  b1[3]=28561./56430.;
  b1[4]=-9./50.;
  b1[5]= 2./55.;
  b2[0]=25./216.;
  b2[1]=0.;
  b2[2]=1408./2565.;
  b2[3]=2197./4104.;
  b2[4]=-1./5.;b2[5]=0.;
  a[1][0]=1./4.;
  a[2][0]=3./32.;a[2][1]=9./32.;
  a[3][0]=1932./2197.;a[3][1]=-7200./2197.;a[3][2]=7296./2197.;
  a[4][0]=439./216.;a[4][1]=- 8.;a[4][2]=3680./513.;a[4][3]=-845./4104.;
  a[5][0]=-8./27.;a[5][1]=2.;a[5][2]=-3544./2565.;a[5][3]=1859./4104.;a[5][4]=-11./40.;

}


void RK12(std::vector<std::vector<double>> & a, std::vector<double> & b1, std::vector<double> & b2, std::vector<double> & c)
{
  a=std::vector<std::vector<double>> (2,std::vector<double> (2,0.));
  b1=std::vector<double> (2,0.);
  b2=std::vector<double> (2,0.);
  c=std::vector<double> (2,0.);
  a[1][0]=1.;
  c[1]=1.;
  b2[0]=1.;
  b1[0]=1./2.;
  b1[1]=1./2.;

}



