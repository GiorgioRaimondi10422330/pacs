#ifndef ciao_My_Thomas_h
#define ciao_My_Thomas_h


#include<iostream>
#include<vector>

template <typename T>
std::vector<T> My_Thomas(const std::vector<T> &a,const std::vector<T> &b,const std::vector<T> &c, const std::vector<T> &f){
  unsigned const int N=b.size();
  std::vector<T> c_new(N),d(N),x(N);
  c_new[0]=c[0]/b[0];
  d[0]=f[0]/b[0];
  for(unsigned int i=1;i<N-1;i++){
    c_new[i]=c[i]/(b[i]-a[i-1]*c_new[i-1]);
    d[i]=(f[i]-a[i-1]*d[i-1])/(b[i]-a[i-1]*c_new[i-1]);
  }
  x[N-1]=(f[N-1]-a[N-2]*d[N-2])/(b[N-1]-a[N-2]*c_new[N-2]);
  for(unsigned int i=1;i<N;i++){
    x[N-1-i]=d[N-1-i]-c_new[N-1-i]*x[N-i];
  }
  return x;
}

#endif
