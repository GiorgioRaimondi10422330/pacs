#include "MyRK_1D.hpp"
#include "method.hpp"
#include <iostream>
#include <functional>
#include <fstream>
#include <cmath>
#include <vector>


int main()
{
  using namespace std;
  using namespace MyODE;
  //auto fun = [](double const & t, double const & y){return -10*y;};
  auto fun = [](double const & t, double const & y){return t;};
  auto exac  = [](double const & t){return t*t/2.+1;};
  double t0=0;
  double y0=1;
  double T=100;
  double h_init=0.2;
  double errorDesired=1.e-7;
  int status;
  vector<double> c, b1, b2;
  vector<vector<double>>  a;
  int ingr;
  cout<<"Se vuoi utilizzare Runge-Kutta 45 digita 1\n"<<"Se vuoi utilizzare Runge-Kutta 12 digita 2\n";
  cin>>ingr;
if(ingr==1)
  RK45(a,b1,b2,c); // imposta il butcher per RK45
else if(ingr==2){
  RK12(a,b1,b2,c);
  h_init=h_init*4;
  errorDesired=0.1;
}
else
{
 cout<<"Non è stato scelto nessun metodo tra quelli proposti\n\n";
 return 0;
}


 MyRK1D<double> rk(c,b1,b2,a);
try{
 auto result= rk(fun,t0,T,y0,h_init,(T-t0)/4.,errorDesired,status,10000);
 
 double error=0.;
 for (std::size_t i=0; i<result.size();i++) 
 {
  error=max(error, abs(result[i].second- exac(result[i].first)));
 }
  ofstream file("result.tex");
  for (auto v : result)
    file<<v.first<<" "<<v.second<<std::endl;
  file.close();
 cout<<"errore voluto "<<errorDesired<<"      errore ottenuto "<<error<<"\n\n";
} catch(...){
cout<<"Il metodo scelto non converge all'errore perchè richiede una\ndiscretizzazione maggiore a quella permessa\n\n";
}
  return 0;
    
}
