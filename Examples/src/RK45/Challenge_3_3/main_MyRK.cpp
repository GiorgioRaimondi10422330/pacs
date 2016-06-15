#include "MyRK_Sys.hpp"
#include "method.hpp"
#include <iostream>
#include <functional>
#include <fstream>
#include <cmath>
#include <vector>
#include "nonLinSys2.hpp"


int main()
{
  using namespace std;
  using namespace MyODE;
  using namespace NonLinearSystems;
  auto fun1 = [](argumentType<double> const & X){
   double t=X[2],x=X[0],y=X[1];
   return y;
  };
  auto fun2 = [](argumentType<double> const & X){
    double t=X[2],x=X[0],y=X[1];
    return 2;
  };
  int scale=10;
  NonLinSys<double> f;
  f.addToSystem(fun1);
  f.addToSystem(fun2);
  auto exac1 = [](double const & t){return t*t;};
  auto exac2 = [](double const & t){return 2*t;};
  double t0=0;
  argumentType<double> y0(2,0.);
  double T=10;
  double h_init=0.2;
  double errorDesired=1.e-2;
  double Tarrivo=0.;
  double LocalErr(0.);
  double tt=T/h_init;
  int status;
  vector<double> c, b1, b2;
  vector<vector<double>>  a;
  int ingr;
  cout<<"Se vuoi utilizzare Runge-Kutta 45 digita 1\n"<<"Se vuoi utilizzare Runge-Kutta 12 digita 2\n";
  cin>>ingr;
  if(ingr==1){
    RK45(a,b1,b2,c); // imposta il butcher per RK45
    cout<<"Metodo scelto RK45\n";
  }
  else if(ingr==2){
    RK12(a,b1,b2,c);
    h_init=h_init*4;
    errorDesired=0.1;
    cout<<"Metodo scelto RK12\n";
  }
  else{
    cout<<"Non è stato scelto nessun metodo tra quelli proposti\n\n";
    return 0;
  }

 MyRKsys<double> rk(c,b1,b2,a);
 try{
   auto result= rk(f,t0,T,y0,h_init/scale,(T-t0)/4.,errorDesired,status,LocalErr,Tarrivo,10000*scale);
 
   double error1=0.;
   double error2=0.;
   cout<<"t0="<<t0<<"  T="<<T<<"  h="<<h_init<<"  h_max="<<(T-t0)/4.<<"  error="<<errorDesired<<"\n\n";
   double error(0.);
   for (std::size_t i=0; i<result.size();i++) 
   {
     error1= abs(result[i].second[0]- exac1(result[i].first));
     error2= abs(result[i].second[0]- exac2(result[i].first));
     error=max(error, sqrt(error1*error1+error2*error2));
   }
   ofstream file("result.tex");
   for (auto v : result){
      file<<v.first;
      for(auto i: v.second) file<<" "<<i;
      file<<std::endl;
   }
   file.close();
 
   cout<<"errore voluto "<<errorDesired<<" the error is "<<error<<"\n"; 
   cout<<"Tempo finale="<<T<<" Tempo raggiunto con la discretizzazione="<<Tarrivo<<"\n\n";

  } catch(...){
    cout<<"Il metodo scelto non converge all'errore perchè richiede una\ndiscretizzazione maggiore a quella permessa\n\n";
  }
  
  return 0;
    
}
