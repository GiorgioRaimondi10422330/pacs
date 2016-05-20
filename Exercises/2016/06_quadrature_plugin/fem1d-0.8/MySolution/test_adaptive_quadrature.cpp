#include <iostream>
#include <iomanip>
#include <cmath>
#include <dlfcn.h>
#include <functional>
#include "Abstract_Integrator.hpp"
#include "Proxy.hpp"
#include "Factory.hpp"
#include <string>
#include <vector>
#include <functional>
#include "GetPot"
using namespace MyAI;
//double integrand (double x){ return (pow (sin (pow (x, 2)), 2)); }
double integrand (double x){ return (x*x); }


int main (int argc, char **argv)
{
std::cout<<"\n\n\n";
  GetPot cl (argc, argv);
try{
MyFactory & fac=MyFactory::Instance();

if(  cl.search(2,"-h","--help")  )
{
std::vector<std::string> temp=fac.registered();
std::cout<<"[-h][--help] per avere l'help in uscita\n"<<"[-a] double   per fissare l'estremo a dell'intervallo [a,b]\n";
std::cout<<"[-b] double   per fissare l'estremo b dell'intervallo [a,b]\n"<<"[-rule] string per scegliere laregola di integrazione\n";
std::cout<<"Per -rule è possibile scegliere tra:    ";
for(int i=0;i<temp.size();i++){
std::cout<<temp[i]<<"      ";
}
std::cout<<"\n\n";
return 0;
}
double a=cl.follow(0.0,"-a");
double b=cl.follow(10.0,"-b");
std::string rule_choosen=cl.follow("Adaptive","-rule");
bool K=1;
std::vector<std::string> temp=fac.registered();
for(int i=0;i<temp.size();i++){
if(temp[i].compare(rule_choosen)) K=0;
}

std::function<double (double)> f(integrand);

if(K){std::cout<<"La rule "<< rule_choosen<<"  non è tra quelle ammissibili\n eseguire con [-h] o [--help] per le rule ammissibili\n\n"; return 0;}
auto  Metod=fac.create(rule_choosen);

double kk=Metod->integrate(f,0,1);
std::cout<<"Metodo di integrazione scelto  "<<rule_choosen<<"\n";
std::cout<<"l'integrale in I("<<a<<","<<b<<")(f)= "<<kk<<";\n\n";

} catch(...){
std::string rule_choosen=cl.follow("Adaptive","-rule");
std::cout<<"La rule "<< rule_choosen<<"  non è tra quelle ammissibili\n eseguire con [-h] o [--help] per le rule ammissibili\n\n"; return 0;
}

/*
  double (*integrate) (std::function<double (double)>, double, double);
    MyProxy<Adaptive>("Adaptive");

 
  void * handle = dlopen ("adaptive_quadrature.so", RTLD_LAZY);
  if (! handle)
    {
      std::cerr << "cannot load object!" << std::endl;
      std::cerr << dlerror () << std::endl;
      return (-1);
    }

  void * sym = dlsym (handle, "integrate");
  if (! sym)
    {
      std::cerr << "cannot load symbol!" << std::endl;
      std::cerr << dlerror () << std::endl;
      return (1);
    }
  
  integrate = reinterpret_cast<double (*) (std::function<double (double)>, double, double)> (sym);
  
  double pi = 4 * atan (1);
  double res = integrate (integrand, 0, pi);
*/
  return 0;
}
