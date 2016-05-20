#ifndef AbInt_hpp
#define ABInt_hpp

#include <functional>
#include "Factory.hpp"
#include "Proxy.hpp"
namespace MyAI{
class Abstract_Integrator
{
  public:
   Abstract_Integrator()=default;
   Abstract_Integrator & operator =(Abstract_Integrator &)=default;
   virtual double integrate (std::function<double (double)> f, double a, double b) const=0;
   
};


using MyFactory=GenericFactory::Factory<Abstract_Integrator,std::string>;

template<class C>
using MyProxy=GenericFactory::Proxy<MyFactory,C>;
}

#endif
