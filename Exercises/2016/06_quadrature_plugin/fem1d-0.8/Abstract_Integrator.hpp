#ifndef AbInt_hpp
#define ABInt_hpp

#include <functional>
class Abstract_Integrator
{
  public:
   ~integrate(std::function<double (double)>, double, double);
   Abstract_Integrator();
  private:
   ~Abstract_Integrator()=delete;
}





#endif
