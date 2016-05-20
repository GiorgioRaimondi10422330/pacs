
#ifndef HAVE_MIDPOINT_H
#define HAVE_MIDPOINT_H
#include <functional>
#include "Abstract_Integrator.hpp"
using namespace MyAI;
class Midpoint final:public Abstract_Integrator
{
  public:  
  double integrate (std::function<double (double)> f, double a, double b) const;
};

#endif
