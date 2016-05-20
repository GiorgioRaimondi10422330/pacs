
#ifndef HAVE_TRAPEZOIDAL_H
#define HAVE_TRAPEZOIDAL_H
#include <functional>
#include "Abstract_Integrator.hpp"
using namespace MyAI;
class Trapezoidal final:public Abstract_Integrator
{
  public:
  double integrate (std::function<double (double)> f, double a, double b) const;
};

#endif
