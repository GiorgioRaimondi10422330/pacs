
#ifndef HAVE_ADAPTIVE_QUADRATURE_H
#define HAVE_ADAPTIVE_QUADRATURE_H
#include <functional>
#include "Abstract_Integrator.hpp"



class Adaptive final: public MyAI::Abstract_Integrator
{
  public:
  Adaptive()=default;
  double integrate (std::function<double (double)> f, double a, double b) const;
};

#endif
