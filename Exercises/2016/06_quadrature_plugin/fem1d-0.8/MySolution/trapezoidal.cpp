#include "trapezoidal.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "Factory.hpp"
#include "Proxy.hpp"
using namespace MyAI;

double
Trapezoidal::integrate (std::function<double (double)> f, double a, double b) const
{
  return ((b - a) * (.5 * f(b) + .5 * f(a)));
};


 __attribute__((constructor))
static void loadFactoryTrapezoidal()
{
    MyProxy<Trapezoidal>("Trapezoidal");
}


