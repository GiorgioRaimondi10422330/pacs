#include "midpoint.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "Factory.hpp"
#include "Proxy.hpp"
using namespace MyAI;
double
Midpoint::integrate (std::function<double (double)> f, double a, double b) const
{
  return ((b - a) * f(.5*b + .5*a));
};

 __attribute__((constructor))
static void loadFactoryMidpoint()
{
    MyProxy<Midpoint>("Midpoint");
}


