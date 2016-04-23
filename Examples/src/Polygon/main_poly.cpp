#include "Polygon.hpp"
#include <iostream>
//! Main program
int main()
{
  using namespace Geometry;

  Grid P;
  P.Build("mesh.dat");
  std::cout<<"L'area totale è "<<P.Full_Area()<<"\n\n";
  return 0;

}
  

