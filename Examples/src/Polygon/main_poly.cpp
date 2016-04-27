#include "Polygon.hpp"
#include <iostream>
//! Main program
int main()
{
  using namespace Geometry;

  std::cout<<"Inizio \n\n";
  Grid G;
  std::cout<<"Grid G\n\n";
  G.Build("mesh.dat");
  std::cout<<"G.Build('mesh.dat')\n\n";
  double P=G.Full_Area();
  std::cout<<"G.Full_Area() "<<P<<" \n\n";

  std::ofstream f("All.dat"),g("External.dat"),h("Internal.dat");

  G.printG(f);
  G.printG(g,1);
  G.printG(h,2);

  return 0;

}
  

