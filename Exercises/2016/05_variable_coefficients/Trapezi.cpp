#include<iostream>
#include "muParserInterface.hpp"
#include<string>
#include"mesh.h"

using namespace std;
using namespace MuParserInterface;
auto Trapezi(const muParserInterface &func,const muParserInterface &phi,const double a,const duoble b)
{
	array<double,2> A={{a,b}},B(0)={{b,b}};
	return (func(a,A.data())*phi(a,A.data())+func(a,B.data())*phi(a,B.data()))*(b-a)/2;
}
