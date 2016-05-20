#include "fem1d.h"
#include "matrix.h"
#include "GetPot"
#include "muParserInterface.hpp"
#include<string>
#include <vector>

using namespace std;

int main (int argc, char **argv)
{

  GetPot cl (argc, argv);

  if (cl.search (2, "-h", "--help"))
    {
      std::cerr << help_text << std::endl;
      return 0;
    }
  const double a = cl.follow (0.0, "-a");
  const double b = cl.follow (1.0, "-b");
  const unsigned int nnodes = cl.follow (100, 2, "-n", "-nnodes");

  mesh m (a, b, nnodes);
  matrix A(nnodes);
  muParserInterface diffusion, forzante;
  vector<muParserInterface> phi(2), phigrad(2);
  string p0="(t-x)/(t-y)", p1="(x-y)/(t-y)",pg0="1/(t-y)*1/(t-y)",pg1="-1/(t-y)*1/(t-y)";
  phi[0].set_expression(p1);
  phi[1].set_expression(p0);
  phigrad[0].set_expression(pg0);
  phigrad[1].set_expression(pg1);


  matrix mloc(2);
  for (unsigned int iel = 0; iel < m.nels; ++iel)
    {

      std::fill (mloc.get_data (),
                 mloc.get_data () + 4,
                 0.0);
      
      for (unsigned int inode = 0; inode < 2; ++inode)
        {
          for (unsigned int jnode = 0; jnode < 2; ++jnode)
            {
              if(inode!=jnode)
		{
                  mloc(inode,jnode) = Trapezi(duffusion,phigrad[1]);
		}
	      else
		{
		  mloc(inode,jnode) = Trapezi(duffusion,phigrad[0]);
		}
              A(m.elements[iel][inode],m.elements[iel][jnode]) += mloc(inode,jnode);
            }
        }
    }

  matrix f(nnodes, 1);
  matrix vloc(2, 1);
  
  for (unsigned int iel = 0; iel < m.nels; ++iel)
    {
      std::fill (mloc.get_data (),
                 mloc.get_data () + 2,
                 0.0);
      
      for (unsigned int inode = 0; inode < 2; ++inode)
        {
          vloc(inode, 0) = Trapezi(forzante,phi[inode]);
          f(m.elements[iel][inode], 0) += vloc(inode, 0);
        }
    }

  f(0, 0) = 0;
  f(nnodes - 1, 0) = 0;

  A(0,0) = 1.0;
  A(nnodes-1,nnodes-1) = 1.0;
  for (unsigned int ii = 1; ii < nnodes; ++ii)
    {
      A(0, ii) = 0.0;
      A(nnodes-1, nnodes-1-ii) = 0.0;
    }

  
  matrix uh(f);
  A.solve (uh);

  for (unsigned int ii = 0; ii < nnodes; ++ii)
    std::cout << uh(ii, 0) << std:: endl;
      
  return 0;
};



