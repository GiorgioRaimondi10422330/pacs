#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
struct parameters
{
  //! max number of iteration for Gauss-Siedel
  int   itermax;
  //! Tolerance for stopping criterion
  double  toler;
  //! Bar length
   double L;
  //! First longitudinal dimension
  double a1;
 //! Second longitudinal dimension
  double a2;
  //! Dirichlet condition
  double To;
  //! External temperature 
  double Te;
  //! Conductivity
  double k;
  //! Convection coefficient
  double hc;
  //! Number of elements
  int M;
  //! Se voglio il grafico 0, se voglio i dati in uscita 1 , se voglio entrambi 2
  int Ex;
  //! Norma scelta: 0 ->L2, 1 ->H1
  int norma;
  //! Metodo di risoluzione 0->GaussSiedel  1->Thomas
  int Metodo;
  //! Constructor takes default values
  

  parameters():
    itermax(1000000),
    toler(1e-8),
    L(40.),
    a1(4.),
    a2(50.),
    To(46.),
    Te(20.),
    k(0.164),
    hc(1.e-6*200.),
    M(100),
    Ex(2),
    norma(0),
    Metodo(0)
  {}
};
//! Prints parameters
std::ostream & operator << (std::ostream &,const parameters &);
#endif
