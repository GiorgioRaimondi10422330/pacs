#ifndef MYRK_HPP
#define MYRK_HPP
#include <functional>
#include "nonLinSys2.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm> // for max
#include <iostream>
using namespace std;
namespace MyODE
{
  //! Defailt max number of iteration
constexpr std::size_t MAX_STEPS=10000; 
using namespace NonLinearSystems;

template<class prec=double>
class MyRKsys
 {
 public:
  MyRKsys(): NC(0),NK(0){};

  MyRKsys(
    std::vector<prec> const & c_in, 
    std::vector<prec> const &b1_in, 
    std::vector<prec> const &b2_in, 
    std::vector<std::vector<prec>> const & a_in
    ):c(c_in),b1(b1_in),b2(b2_in),a(a_in)
   {
    NC=c.size();
    NK=a[0].size();
  };

  std::vector<std::pair<prec,returnType<prec>>> operator ()(
    NonLinSys<prec> const & dy, 
    prec const & t0, prec const & T,	 argumentType<prec> const & y0,	 
    prec const & h_initial, 	 
    prec const & h_max, 	 
    prec const & final_error,	 
    int & status,
    prec & Tar,
    std::size_t const & maxSteps=MAX_STEPS
  );
 private:
  
  std::vector<prec> c,b1,b2;
  std::vector<std::vector<prec>> a;
  unsigned int NC,NK;
  argumentType<prec> rk_step(
    NonLinSys<prec> const & dy,
    argumentType<prec> const & y0,
    prec const & t0,
    prec const & h, 
    prec & error
  );  
};

template<class prec>
  argumentType<prec> MyRKsys<prec>::rk_step(
   NonLinSys<prec> const & f,
   argumentType<prec> const & y0, 
   prec const & t0, 
   prec const & h, 
   prec & error
  )
  {
    std::size_t N=y0.size();
    std::vector<argumentType<prec>> K(NK,argumentType<prec>(N));
    argumentType<prec> sumk;
    argumentType<prec> IngF(N+1);
    for(unsigned int i=0;i<NK;i++){
       sumk=y0;
       for(unsigned int j=0; j<i;j++){
 	       for(unsigned int q=0; q<N;q++)
         	  sumk[q]+=K[j][q]*(a[i][j]*h);
       }
       for(unsigned int j=0; j<N; j++ ) 
          IngF[j]=sumk[i];
       IngF[N]=t0+c[i]*h;
       K[i]=f(IngF);
    }

    argumentType<prec> y4 =y0, y5 =y0;
    for (unsigned int i=0;i<NK;i++){
       for(unsigned int j=0;j<N; j++){
          y4[j]+=K[i][j]*(b2[i])*h;
          y5[j]+=K[i][j]*(b1[i])*h;
       }
    }
    error=0.;
    for(unsigned int i=0; i<N; i++){
        error+= (y5[i]-y4[i])*(y5[i]-y4[i]);
    }
    error=sqrt(error);
    return y5;
  }

template<class prec>
  std::vector<std::pair<prec,returnType<prec>>> 
     MyRKsys<prec>::operator ()(
      NonLinSys<prec> const & dy,
	    prec const & t0,
	    prec const & T,
	    argumentType<prec> const & y0,
	    prec const & h_initial, 
	    prec const & h_max, 
	    prec const & final_error,
	    int & status,
            prec & Tar,
	    std::size_t const & maxSteps
     )
  {
    status=0;
    const std::size_t maxReduction=maxSteps;
    // parameters for decreasing/increasing time step
    prec const c1=3.0;
    // I need to have a sufficient decrease of the local error
    // to allow time step coarsening
    prec const c2=1./64.;

    prec length=T-t0;
    //! Make sure that h allows to reach T
    std::size_t initialNSteps=std::max(static_cast<size_t>(1),static_cast<size_t>(length/h_initial));
    prec h=length/initialNSteps;
    
    prec h_min = length/(128*maxSteps);
    
    std::size_t stepsCounter(0);
    
    prec time(t0);
    argumentType<prec> y(y0);
    prec errorPerTimeStep=final_error/initialNSteps;
    if (initialNSteps>=maxSteps) throw std::runtime_error("RK45: initial time step h too small!");
    std::vector<std::pair<prec,argumentType<prec>>> solution;
    solution.emplace_back(std::make_pair(t0,y0));
    prec localError;
    argumentType<prec> newy;
    while (time<T && stepsCounter <maxSteps)
      {
      	if (time + h > T) h = T-time;
        unsigned int count=0;
      	newy = rk_step(dy,y,time,h,localError);
      	while (h> h_min && localError > c1*errorPerTimeStep)
      	  {
      	    // half time step
            count++;
      	    h /=2;
      	    errorPerTimeStep /=2;
      	    newy = rk_step(dy,y,time,h,localError);
          }
	      if (localError>errorPerTimeStep)status=1;
      	//! advance
      	y = newy;
      	time +=h;
      	++stepsCounter;
      	solution.emplace_back(std::make_pair(time,y));
      	//! check if we reached end
      	if(localError<c2*errorPerTimeStep && h<h_max)
      	  {
      	    // prec step
      	    h *=2;
      	    errorPerTimeStep *=2;
      	  }
      }
  //handle exceptions
    if(0){
    if(stepsCounter>=maxSteps && time < T)          {
       	status=2;
        throw std::runtime_error("RK: Max number of time steps exceeded");
    }
    }
    return solution;
  }


}// end namespace



#endif
