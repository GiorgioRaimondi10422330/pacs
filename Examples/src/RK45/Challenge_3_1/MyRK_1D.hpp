#ifndef MYRK_HPP
#define MYRK_HPP
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm> // for max
namespace MyODE
{
  //! Defailt max number of iteration
constexpr std::size_t MAX_STEPS=10000; 

template<class prec=double>
class MyRK1D
 {
 public:
  MyRK1D(): NC(0),NK(0){};

  MyRK1D(
    std::vector<prec> const & c_in, 
    std::vector<prec> const &b1_in, 
    std::vector<prec> const &b2_in, 
    std::vector<std::vector<prec>> const & a_in
    ):c(c_in),b1(b1_in),b2(b2_in),a(a_in)
   {
    NC=c.size();
    NK=a[0].size();
  };

  std::vector<std::pair<prec,prec>> operator ()(
    std::function<prec (prec const &, prec const &)> const & dy, 
    prec const & t0, prec const & T,	 prec const & y0,	 
    prec const & h_initial, 	 
    prec const & h_max, 	 
    prec const & final_error,	 
    int & status,	 
    std::size_t const & maxSteps=MAX_STEPS
  );
 private:
  
  std::vector<prec> c,b1,b2;
  std::vector<std::vector<prec>> a;
  unsigned int NC,NK;
  prec rk_step(
    std::function<prec (prec const &, prec const &)> const & dy,
    prec const & y0,
    prec const & t0,
    prec const & h, 
    prec & error
  );  
};

template<class prec>
  prec MyRK1D<prec>::rk_step(
   std::function<prec (prec const &, prec const &)> const & f,
   prec const & y0, 
   prec const & t0, 
   prec const & h, 
   prec & error
  )
  {
    std::vector<prec> K(NK,0);
    prec sumk=0;
    for(unsigned int i=0;i<NK;i++)
    {
       sumk=0;
       for(unsigned int j=0; j<i;j++)
       {
         sumk+=a[i][j]*K[j]*h;
       }
       K[i]=f(t0+c[i]*h,y0+sumk);
    }

    prec y4 =y0, y5 =y0;
    for (unsigned int i=0;i<NK;i++)
    {
       y4+=h*K[i]*b2[i];
       y5+=h*K[i]*b1[i];
    }
    
    error = std::abs(y5 - y4);
    return y5;
  }

template<class prec>
  std::vector<std::pair<prec,prec>> 
     MyRK1D<prec>::operator ()(
      std::function<prec (prec const &, prec const &)> const & dy,
	    prec const & t0,
	    prec const & T,
	    prec const & y0,
	    prec const & h_initial, 
	    prec const & h_max, 
	    prec const & final_error,
	    int & status,
	    std::size_t const & maxSteps
     )
  {
    status=0;
    const std::size_t maxReduction=maxSteps;
    // parameters for decreasing/increasing time step
    prec const c1=1.0;
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
    prec y(y0);
    prec errorPerTimeStep=final_error/initialNSteps;
    if (initialNSteps>=maxSteps) throw std::runtime_error("RK45: initial time step h too small!");
    std::vector<std::pair<prec,prec>> solution;
    solution.emplace_back(std::make_pair(t0,y0));
    prec localError;
    prec newy;
    while (time<T && stepsCounter <maxSteps)
      {
      	if (time + h > T) h = T-time;
      	newy = rk_step(dy,y,time,h,localError);
      	while (h> h_min && localError > c1*errorPerTimeStep)
      	  {
      	    // half time step
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
          if(stepsCounter>=maxSteps && time < T)
            {
        	status=2;
        	throw std::runtime_error("RK: Max number of time steps exceeded");
            }
    return solution;
  }


}// end namespace



#endif
