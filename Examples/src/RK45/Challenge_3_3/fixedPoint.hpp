#ifndef HH_FIXEDPOINT_HH_
#define HH_FIXEDPOINT_HH_
#include  "nonLinSys2.hpp"
#include <cmath>
namespace NonLinearSystems{
  //! Struct that holds the option for Fixed Point iterations.
  /*!
    It provides the parameter to guide the iterations.
   */
  template<class T=double>
  struct FixedPointOptions{
    //! @brief Tolerance on residual.
    /*! @detail The iteration stops if \f$||F(x)||<minRes\f$.*/
    T minRes;
    //! MAx. number of iterations.
    unsigned int maxIter;
    //! Relaxation parameter
    T alpha;
  };

  //! Status returned by the Fixed Point Iteration.
  template<class T=double>
  struct FixedPointStatus{
    //! Last residual.
    T residual;
    //! Last iteration.
    unsigned int iterations;
    //! @brief Converged flag.
    /*! @detail true if two successive iterations smaller
     than given tolerance*/
    bool converged;
  };
  
  //! Function implementing a fixed point method for zero of non-linear systems
  /*!

    It solves
    \f[
    x^{k+1}=x^{k}-\alpha r(x)
    \f]
    where \f$ r(x) \f$ is the residual;
 
    The iteration stops under one of the three conditions:
    - number of iterations exceeded
    - residual less or equal opt.minRes
    - norm of two successive iteration below opt.tolerance
    
    @param fSys The non linear system.
    @param x    In input the starting point. In output the found point.
    @param opt  The given options. 
  */

  template<class T> 
  T norm(argumentType<T> const & X){
    T norma=0.;
    for(auto &i : X) norma+=i*i;
    return sqrt(norma);
  }

  template<class T>
  argumentType<T> somma(argumentType<T> const & X,argumentType<T> const & Y){
     argumentType<T> temp(X.size());
     for(unsigned int i=0; i<X.size();i++) temp[i]=X[i]+Y[i];
     return temp;
  }

  template<class T>
  argumentType<T> prod(argumentType<T> const & X,T const & Y){
     argumentType<T> temp(X.size());
     for(unsigned int i=0; i<X.size();i++) temp[i]=X[i];
     return temp;
  }


  template<class T=double>
  FixedPointStatus<T> fixedPoint(NonLinSys<T> const & fSys,
			  argumentType<T> & x,
			  const FixedPointOptions<T> opt={1.e-03,1000,1.0})
  {
    if(fSys.numEq()!=x.size()) return  FixedPointStatus<T>{0,0,false};
    // Check if we are already at zero
    auto res=fSys.residual(x);
    T resNorm=norm(res);
    if(resNorm<=opt.minRes) return FixedPointStatus<T>{resNorm,0,true};
    // Start Newton iterations
    unsigned int iter(0);
    do{
      // compute step=- alpha*F(x)
      // A better code would check if the Jacobian is singular
      // and eventualy regularize it.
      //
      auto step=prod(res,opt.alpha);
      x =somma(x, step);
      res = fSys.residual(x);
      resNorm = norm(res);
      ++iter;
    }  
    while(iter<opt.maxIter && resNorm>opt.minRes);
    // Check what happened
    bool converged(resNorm<=opt.minRes);
    return FixedPointStatus<T>{resNorm,iter,converged};
  }
}
#endif
