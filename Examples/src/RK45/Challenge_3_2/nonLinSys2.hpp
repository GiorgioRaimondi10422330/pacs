#ifndef __NONLINSYS_HPP__
#define __NONLINSYS_HPP__
#include <functional>
#include <vector>
namespace NonLinearSystems{

  template<class T=double>
  using argumentType=std::vector<T>;

  template<class T=double>
  using returnType=std::vector<T> ;



  //! type of the contained functions (C++11 extension)
  template<class T=double>
  using  funType=std::function<T (std::vector<T> const &)> ;
  
  template<class T=double>
  class NonLinSys{
   public:

    NonLinSys()=default;
    //! I use the synthetic copy constructor (C++11).
    NonLinSys(NonLinSys const &)=default;
    //! I use the synthetic copy-assignment.
    NonLinSys & operator =(NonLinSys const &)=default;
    
    //! Compute \f$ -F(x) \f$
    returnType<T> residual(argumentType<T> const & x) const;
    //! Compute \f$ F(x) \f$.
    returnType<T> operator()(argumentType<T> const & x) const;
    //! Add a function to the system.
    void addToSystem(funType<T> const & f);
    //! Number of equations.
    unsigned int numEq()const;
    //! Returns the i-th function.
    funType<T>  & operator [](unsigned int i){return M_funs[i];}
    //! Returns the i-th function.
    funType<T>  operator [](unsigned int i) const {return M_funs[i];}
    //! Returns \f$ ||f(x)|| \f$
    T norm(argumentType<T> const & x) const;
  private:
    std::vector<funType<T>> M_funs; 
  };
    
template<class T>
returnType<T>
NonLinSys<T>::operator()(
  argumentType<T> const & x
  ) const
   {
    returnType<T> tmp(M_funs.size());
    unsigned int j=0;
    for (auto i: M_funs)tmp[j++]=i(x);
    return tmp;
   }


template<class T>
returnType<T>
NonLinSys<T>::residual(argumentType<T> const & x) const{
  returnType<T> tmp(M_funs.size());
  unsigned int j=0;
  for (auto i: M_funs)tmp[j++]=-i(x);
  return tmp;
}

template<class T>
void
NonLinSys<T>::addToSystem(funType<T> const & f){
  M_funs.push_back(f);
}

template<class T>
unsigned int
NonLinSys<T>::numEq()const{
  return M_funs.size();
}

template<class T>
T NonLinSys<T>::norm(argumentType<T> const & x) const{
  auto tmp=this->operator()(x);
  T norm(0);
  for(auto i: tmp) norm+=i*i;
  return sqrt(norm);
}







} // End of namespace
#endif
