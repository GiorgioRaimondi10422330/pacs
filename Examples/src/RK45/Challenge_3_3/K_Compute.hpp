#ifndef K_COMPUTE_HPP
#define K_COMPUTE_HPP

#include "fixedPoint.hpp"
#include "nonLinSys2.hpp"
#include <functional>

using namespace std;
using namespace NonLinearSystems;

template<class T>
vector<argumentType<T>> conv(vector<T> const & K, size_t N, size_t M){
	vector<argumentType<T>> temp (N,argumentType<T> (M));
	for(size_t i=0; i<N; i++){
		for(size_t j=0;j<M; j++){
			temp[i][j]=K[i*M+j];
		}
	}
	return temp;
}


template<class T=double>
void K_compute(
  vector<argumentType<T>> & K,
  NonLinSys<T> const & fun,
  size_t const dimSys,
  size_t const dimK,
  argumentType<T> const & y0,
  T const t0,
  vector<std::vector<T>> const & A,
  std::vector<T> const & c,
  T h
)
{
    argumentType<T> V(dimSys*dimK,0.);
    NonLinSys<T> F;
    for(size_t i=0; i<dimK; i++){
    	for(size_t j=0;j<dimSys; j++){
    		auto fK=[&A, &c , &y0, t0, dimSys, dimK, h, i,j,&fun](argumentType<T> const & INGR){
    			argumentType<T> Y(dimSys+1);
    			Y[dimSys]=t0*c[i];
    			for(size_t q=0;q<dimSys;q++){
    				Y[q]=y0[q];
    			}
    			for(size_t r=0; r<i+1;r++){
    				for(size_t q=0;q<dimSys;q++){
    					Y[q]+=h*A[i][r]*INGR[r*dimSys+q];
    				}
    			}
    			auto ff=fun[j];
         		return ff(Y);
    		};
    		F.addToSystem(fK);
    	}
    }
    auto result=fixedPoint(F,V);
    cout<<"Residuo="<<result.residual<<"  Iterazioni="<<result.iterations;
    if(result.converged) cout<<"    Il metodo CONVERGE\n";
    else cout<<"    Il metodo NON CONVERGE\n";
    K=conv(V,dimK,dimSys);

}









#endif
