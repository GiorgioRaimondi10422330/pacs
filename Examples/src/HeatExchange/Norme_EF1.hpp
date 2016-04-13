#ifndef Norm_EF1_h
#define Norm_EF1_h

template <typename T>
T norm_EF1(const T dPrec, const T dNext,const T h, const int Tipo){ 
//se Tipo è 0 sto usando la norma L2 altrimenti la norma H1 nell'intervallo [x(i),x(i+1)]
    return (dPrec*dPrec+dNext*dNext+dPrec*dNext)*h/3+Tipo*(dPrec*dPrec+dNext*dNext-2*dPrec*dNext)/h;
	  // norma L2 al quadrato                                  seminorma H1 al quadrato                                       
}





#endif

