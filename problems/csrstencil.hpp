#ifndef CSRSTENCIL_HPP
#define CSRSTENCIL_HPP

#include "../Mat/COOBuffer.hpp"
#include "../Mat/csrmatric.hpp"
#include "problem.hpp"
 
template<typename T>
class PoissonBVP2D : public Problem<T> {
   
public :

 PoissonBVP2D(int _n1d)
   :Problem<T>(_n1d)
   ,n1d(_n1d)   
   {}  
 
private :    
 void
 setMatric( COOBuffer<T> & A  ) const override {
  
  for(int k = 0; k < n1d*n1d ; ++k){
     A.setEntry(k, k, -4.0);
     
         
    	if(k%n1d > 0){
      	A.setEntry(k, k - 1, 1.0);
    	}
    	if(k%n1d < n1d - 1) {
      	A.setEntry(k, k + 1 , 1.0);
    	}
   	    if(k%(n1d*n1d) >= n1d) {
        A.setEntry(k, k - n1d, 1.0);
        }
        if(k%(n1d*n1d) < n1d*n1d - n1d) {
        A.setEntry(k, k + n1d , 1.0);
        }
       }

       A.finalize();
     
     }


   
 void
  setRhsVector
    ( Vector<T> & _rhs
    , Vector<T> & _exactSol ) const override {
   
    double x, y; double phi;
    double sol;

    double h = 1.0/(n1d + 1);
    double h2inv = (h*h);

	for(int k=0; k <  n1d*n1d; k++){  
    	 x = (k%n1d + 1)*h;
   	 y = (k%(n1d*n1d)/n1d + 1)*h;
       //std::cout << k << " :  " << x << " " << y << std::endl;
      	sol = sin(M_PI*x)*sin(M_PI*y);
      	phi = ( -2*(M_PI)*(M_PI) )*sin(M_PI*x)*sin(M_PI*y);
      	_exactSol[k] = sol; 
      	_rhs[k] = h2inv*phi;
      
        }

     }

  private :
      int  n1d;
 
 }; 


#endif
