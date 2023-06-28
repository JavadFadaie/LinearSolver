#ifndef ANISOTROPICPOISSON_BVP_2D_HPP
#define ANISOTROPICPOISSON_BVP_2D_HPP

#include "../Mat/COOBuffer.hpp"
#include "../Mat/csrmatric.hpp"
#include "problem.hpp"

template<typename T>
class AnisotropicPoissonBVP2D : public Problem<T> {

  public :

 AnisotropicPoissonBVP2D
  ( int const _n1d
  , double _eps)
  :Problem<T>(_n1d)
  , n1d(_n1d)
  , eps(_eps)
  {}

 private :

 void
 setMatric( COOBuffer<T> & A  ) const override {

  for(int k = 0; k < n1d*n1d ; ++k){
    
       A.setEntry(k, k, -2.*(1. + eps) );

        if(k%n1d > 0){
        A.setEntry(k, k - 1, 1.0 * eps );
        }
        if(k%n1d < n1d - 1) {
        A.setEntry(k, k + 1 , 1.0 * eps );
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

        sol = sin(M_PI*x)*sin(M_PI*y);
        phi = (-(1. + eps)*(M_PI)*(M_PI) )*sin(M_PI*x)*sin(M_PI*y) ;
        _exactSol[k] = sol;
        _rhs[k] = h2inv*phi;

        }

     }
  
 
  private :
    int  n1d;
    double eps;


 };
#endif
