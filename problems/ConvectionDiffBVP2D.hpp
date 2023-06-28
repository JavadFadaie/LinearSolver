#ifndef CONVECTIONDIFFUSION_BVP_2D_HPP
#define CONVECTIONDIFFUSION_BVP_2D_HPP
 
#include "../Mat/COOBuffer.hpp"
#include "../Mat/csrmatric.hpp"
#include "problem.hpp"


template<typename T>
class ConvectionDiffBVP2D : public Problem<T> {

  public :

 ConvectionDiffBVP2D 
  ( int const _n1d
  , double _eps)
  :Problem<T>(_n1d)
  , n1d(_n1d)
  , eps(_eps)
  {}

  private :

  void
  setMatric( COOBuffer<T> & A  ) const override {
   
    double h = 1.0/( n1d + 1);
    double mu_x=0;
    double mu_y=0;
    double Sigma=0;

   if( a*h > eps ){ mu_x = eps/2*a*h; }
   if( b*h > eps ){ mu_y = eps/2*b*h; }

   if( a*h < - eps ){ mu_x = 1 +  eps/2*a*h; }
   if( b*h < - eps ){ mu_y = 1 +  eps/2*b*h; }

   if(  a*h < eps ){ mu_x = 1/2; }
   if(  b*h < eps ){ mu_y = 1/2; }

   double rightNeighbour = - eps + a*h*(mu_x - 0);
   double leftNeighbour  = - eps + a*h*(mu_x - 1);
   double downNeighbour  = - eps + b*h*(mu_y - 1);
   double upNeighbour    = - eps + b*h*(mu_y - 0);

   Sigma = rightNeighbour +leftNeighbour + downNeighbour + upNeighbour;

  for(int k = 0; k < n1d*n1d ; ++k){

       A.setEntry(k, k, -Sigma );

        if(k%n1d > 0){
        A.setEntry(k, k - 1, leftNeighbour );
        }
        if(k%n1d < n1d - 1) {
        A.setEntry(k, k + 1 , rightNeighbour );
        }
        if(k%(n1d*n1d) >= n1d) {
        A.setEntry(k, k - n1d, downNeighbour  );
        }
        if(k%(n1d*n1d) < n1d*n1d - n1d) {
        A.setEntry(k, k + n1d , upNeighbour );
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
        phi = ( 2*( eps)*(M_PI)*(M_PI)*sin(M_PI*x)*sin(M_PI*y)
                            + a*(M_PI)*cos(M_PI*x)*sin(M_PI*y)
                            + b*(M_PI)*sin(M_PI*x)*cos(M_PI*y) );
        _exactSol[k] = sol;
        _rhs[k] = h2inv*phi;

        }

     }

  private:
  double a = sin(M_PI/8);
  double b = cos(M_PI/8);


  int  n1d;
  double eps;
};
#endif
