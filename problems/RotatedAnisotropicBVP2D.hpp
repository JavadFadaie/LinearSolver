#ifndef ROTATEDANISOPTROPIC_BVP_2D_HPP
#define ROTATEDANISOPTROPIC_BVP_2D_HPP

#include "../Mat/COOBuffer.hpp"
#include "../Mat/csrmatric.hpp"
#include "problem.hpp"

template<typename T>
class RotatedAnisotropicBVP2D : public Problem<T> {

  public :

 RotatedAnisotropicBVP2D
  ( int const _n1d
  , double _eps)
  :Problem<T>(_n1d)
  , n1d(_n1d)
  , eps(_eps)
  {
  std::cout << "this is problem " << std::endl; 
   double epsi(0.001);
   double s = sin(M_PI/4);
   double c = cos(M_PI/4);
   c1= -( c*c + epsi*s*s ) ;
   c2= 2*(1-epsi)*s*c;
   c3= -(s*s +epsi*c*c);  
   
   c1= 1 ;
   c2= 1;
   c3= 1; 
  }

 private :

 void
 setMatric( COOBuffer<T> & A  ) const override {

  for(int k = 0; k < n1d*n1d ; ++k){

       A.setEntry(k, k, -2*c1 -2*c3 );
          
        if(k%n1d > 0){   //Left stencil
        A.setEntry(k, k - 1, c1 );
        }

        if(k%n1d < n1d - 1) {//Right stencil 
        A.setEntry(k, k + 1 , c1 );
        }

        if(k%(n1d*n1d) >= n1d) {
        A.setEntry(k, k - n1d, c3);
        }
        if(k%(n1d*n1d) < n1d*n1d - n1d) {
        A.setEntry(k, k + n1d , c3);
        }

        
        if( (k%n1d > 0) && (k%(n1d*n1d) >= n1d) ){ //left down stencil 
        A.setEntry(k, k - n1d - 1 , c2/4.0  );
        }

        if( (k%n1d < n1d - 1) && (k%(n1d*n1d) >= n1d) ) {//Right down stencil
        A.setEntry(k, k - n1d + 1 , -c2/4.0  );
        }
     
        if((k%(n1d*n1d) < n1d*n1d - n1d) && (k%n1d > 0) ) {//Left up stencil  
        A.setEntry(k, k + n1d - 1 , -c2/4.0);
        }

        if((k%(n1d*n1d) < n1d*n1d - n1d) && (k%n1d < n1d - 1) ) {//Right Up stencil 
        A.setEntry(k, k + n1d + 1 , c2/4.0);
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
        phi = (-(c1)*(M_PI)*(M_PI) )*sin(M_PI*x)*sin(M_PI*y) + (-(c3)*(M_PI)*(M_PI) )*sin(M_PI*x)*sin(M_PI*y) + ( (c2)*(M_PI)*(M_PI) )*cos(M_PI*x)*cos(M_PI*y) ;
        _exactSol[k] = sol;
        _rhs[k] = h2inv*phi;

        }

     }


  private :
    int  n1d;
    double eps;
    double c1;
    double c2;
    double c3;

 };

#endif
