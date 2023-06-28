#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include "../Mat/csrmatric.hpp"
#include <iostream>
#include "IterativeSolver.hpp"

namespace Solver {

template <typename T>
class JacobiSolver
       : public IterativeSolver<T> {

  public :
    
  JacobiSolver(CSRMatric<T>  & A, Vector<T> &x , Vector<T> const & b);


  void solve();

  

   private :

	Vector<T> _DiagInverse;
	bool isProblemset; 

    };

} //namespace Solver

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////

namespace Solver {

template <typename T>
JacobiSolver<T> :: JacobiSolver
                  (CSRMatric<T>  & A, Vector<T> &x , Vector<T> const & b)
:  IterativeSolver<T>(A, x ,b)
, _DiagInverse(A.get_nrows())
{ 
  
   for(int icsr(0); icsr<A.get_nrows(); ++icsr){
        for(size_t  jcsr=A.RowPtr[icsr+0] 
            ;       jcsr<A.RowPtr[icsr+1]
            ;       ++jcsr)
	  {
		    if( A.ColIndices[jcsr] == icsr ){
      			    _DiagInverse[icsr] = 1./A.Values[jcsr];
            	      }
	   }
     }

}
                  
template <typename T>
void
JacobiSolver<T> :: solve(){

    int iter =0; T Norm; 
    T W = 0.666; 

    do{
 
    this->_Residual = this->_b - (this->_A*this->_x);
    for(int icsr(0); icsr<this->_A.get_nrows(); ++icsr){
		   this->_x[icsr] = this->_x[icsr] + W*this->_Residual[icsr]*_DiagInverse[icsr];
		 }

	  Norm = this->_Residual.MaxNorm();
	  iter++;

      }while( this->stoppingCtriterion( Norm , iter ) == false );
        std::cout << "Iteration  " << iter << "  Residual norm " << Norm << std::endl;

}

} //namespace Solver

#endif
