#ifndef GSiedel_SOLVER_HPP
#define GSiedel_SOLVER_HPP

#include "../Mat/csrmatric.hpp"
#include <iostream>
#include "IterativeSolver.hpp"


namespace Solver{

template <typename T>
class GSiedel
       : public IterativeSolver<T> {

  public :
    
  GSiedel(CSRMatric<T> const & A, Vector<T> &x , Vector<T> const & b);

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
GSiedel<T> :: GSiedel
                  (CSRMatric<T> const & A, Vector<T> &x , Vector<T> const & b)
:  IterativeSolver<T>( A, x ,  b)
, _DiagInverse(A.get_nrows())
{ 
  
   for(int icsr(0); icsr<A.get_nrows(); ++icsr){
        for(size_t  jcsr=A.RowPtr[icsr+0] 
            ;       jcsr<A.RowPtr[icsr+1]
            ;       ++jcsr)
	    	{
		    if( A.ColIndices[jcsr] == icsr ){
      			 _DiagInverse[icsr] = A.Values[jcsr];
              }
	        }
	    }

}
                  
template <typename T>
void
GSiedel<T> :: solve(){

    int iter =0; T Norm; 
    T C;
    do{
 
      for(int icsr(0); icsr<this->_A.get_nrows(); ++icsr){
          C = this->_b[icsr] ;
          for(size_t  jcsr=this->_A.RowPtr[icsr+0] 
              ;       jcsr<this->_A.RowPtr[icsr+1]
              ;       ++jcsr)
	    {
                 int j = this->_A.ColIndices[jcsr];
                 if( icsr != j ){
                     C = C - this->_A.Values[jcsr]*this->_x[j];  
                    }
	    }

         this->_x[icsr] = C/(_DiagInverse[icsr]);
        }
 

      this->_Residual = this->_b - (this->_A*this->_x);
       Norm = this->_Residual.MaxNorm();
       iter++;
    
      }while( this->stoppingCtriterion( Norm , iter ) == false );
 }

} //namespace Solver

#endif
