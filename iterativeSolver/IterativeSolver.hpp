#ifndef ITERATIVESOLVER_HPP
#define ITERATIVESOLVER_HPP

#include "../Mat/csrmatric.hpp"
#include <iostream>

namespace Solver{

template <typename T>
class IterativeSolver {

  public :
    IterativeSolver(CSRMatric<T> const & A, Vector<T> &x , Vector<T> const & b)
    : _A(A)
    , _x(x)
    , _b(b)
    , _Residual(_A.get_nrows())
    {}

    ~IterativeSolver()
    {}

   virtual void solve() = 0; 

    void
    setMaxIteration (int  _maxIteration){
      maxIteration = _maxIteration;
    } 

    void
    setMinIteration (int  _minIteration){
      minIteration = _minIteration; 
    } 
   
    void
    setRelativeTolorance ( T  _relativeTolorance ){
      relativeTolorance = _relativeTolorance;
    } 
   

    bool
    stoppingCtriterion ( T & _residualNorm, int iter ){
         
      iteration = iter;

       if(iter == 0){
          norm0 = _residualNorm;
          }
  
     if( (_residualNorm/norm0 < relativeTolorance ) 
          || ( iter > maxIteration ) ){
             isStoppingCondionStatified = true;  
        }

      return isStoppingCondionStatified;

     } 
   

    bool
    stoppingCtriterionError ( T & _ErrorNorm, int iter ){

      iteration = iter;

     if( (_ErrorNorm < relativeTolorance )
          || ( iter > maxIteration ) ){
             isStoppingCondionStatified = true;
        }

      return isStoppingCondionStatified;

     }

   
   T 
   ErrorInfNorm ( Vector<T> & ExactSol, Vector<T> & Approx ){
      
        Vector<T> Error( ExactSol.Size() );
      
        Error = ExactSol - Approx ;

       return Error.MaxNorm();
    }


   int const  & numInteration(){

     if( isStoppingCondionStatified ){
        return iteration;
        }
    }

  protected :
 
    CSRMatric<T> const & _A;
	Vector<T> & _x;
	Vector<T> const & _b;

	Vector<T> _Residual;
	bool isProblemset; 
   
    int maxIteration;
    int minIteration;
    T relativeTolorance;
    T norm0 ;
    int  iteration;
    bool isStoppingCondionStatified = false; 
   };

} //namespace Solver

#endif
