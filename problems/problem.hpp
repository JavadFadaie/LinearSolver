#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include "../Mat/Vector.hpp"
#include "../Mat/csrmatric.hpp"
#include "../Mat/COOBuffer.hpp"
#include <cstdlib>


template< typename T>
class Problem {

public : 

 Problem
   ( int n ){
    _n = n;
    _rhs.resize(n*n);
    _exactSol.resize(n*n);
    _initial.resize(n*n);
  }


  virtual
   ~Problem() {}

  void 
  generateDiscretization(){
     setMatric( _matricBuffer );
     setRhsVector( _rhs, _exactSol );
   }

  void
  randomApproxVectorInit() {
  
    for(int i(0);i<_initial.Size();++i) {
      _initial[i] = static_cast<double>(rand())
                      / static_cast<double>(RAND_MAX);
    }
  }

   CSRMatric<T> & getMatric() {
     return _matricBuffer.CSRFormat();
   }

   Vector<T> & getRhsVector(){
     return _rhs;
   }

   Vector<T> & getExactSol(){
     return _exactSol ;
   }

   Vector<T> & getApproxVector(){
    return _initial;
   }

private :

 virtual void 
 setMatric(COOBuffer<T> &  Matric ) const = 0 ;

 virtual void
 setRhsVector( Vector<T> &  RHS
            , Vector<T>  &  Exact_Sol ) const =0 ;


 

private :

  int _n;
  Vector<T> _rhs;
  Vector<T> _exactSol;
  Vector<T> _initial;    
  
  COOBuffer<T> _matricBuffer;

};

#endif

