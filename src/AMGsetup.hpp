#ifndef AMGSETUP_HPP
#define AMGSETUP_HPP

#include <iostream>
#include "StrengthMatrix.hpp"
#include "RugeCoarsen.hpp"
#include "Interpolation.hpp"
#include "Transpose.hpp"
#include "Multiplication.hpp"
#include <memory>

template<typename T>
class AMGSetUp{

     public:
     AMGSetUp(CSRMatric<T> const & A);

     void setNumThread( int const & nThread);

     void apply();

     CSRMatric<T> const & Prolongation();

     CSRMatric<T> const & Restrict();

     CSRMatric<T> const & CoarseOperator();
                 
     int  & CoarseSize();
  
     int   FineSize();

     ~AMGSetUp(){}

     
    private:

     CSRMatric<T> const & _A;
     int  _nThread;
     std::unique_ptr< Interpolation<T> >  interpolation;
     std::unique_ptr< Transpose<T> >      restriction;
     std::unique_ptr< Multiplication<T> > CoarseOp;

};

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////

template<typename T>
AMGSetUp<T> :: AMGSetUp(CSRMatric<T> const & A):
_A(A){}

template<typename T>
void
AMGSetUp<T> :: setNumThread(int const & nThread){
    _nThread = nThread;
}

template<typename T>
void 
AMGSetUp<T>:: apply(){
    
   
    StrengthMatrix<T> strengthMatrix(_A, _nThread);
    strengthMatrix.strengthThreshold(0.25);
    strengthMatrix.apply();
    
      
    RugeCoarsen<T> rugeCoarsen(strengthMatrix.Get_Strength(), _nThread);

    rugeCoarsen.apply(); 
    rugeCoarsen.Second_Phase_Ruge(); 
   
    interpolation = std::unique_ptr< Interpolation<T> >
                         ( new Interpolation<T> 
                                 (  _A
                                 , strengthMatrix.Get_Strength()
                                 , rugeCoarsen.CF_Marker()
                                 , _nThread  )
                          );

    
    interpolation->apply();
    
  

    restriction  = std::unique_ptr< Transpose<T> > 
                        ( new Transpose<T> 
                                  ( interpolation->Get_Interpolation()
                                  , interpolation->NumCoarse()
                                  , _nThread )
                        );
    restriction->apply();

     Multiplication<T> M1( _A
                         , interpolation->Get_Interpolation()
                         , interpolation->NumCoarse() );
     M1.setNumThreads(_nThread );
     M1.apply();
 

     CoarseOp = std::unique_ptr< Multiplication<T> > 
                    ( new Multiplication<T> 
                         ( restriction->Get_Transpose()
                         , M1.Get_Multiplication()
                         , interpolation->NumCoarse() ) 
                    );
    
      CoarseOp->setNumThreads(_nThread );
      CoarseOp->apply();
   
  
}

template<typename T>
CSRMatric<T> const &
AMGSetUp<T> :: Prolongation(){

   return interpolation->Get_Interpolation();
}

template<typename T>
CSRMatric<T> const & 
AMGSetUp<T> :: Restrict(){

   return restriction->Get_Transpose();
}

template<typename T>
CSRMatric<T> const &
AMGSetUp<T> :: CoarseOperator(){

   return CoarseOp->Get_Multiplication();
}

template<typename T>
int &
AMGSetUp<T> :: CoarseSize(){

   return interpolation->NumCoarse();
}

template<typename T>
int 
AMGSetUp<T> :: FineSize(){

   return _A.get_nrows();
}



#endif
