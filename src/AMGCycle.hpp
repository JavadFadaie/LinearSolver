#ifndef AMG_CYCLE_HPP
#define AMG_CYCLE_HPP
 
#include "AMGsetup.hpp"
#include "../iterativeSolver/IterativeSolver.hpp"
#include "../Mat/Vector.hpp"
#include <cmath>
#include <limits>

template <typename T>
class AMGCycle 
         {

      using IterativeSolver = Solver::IterativeSolver<T>;
      using GSiedel = Solver::GSiedel<T>;

       public :
        AMGCycle
           ( CSRMatric<T> const & A
           , T         &  strength_threshold
           , int const & num_Pre_Smoothing
           , int const & num_Post_Smoothing
           , int const & gridSizeThreshold
           , int iLevel
           , int nLevelMax
           , int nThread
           );

        ~AMGCycle(){
            delete _postSmoother;
            delete _preSmoother;
         }
        void V_Cycle
              ( Vector<T> & x
              , Vector<T> const & b );

       private :

       CSRMatric<T> const & _operator;
       AMGCycle<T>  * _nextLevelAMG;
       std::unique_ptr < AMGSetUp<T> > _amgSetUp;

       IterativeSolver *  _finalSolver;
       IterativeSolver *  _preSmoother;
       IterativeSolver *  _postSmoother;
 
       Vector<T> _fineResidual;
       Vector<T> _coarseResidual;
       Vector<T> _coarseError;
       Vector<T> _fineError;
       Vector<T> _Residual; 

       T   _strength_threshold;
       int _gridSizeThreshold;
       int _num_Pre_Smoothing; 
       int _num_Post_Smoothing;
       int _iLevel;    
       int _iLevelMax;
       bool _isLastLevel; 
       };


////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////

template<typename T>
AMGCycle<T> :: AMGCycle
           (CSRMatric<T> const & A
           , T         &  strength_threshold
           , int const & num_Pre_Smoothing 
           , int const & num_Post_Smoothing
           , int const & gridSizeThreshold 
           , int iLevel  
           , int nLevelMax
           ,  int nThread 
           )
: _operator(A)     
, _strength_threshold(strength_threshold)      
, _num_Pre_Smoothing(num_Pre_Smoothing)
, _num_Post_Smoothing(num_Post_Smoothing)
, _gridSizeThreshold(gridSizeThreshold )
, _iLevel(iLevel)
, _iLevelMax(std::max(0,nLevelMax-1))
, _isLastLevel(false) 
{


   std::cout <<"AMG Level:   " << _iLevel << std::endl;
   std::cout <<"number of Var:  "<< A.get_nrows() << std::endl;

  if(A.get_nrows()<_gridSizeThreshold 
         || _iLevel == _iLevelMax ) {
    _isLastLevel =true;
   }
   else{

    _amgSetUp = std::unique_ptr < AMGSetUp<T> > 
                          ( new AMGSetUp<T>(A) );

    _amgSetUp-> setNumThread(  nThread);
    _amgSetUp->apply();

    _nextLevelAMG = new AMGCycle<T>
          ( _amgSetUp->CoarseOperator()
          , _strength_threshold
          , _num_Pre_Smoothing
          , _num_Post_Smoothing
          , _gridSizeThreshold
          , _iLevel+1
          , nLevelMax
          ,  nThread
          );

        _fineResidual.resize(_amgSetUp->FineSize());
        _coarseResidual.resize(_amgSetUp->CoarseSize());
        _coarseError.resize(_amgSetUp->CoarseSize());
        _fineError.resize(_amgSetUp->FineSize());
    }
   
}

template<typename T>
void
AMGCycle<T> 
   ::V_Cycle
   ( Vector<T> & x
   , Vector<T> const & b) {


    if(_isLastLevel){
        _finalSolver = new  GSiedel (_operator,x,b);
          _finalSolver->setRelativeTolorance(0.0001);
          _finalSolver->setMaxIteration(50);
         _finalSolver->solve();
     }
     else
     {
       try
      {
        {
          _preSmoother = new GSiedel  (_operator,x,b);
          _preSmoother->setMaxIteration(3);
          _preSmoother->solve();
         }

        {
          _fineResidual = b - (_operator*x);  
        }
       

        _coarseResidual =(_amgSetUp->Restrict())*_fineResidual;
      
        _coarseError.set(0);
       }

        catch(std::exception & e){
         std::cerr << "Level " 
                   << _iLevel << ": exception: " 
                   << e.what() << std::endl;
          throw;
          }

       _nextLevelAMG->V_Cycle(_coarseError, _coarseResidual);

      try 
      { 
         
       _fineError =(_amgSetUp->Prolongation())*_coarseError;

        x = x + _fineError;
        
        {
          _postSmoother = new GSiedel (_operator,x,b);
          _postSmoother->setMaxIteration(3);
          _postSmoother->solve();
        } 
        
      }

         catch(std::exception & e) {
          std::cerr << "Level " << _iLevel << ": exception: " << e.what() << std::endl;
          throw;
           }
   
      }


   }
          
 


#endif
