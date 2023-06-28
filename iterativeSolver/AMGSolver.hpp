#ifndef AMGSOLVER_HPP
#define AMGSOLVER_HPP

 #include "../src/AMGCycle.hpp"
 #include <memory>
 
namespace Solver{

template<typename T>
class AMGSolver
          : public IterativeSolver<T> {


    public :
      AMGSolver(CSRMatric<T> const & A , Vector<T> &x , Vector<T> const & b);


      ~AMGSolver(){}
      
      void strengthThreshlod( T  const  strenght);
   
      void setNumThreads( int nThread);

      void numPreSmoothing( int const  num);
 
      void numPostSmoothing(int const  num);

      void gridThresholdSize( int const threshold);

      void maxLevelAMG(int const  maxlevel);

      void  solve( );

    private :


     T strength_threshhld;
     int num_pre_smoothing;
     int num_post_smoothing;
     int gridSizeThreshold;
     int nLevelMax;
     int nThread;
     std::unique_ptr < AMGCycle<T> > amg;
 
};

} //namespace Solver

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////

namespace Solver {

template <typename T>
AMGSolver<T> :: AMGSolver
                 (CSRMatric<T> const & A , Vector<T> &x , Vector<T> const & b)
: IterativeSolver<T>( A, x, b)
{ 


  }

template< typename T>
void
AMGSolver<T> :: strengthThreshlod(T  const  strength){
   strength_threshhld = strength;
}

template<typename T> 
void 
AMGSolver<T> :: numPreSmoothing( int const  num){
    num_pre_smoothing = num;
}

template<typename T>
void 
AMGSolver<T> :: numPostSmoothing(int const  num){
    num_post_smoothing = num;
}

template<typename T>
void 
AMGSolver<T> :: gridThresholdSize( int const  threshold){
    gridSizeThreshold = threshold;
}

template<typename T>
void
AMGSolver<T> :: maxLevelAMG(int const  maxlevel){
    nLevelMax = maxlevel;
}

template<typename T>
void
AMGSolver<T> :: setNumThreads( int _nThread ){
     nThread = _nThread;
}


template<typename T>
void 
AMGSolver<T> :: solve(  ){

   int iLevel(0);
 

  amg = std::unique_ptr< AMGCycle<T> >
          ( new AMGCycle<T>(  this->_A
                        , strength_threshhld
                        , num_pre_smoothing
                        , num_post_smoothing
                        , gridSizeThreshold
                        , iLevel
                        , nLevelMax
                        , nThread
                        )
           );


    int iter(0); double norm(1);
    do{
     amg->V_Cycle( this->_x , this->_b );
     this->_Residual = this->_b - (this->_A*this->_x);
     norm =  this->_Residual.MaxNorm();
     std::cout << iter << " : "<< norm  << std::endl;
     iter++;

    }while(this->stoppingCtriterionError( norm , iter ) == false );


}

} //namespace Solver








#endif
