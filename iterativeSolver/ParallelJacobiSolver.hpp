#ifndef PARALLEL_JACOBI_SOLVER_HPP
#define PARALLEL_JACOBI_SOLVER_HPP

namespace Solver{


template<typename T>
class ParallelJacobiSolver
        :public IterativeSolver<T>{

   public : 
   
   ParallelJacobiSolver(CSRMatric<T> const & A, Vector<T> &x , Vector<T> const & b);

   void setNumThread(int const & nThread);
   
   void solve();
   
   
   private :
   
   static void 
   computeResidual
      ( CSRMatric<T> const & A
      , Vector<T> const    & x
      , Vector<T> const    & b
      , Vector<T> & Residual
      , int const Offset
      , int const Endset );
      
   static void
   updateInitialSolution
      ( Vector<T> & x
      , Vector<T> const & Residual
      , Vector<T> const & DiagInverse
      , int const Offset
      , int const Endset );
         
   Vector<T> DiagInverse;
   bool isProblemset;
   int _nThread;  
   T Norm;
   int iter;
 };
}//namespace Solver

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////

namespace Solver{

template<typename T>
ParallelJacobiSolver<T>::ParallelJacobiSolver
                      (CSRMatric<T> const & A, Vector<T> & x, Vector<T> const & b)
: IterativeSolver<T>(A,x,b)
, DiagInverse(A.get_nrows())
,_nThread(1)
,Norm(0)
,iter(0)
{

 

}

template<typename T>
void 
ParallelJacobiSolver<T>::setNumThread(int const & nThread){
   _nThread = nThread;
}

template<typename T>
void
ParallelJacobiSolver<T>::solve(){

omp_set_num_threads(_nThread);
#pragma omp parallel
{

  int tid = omp_get_thread_num();
  CSRMatric<T> const & A (this->_A);
  
  UniformPartition<int> P(0, A.get_nrows(), _nThread);
  int Offset = P.begin(tid+0);
  int Endset = P.begin(tid+1);
  
  for(int icsr(Offset); icsr<Endset; icsr++){
     for(size_t  jcsr=A.RowPtr[icsr+0]
         ;       jcsr<A.RowPtr[icsr+1]
         ;       jcsr++)
      {  
             if(A.ColIndices[jcsr] == icsr ){
                DiagInverse[icsr] = 1./A.Values[jcsr];
             }
      }
  }

}





 
  
 omp_set_num_threads(_nThread);
 CSRMatric<T> const & A (this->_A);
 Vector<T> const & b    (this->_b);
 Vector<T> & Residual   (this->_Residual); 
 Vector<T> & x          (this->_x);
 T Norm                 (this->Norm);
 int iter               (this->iter);
 
do{

  #pragma omp parallel
  {
  
 
 int tid = omp_get_thread_num();
  UniformPartition<int> P(0, A.get_nrows(), _nThread);
  int Offset = P.begin(tid+0);
  int Endset = P.begin(tid+1);
  
  computeResidual
     ( A
     , x
     , b
     , Residual
     , Offset
     , Endset ); 
  
  #pragma omp barrier

  UniformPartition<int> Q(0, A.get_nrows(), _nThread);
  int offset = Q.begin(tid+0);
  int endset = Q.begin(tid+1);
   
    updateInitialSolution
      ( x
      , Residual
      , DiagInverse
      , offset
      , endset
      );
   }
 
 
   Norm = this->_Residual.MaxNorm(); 
   iter++;
   } while( this->stoppingCtriterion( Norm , iter ) == false );
        
   std::cout << "Iteration  " << iter << "  Residual norm " << Norm << std::endl;   

}

template<typename T>
void 
ParallelJacobiSolver<T>
::computeResidual
( CSRMatric<T> const & A
, Vector<T> const    & x
, Vector<T> const    & b
, Vector<T> & Residual
, int const Offset
, int const Endset )
{
 
    		
  for(int icsr(Offset); icsr<Endset; icsr++){
     T sum(0);
     for(size_t  jcsr=A.RowPtr[icsr+0]
         ;       jcsr<A.RowPtr[icsr+1]
         ;       jcsr++)
      {  
         sum += A.Values[jcsr]*x[A.ColIndices[jcsr]];
      }
      Residual[icsr] = sum;
      Residual[icsr] = b[icsr] - Residual[icsr]; 
  }   		
    		
}

template<typename T>
void 
ParallelJacobiSolver<T>
::updateInitialSolution
( Vector<T> & x
, Vector<T> const & Residual
, Vector<T> const & DiagInverse
, int const Offset
, int const Endset )
{

   for(int icsr(Offset); icsr<Endset; icsr++){
      x[icsr] = x[icsr] + (0.666)*Residual[icsr]*DiagInverse[icsr];
   }

}    		

}//namespace Solver

#endif
