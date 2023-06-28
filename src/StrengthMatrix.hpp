#ifndef STRENGTHMATRIX_HPP
#define STRENGTHMATRIX_HPP

#include <limits>
#include "../Mat/csrmatric.hpp"
#include "UniformPartition.hpp"
#include <iostream>
#include <omp.h>

template <typename T>
class StrengthMatrix
      {

  public :
    
  StrengthMatrix(CSRMatric<T> const & A
                , int const & nThread );

  void apply();

  CSRMatric<T> const & Get_Strength(); 
   
  void
  strengthThreshold
      (T const & _strength_threshold);


  protected :

    static void 
    Strength_thread_func( int tid , void *arg);

    static void 
    Strength_Seq
          ( int tid
          , CSRMatric<T> const & A
          , std::vector<int> &  S_temp_j
          , std::vector<size_t> & tid_nnz
          , int Offset
          , int Endset
          , double const & theta); 

    CSRMatric<T> const & _A;
    CSRMatric<T>   _strengthMatrix;
    
    std::vector<int>  _S_temp_j;
    double strength_threshold;
    std::vector<size_t> tid_nnz;
    std::vector<size_t> work_space;
 
    int _nRows;
    int const & _nThread;
    };

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////


template <typename T>
StrengthMatrix<T> 
      :: StrengthMatrix(CSRMatric<T> const & A
                       , int const & nThread )
: _A(A)
, _nRows(A.get_nrows())
,_nThread(nThread)
{ }


template <typename T>
void
StrengthMatrix<T> 
     :: strengthThreshold( T  const & _strength_threshold ) 
  {
       strength_threshold = _strength_threshold;
  }
                  
template <typename T>
void
StrengthMatrix<T> 
     :: apply(){
  
  _S_temp_j.resize(_A.get_nnzs(),-1);
  
  tid_nnz.resize(  _nThread ,0);
  work_space.resize( _nThread + 1 ,0 ); 
  
  omp_set_num_threads(_nThread);  

  #pragma omp parallel
  {

    int tid = omp_get_thread_num();
  
    Strength_thread_func( tid , this); 
  }

}


template <typename T>
CSRMatric<T> const & 
StrengthMatrix<T> 
     :: Get_Strength(){

    return _strengthMatrix; 
}

template< typename T >
void
StrengthMatrix<T>
     :: Strength_thread_func( int tid , void *arg){

 
   StrengthMatrix< T > & strengthmatrix( *reinterpret_cast<StrengthMatrix<T>*> (arg) );
   CSRMatric<T> const & A ( strengthmatrix._A);
   int const & numThreads = strengthmatrix._nThread;   

   std::vector<int> &  S_temp_j  ( strengthmatrix._S_temp_j);     
   std::vector<size_t> & tid_nnz ( strengthmatrix.tid_nnz);
   std::vector<size_t> & work_space   ( strengthmatrix.work_space );
 
   double   const theta    =  strengthmatrix.strength_threshold;
 
   UniformPartition<int> P(0, A.get_nrows(), numThreads);
     
   int Offset = P.begin(tid+0);
   int Endset = P.begin(tid+1);


    Strength_Seq
       ( tid
       , A
       , S_temp_j
       , tid_nnz
       , Offset
       , Endset
       , theta );
  
             // Allocate the Strength Matrix
   #pragma omp barrier
    
   if(tid == 0){

      work_space[0] = 0;
      size_t s1 = 0;
      
      for(int i=1; i<numThreads+1; i++ ){
         s1 += tid_nnz[i-1];
         work_space[i]  = s1;
         }
     
      strengthmatrix._strengthMatrix.init(A.get_nrows(), work_space[numThreads] ,0);
      strengthmatrix._strengthMatrix.RowPtr[A.get_nrows()] = s1;   
    }
    
    #pragma omp barrier
   
    CSRMatric<T> & S ( strengthmatrix._strengthMatrix ); 
    

    int js=work_space[tid];

    for(int i=Offset
       ;    i<Endset
       ;    ++i)
       {
             S.RowPtr[i] = js;
             for(size_t jA=A.RowPtr[i+0]
                ;       jA<A.RowPtr[i+1]
                ;       ++jA)
                {
                  if( S_temp_j[jA] > -1 ){
                      S.ColIndices[js] = S_temp_j[jA];
                      js++; 
                     }
                }
              
        }


  }


template< typename T >
void
StrengthMatrix<T>
     :: Strength_Seq
          ( int tid
          , CSRMatric<T>  const & A
          , std::vector<int> &  S_temp_j
          , std::vector<size_t> & tid_nnz
          , int Offset
          , int Endset
          , double const & theta)
{


for(int i= Offset
    ;   i< Endset
    ;   ++i )
    {
       
    
    double row_scale = -std::numeric_limits<T>::max();
    for(size_t jA=A.RowPtr[i+0]
        ;      jA<A.RowPtr[i+1]
        ;    ++jA)
       {
        if( A.ColIndices[jA] != i ){ 
           row_scale = std::max( row_scale
                               , std::abs(A.Values[jA]) );
           }
        } 
 
    tid_nnz[tid] += A.RowPtr[i+1]
                  - A.RowPtr[i+0]
                  - 1;
  
    row_scale *= theta;

    for(size_t jA=A.RowPtr[i+0]
       ;       jA<A.RowPtr[i+1]
       ;      ++jA )
       {
        if( A.ColIndices[jA] != i ){
          if(std::abs( A.Values[jA])> row_scale ){
             S_temp_j[jA] = A.ColIndices[jA];
            }
            else{
                S_temp_j[jA] = -1;
                --tid_nnz[tid];
             }
           }
        }        

    }

}


#endif
