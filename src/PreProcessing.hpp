#ifndef _PREPROCESS_HPP
#define _PREPROCESS_HPP

#include<iostream>
#include <vector>



template< typename T >
class PreProcess {

public :

  PreProcess();
        
  PreProcess
    ( CSRMatric<T> & B 
    , int const  _num_B_partition );

  void
  init
    ( CSRMatric<T> & B 
    , int const  _num_B_partition );

  void
  apply();

  size_t const * const
  row_begin_ptr
    ( const int iPartition) const;

  size_t const * const
  row_end_ptr
    ( const int iPartition) const;
  
protected:

  static void
  Row_Pointer_Thread
    ( int tid
    , void *arg);



  CSRMatric<T> & B; 
 
  std::vector<size_t> value;
  int nrows;
  int num_B_partition;
};


////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////


template< typename T >
PreProcess<T>
  ::PreProcess()
{
  value.resize(0);
  nrows = -1;
  num_B_partition = -1;
}

template< typename T >
PreProcess<T>
  ::PreProcess
      ( CSRMatric<T> & _B
      , int const  _num_B_partition )
: B(_B)
, nrows( B.get_nrows() )
, num_B_partition( _num_B_partition)
, value( (nrows)*(num_B_partition+1), 0)
{
 //std::cout << "this" <<  std::endl;
 // B.Disply();
}


template< typename T >
void
PreProcess<T>
  ::apply()
{

    value.resize( (nrows)*(num_B_partition+1), 0);
  
    omp_set_num_threads(num_B_partition);
   
   #pragma omp parallel
   {
     int tid = omp_get_thread_num();  
    
     Row_Pointer_Thread(tid , this );
   }

}

template< typename T >
size_t const * const
PreProcess<T>
  ::row_begin_ptr( const int iPartition) const
{
  return & value[iPartition*(nrows)];
}

template< typename T >
size_t const * const
PreProcess<T>
  ::row_end_ptr( const int iPartition) const
{
   return & value[ (iPartition + 1)*(nrows) ];
}

template< typename T >
void
PreProcess<T>
  ::Row_Pointer_Thread(int tid, void *arg)
{
   
   PreProcess< T >  & _preProcess( * reinterpret_cast < PreProcess< T > *> (arg) );
   
   int  num_B_partition     =  _preProcess.num_B_partition;

   CSRMatric<T> & B ( _preProcess.B ) ;
   size_t nrows = B.get_nrows();

   UniformPartition<int> P(0, nrows, num_B_partition );
   int row_offset = P.begin(tid+0);
   int row_endset = P.begin(tid+1);

   UniformPartition<int> Q(0, nrows, num_B_partition );

   for(int i(row_offset); i<row_endset; ++i){
      _preProcess.value[i] = B.RowPtr[i];
      for(int id(0); id < num_B_partition; ++id ){

         size_t nnz(0);
         for(size_t iA = B.RowPtr[i+0]; iA < B.RowPtr[i+1]; iA++ ){
              int jj =  B.ColIndices[iA];
              if( (Q.begin(id+0) <= jj) && (jj<Q.begin(id+1)) ){
                      nnz++;
               }
         }
       _preProcess.value[(id+1)*nrows + i] = _preProcess.value[id*nrows + i] + nnz;
      }
   }

}


#endif
