#ifndef TRANSPOSE_HPP
#define TRANSPOSE_HPP

#include <iostream>
#include <omp.h>

template<typename T>
class Transpose{

     public:
      Transpose(CSRMatric<T> const & A
              , int  const & numCols
              , int  const & nThread );

      void apply();

      CSRMatric<T> const & Get_Transpose();

      private:

      static void 
      Transpose_Thread_Func(int tid, void *arg);

      static void
      Transpose_func
         ( CSRMatric<T> const & A
         , CSRMatric<T> & AT
         , int Offset
         , int Endset
         );

     static void
     InsetTranspose
         ( std::vector< CSRMatric<T> > & TransChunk
         , CSRMatric<T> & Trans
         , int Offset
         , int Endset );


      CSRMatric<T> const & _A;
      CSRMatric<T>  _A_T;

      std::vector< CSRMatric<T> > _TransChunk;
      int const & _numCols;
      int const & _nThread;

};

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////

template<typename T>
Transpose<T> 
     :: Transpose
           (CSRMatric<T> const & A
           ,int const &  numCols
           ,int const & nThread )
:_A(A)
,_numCols(numCols)
,_nThread(nThread)
{}

template<typename T>
void 
Transpose<T>::apply(){
    _TransChunk.resize( _nThread );  
    _A_T.init(_numCols, _A.get_nnzs(), 0);
    omp_set_num_threads(_nThread);
  
   #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      Transpose_Thread_Func(tid ,this);
    }
  
  
}          

template<typename T>
CSRMatric<T> const & 
Transpose<T>::Get_Transpose(){

    return  _A_T;
}

template<typename T> 
void 
Transpose<T> :: Transpose_Thread_Func(int tid, void *arg){

   Transpose<T> & transpose( * reinterpret_cast< Transpose<T> *> (arg) );
   int const & numThreads = transpose._nThread;
   int const & numCols = transpose._numCols;
  
   CSRMatric<T> & Trans   ( transpose._A_T);
   CSRMatric<T> const & A       ( transpose._A);
     
   std::vector< CSRMatric<T> > & TransChunk (transpose._TransChunk);

   UniformPartition<int> P(0, A.get_nrows(), numThreads);

   int Offset = P.begin(tid+0);
   int Endset = P.begin(tid+1);
    
    
   TransChunk[tid].init( numCols
                       , A.RowPtr[Endset]
                       - A.RowPtr[Offset]
                       , 0 );
             // Initialize the Transpose matrics 
     
   Transpose_func
     ( A
     , TransChunk[tid]
     , Offset
     , Endset
     );
  
     #pragma omp barrier
  

    UniformPartition<int> Q(0, Trans.get_nrows(), numThreads);   

    InsetTranspose
     ( TransChunk
     , Trans
     , Q.begin(tid)
     , Q.end(tid) );
    
}


template <typename T>
void
Transpose<T>
      ::Transpose_func
         ( CSRMatric<T> const & A
         , CSRMatric<T>  & AT
         , int Offset
         , int Endset
         )
{
      
     int nRowsAT = AT.get_nrows();

       for(int  i=0
        ;       i<= nRowsAT
        ;       i++)
          AT.RowPtr[i] = 0;


    for(int  i=A.RowPtr[Offset]
         ;   i<A.RowPtr[Endset] 
         ;   i++ )
          AT.RowPtr[A.ColIndices[i]+1]++;



    for(int  i=0
         ;   i< nRowsAT
         ;   i++)
          AT.RowPtr[i+1] += AT.RowPtr[i];


     for(int iRow=Offset
        ;    iRow<Endset
        ;    iRow++){
          for(size_t csrIndex=A.RowPtr[iRow+0]
              ;      csrIndex<A.RowPtr[iRow+1]
              ;      ++csrIndex ){
               int iCol = A.ColIndices[csrIndex];
               AT.ColIndices[ AT.RowPtr[iCol] ] = iRow;
               AT.Values[AT.RowPtr[iCol]] = A.Values[csrIndex];
               AT.RowPtr[iCol]++;
              }
       }

   
     for(int i=nRowsAT
        ;    i>0
        ;    i--)
        AT.RowPtr[i] = AT.RowPtr[i-1];

        AT.RowPtr[0]=0;

}

template <typename T>
void
Transpose<T>
  ::InsetTranspose
     ( std::vector< CSRMatric<T> > & TransChunk
     , CSRMatric<T> & Trans
     , int Offset
     , int Endset )
 {

     size_t Trans_csrIndex(0);
     for(size_t i(0)
        ;       i<TransChunk.size()
        ;       i++)
      {
        Trans_csrIndex += TransChunk[i].RowPtr[Offset];
      }

     
     for(int irow(Offset)
        ;    irow<Endset
        ;    irow++)
      { 
         Trans.RowPtr[irow] = Trans_csrIndex;
         for(size_t i(0); i<TransChunk.size(); i++)
         {
           CSRMatric<T>  & AT(TransChunk[i]);

           for(size_t csrIndex(AT.RowPtr[irow+0]) 
              ;       csrIndex<AT.RowPtr[irow+1]
              ;       csrIndex++ )
            {     
              Trans.ColIndices[Trans_csrIndex] = AT.ColIndices[csrIndex];
              Trans.Values[Trans_csrIndex]     = AT.Values[csrIndex];       
              Trans_csrIndex++;
            }
         }
      } 

      Trans.RowPtr[Endset] = Trans_csrIndex; 

 }

#endif
