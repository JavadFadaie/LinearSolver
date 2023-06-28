#ifndef MULTIPLICATION_HPP
#define MULTIPLICATION_HPP

#include "PreProcessing.hpp"
#include <memory>

template<typename T>
class Multiplication{

     public : 
 
      Multiplication( CSRMatric<T> const & A
                    , CSRMatric<T> const & B
                    , int const & _C_ncols );

      void setNumThreads( int const & nThread );
                
      void apply();

      CSRMatric<T> const & Get_Multiplication();

     private: 

     static
     void Multiplication_Thread_Func(int tid, void *arg);

     static size_t
     SpGEMM_estimate_nnz
           ( CSRMatric<T> const & A
           , CSRMatric<T> const & B
           , int row_offset
           , int row_endset
           );


     static void
     SpGEMM_Sequential
          ( CSRMatric<T> const & A
     	  , CSRMatric<T> const & B
    	  , CSRMatric<T> & C
    	  , int row_offset
     	  , int row_endset
     	  , int const & C_ncols );

     static int 
     block_row_idx
            ( int const iBlock
            , int const nBlockCols);

     static  int
     block_col_idx
            ( int const iBlock
            , int const nBlockCols );

     static  int
     block_idx
            ( int const block_row_idx
            , int const block_col_idx
            , int const nBlockCols );



     std::vector< CSRMatric<T> > C_P;
     CSRMatric<T> const & _A;
     CSRMatric<T> const & _B;
     CSRMatric<T>   _C;

     std::vector<int> tid_nnz;
     std::vector<int> work_space;

     int    _nThread;
     int    num_A_partition;
     int    num_B_partition;
     size_t atomicCounter;
     int const &  C_ncols;

    std::unique_ptr < PreProcess<T> > preProcess;
    std::vector<size_t> IBlock;
};

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////


template<typename T>
Multiplication<T> :: Multiplication( CSRMatric<T> const & A
                                   , CSRMatric<T> const & B
                                   , int const & _C_ncols )
: _A(A)
, _B(B)
, C_ncols(_C_ncols )
{}

template<typename T>
void
Multiplication<T> :: setNumThreads(int const & nThread ){
   
     _nThread = nThread;
    
     num_A_partition = _nThread;
     num_B_partition = _nThread ;
}

template<typename T>
void 
Multiplication<T> :: apply(){

  

   int NumBlock = num_B_partition;

   C_P.resize(_nThread);
   tid_nnz.resize(num_A_partition);
   work_space.resize(num_A_partition+1);
   IBlock.resize(NumBlock);

   omp_set_num_threads(_nThread);
   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      Multiplication_Thread_Func( tid ,this);
   }


}


template<typename T>
CSRMatric<T> const & 
Multiplication<T> :: Get_Multiplication(){
   return _C;
}


template<typename T>
void 
Multiplication<T>
:: Multiplication_Thread_Func(int tid, void *arg){

   Multiplication<T> & multiplication( * reinterpret_cast< Multiplication<T> *> (arg) ); 

   CSRMatric<T> const & A ( multiplication._A);
   CSRMatric<T> const & B ( multiplication._B);

   std::vector< CSRMatric<T> > & C_P1 (multiplication.C_P );

   int const & numThreads        (multiplication._nThread );
   std::vector<int> & tid_nnz    (multiplication.tid_nnz  );
   std::vector<int> & work_space (multiplication.work_space);
   int nrows = A.get_nrows();   
    
   int const nBlockRows ( multiplication.num_A_partition );
   int const nBlockCols ( multiplication.num_B_partition );
   int const nBlocks    ( nBlockRows*nBlockCols );
   
   int const &  C_ncols   (multiplication.C_ncols);
   UniformPartition<int> P(0, nrows, nBlockRows);

   size_t atomicCounter  (multiplication.atomicCounter );    

   int  row_offset = P.begin(tid+0);
   int  row_endset = P.begin(tid+1);
  
    size_t  nnz_estimate
       ( SpGEMM_estimate_nnz
          ( A
          , B
          , row_offset
          , row_endset
          )
       );

   multiplication.C_P[tid].init(row_endset - row_offset, nnz_estimate, 0);

    SpGEMM_Sequential
      ( A
      , B
      , multiplication.C_P[tid]
      , row_offset
      , row_endset
      , C_ncols
      );

     
    {
     tid_nnz[tid] += multiplication.C_P[tid].RowPtr[ multiplication.C_P[tid].get_nrows() ];
    }

    #pragma omp barrier

   if(tid == 0){

      work_space[0] = 0;
      size_t s1 = 0;

      for(int i=1; i<numThreads+1; i++ ){
         s1 += tid_nnz[i-1];
         work_space[i]  = s1;
         }
   
   
      multiplication._C.init(A.get_nrows(), work_space[numThreads] ,0);
      multiplication._C.RowPtr[A.get_nrows()] = s1;
    }
   
     #pragma omp barrier
 
     CSRMatric<T> & C   (multiplication._C);

     UniformPartition<int> Q(0, nrows, nBlockRows);
     int  Offset = Q.begin(tid+0);
     int  Endset = Q.begin(tid+1);

 
     size_t js=work_space[tid];
  
    for(int i=0
       ;    i<Endset - Offset
       ;    ++i)
       {
             multiplication._C.RowPtr[i+Offset] = js;
             for(size_t jA= multiplication.C_P[tid].RowPtr[i+0]
                ;       jA<  multiplication.C_P[tid].RowPtr[i+1]
                ;       ++jA)
                {
                  
                      multiplication._C.ColIndices[js] =  multiplication.C_P[tid].ColIndices[jA];
                      multiplication._C.Values[js] =  multiplication.C_P[tid].Values[jA];
                      js++;                  
                }

        }
   
  }




template<typename T>
void
Multiplication<T>
::SpGEMM_Sequential
      ( CSRMatric<T> const & A
      , CSRMatric<T> const & B
      , CSRMatric<T> & C
      , int row_offset
      , int row_endset 
      , int const & C_ncols )
 {

  int nrows = A.get_nrows();


  Vector<T> X(C_ncols+1);
  int irow(0); int nnz(0); 

  C.RowPtr[irow] = 0;
  T A_ik; T B_kj;
 
  
  for(int i=row_offset
     ;    i<row_endset
     ;    i++)
   {
     
     for(size_t iA = A.RowPtr[i+0]
        ;       iA < A.RowPtr[i+1]
        ;       iA++)
      {
           int jj = A.ColIndices[iA]; A_ik = A.Values[iA];

           for(size_t jA=B.RowPtr[jj+0]
              ;       jA<B.RowPtr[jj+1]
              ;       jA++)
            {          
                B_kj = B.Values[jA];
                X[ B.ColIndices[jA] ] += A_ik*B_kj ;
            }  
       }
       
      for(int it=0
         ;    it<X.Size()
         ;    it++)
       {
            if( X[it] != 0){
                C.ColIndices[nnz]=it;
                C.Values[nnz] = X[it];
                ++nnz; 
              }
       }  
       
       X.clear();            
       ++irow;
       C.RowPtr[irow] = nnz;
     }

}

template<typename T>
size_t
Multiplication<T> 
:: SpGEMM_estimate_nnz
     ( CSRMatric<T> const & A
     , CSRMatric<T> const & B
     , int row_offset
     , int row_endset
     )
{

    size_t nnz_estimate(0);
   for(int i(row_offset)
      ;    i<row_endset
      ;    i++)
    {
       size_t row_nnz_estimate(0);
       for(int icsr=A.RowPtr[i+0]
          ;    icsr<A.RowPtr[i+1]
          ;    icsr++){
       
           int jcols = A.ColIndices[icsr];
           row_nnz_estimate += B.RowPtr[jcols+1]
                             - B.RowPtr[jcols+0];
          }
         
         nnz_estimate += row_nnz_estimate;   
     }

 
   return nnz_estimate;
}

template<typename T>
int
Multiplication<T>
:: block_row_idx
      ( int const iBlock
      , int const nBlockCols)
{

    return iBlock/nBlockCols;
}

template<typename T>
int
Multiplication<T>
:: block_col_idx
     ( int const iBlock
     , int const nBlockCols )
{
     return iBlock % nBlockCols;
}


template<typename T>
int
Multiplication<T>
:: block_idx
     ( int const block_row_idx
     , int const block_col_idx
     , int const nBlockCols )
{
    return block_row_idx* nBlockCols + block_col_idx;
}

/* 
template<typename T>
void 
Multiplication<T>
:: Multiplication_Thread_Func(int tid, void *arg){

   Multiplication<T> & multiplication( * reinterpret_cast< Multiplication<T> *> (arg) ); 

   CSRMatric<T> & A ( multiplication._A);
   CSRMatric<T> & B ( multiplication._B);

   std::vector< CSRMatric<T> > & C_P1 (multiplication.C_P );

   int const nThread             (multiplication._nThread );
   std::vector<int> & tid_nnz    (multiplication.tid_nnz  );
   std::vector<int> & work_space (multiplication.work_space);
   int nrows = A.get_nrows();   
    
   int const nBlockRows ( multiplication.num_A_partition );
   int const nBlockCols ( multiplication.num_B_partition );
   int const nBlocks    ( nBlockRows*nBlockCols );
   
   size_t  C_ncols       (multiplication.C_ncols);
   PreProcess<T> const & preProcess ( * multiplication.preProcess );
   UniformPartition<int> P(0, nrows, nBlockRows);

   size_t atomicCounter  (multiplication.atomicCounter );    
     
   #pragma omp barrier
   if(tid == 0){
       atomicCounter = 0;
   }
   #pragma omp barrier

   int iBlock(0);
   bool done =false;

   while(!done) {

    
   #pragma omp barrier
   if(tid == 0){
       atomicCounter = 0;
   }
   #pragma omp barrier

   int iBlock(0);
   bool done =false;

   while(!done) {

   #pragma omp critical
   {
     iBlock = atomicCounter;
     atomicCounter += 1;
   }
     

   if( iBlock >= nBlocks )
   {
     done =true;
   } 
   else
   {
   
   int const block_idx_i (block_row_idx(iBlock,nBlockCols));
   int const block_idx_j (block_col_idx(iBlock,nBlockCols));

   int  row_offset = P.begin(block_idx_i+0);
   int  row_endset = P.begin(block_idx_i+1);

   size_t const * const  row_begin_ptr ( preProcess.row_begin_ptr(block_idx_j) );
   size_t const * const  row_end_ptr   ( preProcess.row_end_ptr(block_idx_j) );

  

   size_t  nnz_estimate
       ( SpGEMM_estimate_nnz
          ( A
          , B
          , row_offset
          , row_endset
          , row_begin_ptr
          , row_end_ptr 
          )
       );

   if (tid == 0){
   //  std::cout <<tid << " :  " << iBlock <<"  " << C_ncols <<  std::endl;
    }   

   // if (tid == 1)
    // std::cout <<tid << " :  " << iBlock << std::endl;

   //CSRMatric<T> & C_P1 ( C_P1[iBlock] );

   multiplication.C_P[iBlock].init(row_endset - row_offset, nnz_estimate, 0);
  // C_P1[iBlock].init(row_endset - row_offset, nnz_estimate, 0);
  
    if (tid == 0){
       std::cout <<tid << " :  " << iBlock <<"  " <<  C_P1[iBlock].get_nnzs() <<  std::endl;
    }


   SpGEMM_Sequential
      ( A
      , B
      , multiplication.C_P[iBlock]
      , row_offset
      , row_endset
      , row_begin_ptr
      , row_end_ptr
      , C_ncols
      );
  //   std::cout << "This is section " << std::endl; 
  
    #pragma omp barrier 
    {
     tid_nnz[block_idx_i] += multiplication.C_P[iBlock].RowPtr[ multiplication.C_P[iBlock].get_nrows() ];
    }
    
        
    }

   }
 

    CSRMatric<T> & C ( multiplication._C );

     #pragma omp barrier
   if( tid == 0 ){

    work_space[0] = 0;
    size_t s1 (0);

   for( int block_idx_i(1)
      ;     block_idx_i< nBlockRows+1
      ;     ++block_idx_i ){
      s1 += tid_nnz[block_idx_i-1];
      work_space[block_idx_i] = s1;
    }

    // std::cout <<"This is the nnzs : " << work_space[nBlockRows] << std::endl;
     multiplication._C.init( nrows, work_space[nBlockRows],0);
     multiplication._C.RowPtr[nrows] = work_space[nBlockRows];
   }



     #pragma omp barrier

  {

    UniformPartition<int> PBlockRows(0, nBlockRows, nThread );
    
    UniformPartition<int> P(0, nrows, nBlockRows );
  
    for(int block_idx_i(PBlockRows.begin(tid))
       ;    block_idx_i<PBlockRows.end(tid)
       ;    block_idx_i++)
     {
       int const row_offset ( P.begin(block_idx_i) );
       int const row_endset ( P.end  (block_idx_i) );

       int js = work_space[block_idx_i];
       for( int i(0)  
          ;     i<row_endset-row_offset
          ;     ++i)
        {
           C.RowPtr[i+row_offset] = js;

          for( int block_idx_j(0)
             ;     block_idx_j < nBlockCols
             ;     ++block_idx_j)      
           {
             int const iBlock ( block_idx(block_idx_i, block_idx_j, nBlockCols) );
           // CSRMatric<T> const & CP = multiplication.C_P[iBlock];

             for( int j= multiplication.C_P[iBlock].RowPtr[i+0]
                ;     j< multiplication.C_P[iBlock].RowPtr[i+1]
                ;     j++)
              {
                 C.ColIndices[js] =  multiplication.C_P[iBlock].ColIndices[j];          
                 C.Values[js]     =  multiplication.C_P[iBlock].Values[j];
                ++js;
              }
           }
        } 
     }
  }

 }
*/



#endif
