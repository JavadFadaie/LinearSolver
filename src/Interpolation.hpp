#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <iostream>

template<typename T>
class Interpolation{

   public : 

   Interpolation
       ( CSRMatric<T> const & A
       , CSRMatric<T> const & S
       , Vector<int> & CF_marker
       , int const & nThread );



   void 
     apply();     

   CSRMatric<T> const & 
      Get_Interpolation(); 
 
   int &
     NumCoarse();

   private:

    static void 
      Interpolation_thread_Func(int tid, void * arg);

    static void
      Interpolation_Structure
            ( int tid
            , CSRMatric<T> const & S
	    , Vector<int>  & CF_marker
            , std::vector<int>  & fine_to_coarse
            , std::vector<int>  & coarse_counter
	    , std::vector<int> & nnz_interp
            , int Offset  
            , int Endset);

    static void 
      Interpolation_Weigth_Seq
            ( CSRMatric<T> const & S
            , CSRMatric<T> const & A
            , CSRMatric<T> & P
            , std::vector<int> & fine_to_coarse
            , Vector<int> & CF_marker
            , size_t jj_counter_1
            , int Offset
            , int Endset );


     CSRMatric<T> const & _A;
     CSRMatric<T> const & _S;
     CSRMatric<T>   _P;
 
     Vector<int>  & _CF_marker;
   

     std::vector<int> _fine_to_coarse;
     std::vector<int> _coarse_shift;
     std::vector<int> _coarse_counter;
     std::vector<int> _nnz_interp;
     std::vector<int> _nnz_thread;
     int const & _nThread;    

};

////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////


template<typename T>
Interpolation<T> 
       :: Interpolation
           ( CSRMatric<T> const & A
           , CSRMatric<T> const & S
           , Vector<int> & CF_marker
           , int const & nThread )
:_A(A)
,_S(S)
,_CF_marker(CF_marker)
,_nThread(nThread)
{}

template<typename T>
void 
Interpolation<T>::apply()
{
   _coarse_counter.resize(_nThread,0);
   _nnz_interp.resize(_nThread,0);
   _fine_to_coarse.resize(_A.get_nrows(),-1);
   _coarse_shift.resize(_nThread+1,0);
   _nnz_thread.resize(_nThread+1,0);
 
   omp_set_num_threads(_nThread);

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    
    Interpolation_thread_Func(tid , this);
  }
} 

template<typename T>
int & 
Interpolation<T>::NumCoarse()
{
  return _coarse_counter[_nThread-1];
}

template<typename T>
CSRMatric<T> const &
Interpolation<T> :: Get_Interpolation(){

    return _P;
}

template<typename T>
void 
Interpolation<T>::Interpolation_thread_Func(int tid, void *arg)
{

   Interpolation<T> & interpolation(* reinterpret_cast <Interpolation<T> *> (arg) );

   CSRMatric<T> const & A              (interpolation._A); 
   CSRMatric<T> const & S              (interpolation._S);
   Vector<int>  & CF_marker      (interpolation._CF_marker);
   int const & numThreads          (interpolation._nThread); 
  
   std::vector<int> & coarse_counter   (interpolation._coarse_counter);
   std::vector<int> & fine_to_coarse   (interpolation._fine_to_coarse);
   std::vector<int> & nnz_interp       (interpolation._nnz_interp);
             
   UniformPartition<int> Q(0, A.get_nrows(), numThreads);

   int Offset = Q.begin(tid+0);
   int Endset = Q.begin(tid+1);


 
 
   Interpolation_Structure
       ( tid
       , S
       , CF_marker
       , fine_to_coarse
       , coarse_counter
       , nnz_interp
       , Offset  
       , Endset );

     std::vector<int> & coarse_shift  (interpolation._coarse_shift); 
     std::vector<int> & nnz_thread    (interpolation._nnz_thread); 

      // Allocate the Interpolation Matrix
      #pragma omp barrier

      if(tid == 0){
       
       for(int i(0); i < numThreads-1 ; ++i){
          coarse_counter[i+1] += coarse_counter[i];
          nnz_interp[i+1]     += nnz_interp[i];
          coarse_shift[i+1]    = coarse_counter[i]; 
          nnz_thread[i+1]      = nnz_interp[i];
        }

        int i = numThreads - 1;
        size_t P_diag_size = nnz_interp[i];
        coarse_shift[numThreads] = coarse_counter[i];
        nnz_thread[numThreads]  = nnz_interp[i];

     
        interpolation._P.init( A.get_nrows(),P_diag_size , 0 );
        interpolation._P.RowPtr[A.get_nrows()] = P_diag_size;
     
      }

       #pragma omp barrier

       size_t shift = coarse_shift[tid];
       for (int i = Offset; i < Endset; i++){
        if( fine_to_coarse[i] != F_PT ){
           fine_to_coarse[i] += shift;
         }
       }

       #pragma omp barrier

      CSRMatric<T>  & P ( interpolation._P );
      size_t jj_counter=0;
      jj_counter = nnz_thread[tid];



     Interpolation_Weigth_Seq
            ( S
            , A 
            , interpolation._P
            , fine_to_coarse
            , CF_marker
            , jj_counter
            , Offset
            , Endset );


}


template<typename T>
void 
Interpolation<T>
        :: Interpolation_Weigth_Seq
            ( CSRMatric<T> const & S
            , CSRMatric<T> const & A
            , CSRMatric<T> & P
            , std::vector<int> & fine_to_coarse
            , Vector<int> & CF_marker
            , size_t jj_counter_1
            , int Offset
            , int Endset )
{

  long  jj_counter(jj_counter_1);
  int nRow = A.get_nrows();
  std::vector<int> P_marker;
  P_marker.resize(nRow,-1);
  long strong_fine_point = -2;

 for(int i=Offset
    ;    i<Endset
    ;    i++)
{


   if(CF_marker[i] == C_PT ){
     P.RowPtr[i] = jj_counter;
     P.ColIndices[jj_counter] = fine_to_coarse[i];
     P.Values[jj_counter] = 1;
     jj_counter++; 
     }
    else
     {
   
      P.RowPtr[i] = jj_counter;
      int jj_begin_row = jj_counter;
  
         for(size_t  jj=S.RowPtr[i+0]
            ;        jj<S.RowPtr[i+1]
            ;        ++jj)
            {
              int i1 = S.ColIndices[jj];

                if(CF_marker[i1] == C_PT){
                   P_marker[i1] = jj_counter;
                   P.ColIndices[jj_counter] = fine_to_coarse[i1];
                   P.Values[jj_counter] = 0;
                   jj_counter++; 
                   }
                 else if(CF_marker[i1] != F_ST_PT )
                  {
                   P_marker[i1] = strong_fine_point;        
                  }
            }

      
       int jj_end_row = jj_counter;
       T diagonal(0);
       for(size_t jj=A.RowPtr[i+0] 
           ;       jj<A.RowPtr[i+1]
           ;       ++jj)   
         {              
            
               int i1 = A.ColIndices[jj];
   
               // Finding the diagonal elements
                if(i1 == i){
                   diagonal += A.Values[jj];
                   continue;
                  }

               // i1 is C point and strongly influence i which is F_PT
                 if(P_marker[i1] >= jj_begin_row){
                     P.Values[ P_marker[i1] ] += A.Values[jj];
                   } 
               

                else if(P_marker[i1] == strong_fine_point ){
                    T sum(0);
                    for(size_t jj1=A.RowPtr[i1+0]
                        ;      jj1<A.RowPtr[i1+1]
                        ;      ++jj1){
                       
                         int i2 = A.ColIndices[jj1]; 
                         if( P_marker[i2] >= jj_begin_row ){
                             sum += A.Values[jj1];
                         }

                       }


                   if(sum != 0){
                     
                       T distribute = A.Values[jj]/sum;
                       for(size_t  jj1=A.RowPtr[i1+0]
                          ;        jj1<A.RowPtr[i1+1]
                          ;        ++jj1){
                       
                         int i2 = A.ColIndices[jj1]; 
                         if( P_marker[i2] >= jj_begin_row ){
                             P.Values[ P_marker[i2] ] += distribute*A.Values[jj1];
                             }

                          }
                    }else
                        {
                           diagonal += A.Values[jj]; 
                        } 

                   }
                 else if ( CF_marker[i1] != F_ST_PT  ) {
                    diagonal += A.Values[jj];
                  }

             }

          
            if(diagonal == 0.0){
              for(size_t jj=jj_begin_row
                 ;       jj<jj_end_row
                 ;       ++jj)
               {     
                  P.Values[jj] =0.0;
               } 
            }
            else{
              for(size_t jj=jj_begin_row
                 ;       jj<jj_end_row
                 ;       ++jj)  
               {   
                  P.Values[jj] /= -diagonal;
               }
            }
         }
     
      strong_fine_point--;
      }


}     


template<typename T>
void 
Interpolation<T> 
          ::  Interpolation_Structure
            ( int tid
            , CSRMatric<T> const & S
            , Vector<int>  & CF_marker
            , std::vector<int>  & fine_to_coarse
            , std::vector<int>  & coarse_counter
            , std::vector<int> & nnz_interp
            , int Offset
            , int Endset){

   
   for(int i= Offset
      ;    i< Endset
      ;    i++)
     {
     
      if( CF_marker[i] == C_PT ){
           nnz_interp[tid]++;
           fine_to_coarse[i] =  coarse_counter[tid];
          coarse_counter[tid]++;   
         }
         else
          {         
             for(size_t jj=S.RowPtr[i+0]
               ;        jj<S.RowPtr[i+1]
               ;        ++jj
               )
              {   
                int i1 = S.ColIndices[jj];
                if( CF_marker[i1] == C_PT ){
                    nnz_interp[tid]++;
                  }
               }
           }

     }
  

}







#endif
