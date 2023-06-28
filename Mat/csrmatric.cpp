#ifndef CSRMATRIC_CPP
#define CSRMATRIC_CPP
#include "csrmatric.hpp"

template<typename T>
CSRMatric<T> :: CSRMatric() {

      ColIndices.resize(0);
      RowPtr.resize(0);
      Values.resize(0);
    
      nnz=0;
      nrows=0;
      }

template<typename T>
CSRMatric<T> :: ~CSRMatric() {

      ColIndices.resize(0);
      RowPtr.resize(0);
      Values.resize(0);
    
      nnz=0;
      nrows=0;
      }

template<typename T>
CSRMatric<T> :: CSRMatric( int _nrows, size_t _nnz, T _initial ){
      
      nnz =_nnz;
      nrows = _nrows;

      Values.init(nnz, _initial );
      ColIndices.resize(nnz);
      RowPtr.resize(nrows + 1);
      }

template<typename T>
void CSRMatric<T> ::  init( int _nrows, size_t _nnz, T _initial ){
      
      nnz=_nnz;
      nrows=_nrows;
   

      Values.init(nnz,_initial);
      ColIndices.init(nnz, _initial);
      RowPtr.init(nrows + 1, _initial);
      }

template<typename T>
void CSRMatric<T> ::  Display(  ){
      
      std::cout << std::endl; std::cout << "Row Ptr" << std::endl;
     for(size_t i(0); i < this->get_nrows()+1; i++ ){
         std::cout << this->RowPtr[i] << "   ";
        }
      std::cout << std::endl; std::cout << "Col Indices" << std::endl;
      for(int j(0); j < this->get_nnzs(); j++ ){
        std::cout <<  this->ColIndices[j] << "   ";
        }
 
       std::cout << std::endl; std::cout << "Values " << std::endl;
      for(int j(0); j < this->get_nnzs(); j++ ){
        std::cout <<  Values[j] << "   ";
        }
       std::cout << std::endl;  
      }

template<typename T>
void CSRMatric<T> ::  Display() const {
      
      std::cout << std::endl; std::cout << "Row Ptr" << std::endl;
     for(size_t i(0); i < this->get_nrows()+1; i++ ){
         std::cout << this->RowPtr[i] << "   ";
        }
      std::cout << std::endl; std::cout << "Col Indices" << std::endl;
      for(int j(0); j < this->get_nnzs(); j++ ){
        std::cout <<  this->ColIndices[j] << "   ";
        }
 
       std::cout << std::endl; std::cout << "Values " << std::endl;
      for(int j(0); j < this->get_nnzs(); j++ ){
        std::cout <<  Values[j] << "   ";
        }
       std::cout << std::endl;  
      }


// = operator 
template<typename T>
CSRMatric<T> & CSRMatric<T> :: operator=(CSRMatric<T>  & rhs ){
      nrows = rhs.get_nrows();
      nnz = rhs.get_nnzs(); 

      for(int i(0); i < nrows+1; i++ ){
        this->RowPtr[i] = rhs.RowPtr[i];
        }
 
      for(int j(0); j < nnz; j++ ){
        this->ColIndices[j] = rhs.ColIndices[j];
        this->Values[j] = rhs.Values[j]; 
        }

       return *this;
       }

//matrix vector multiplication
template<typename T>
Vector<T>  CSRMatric<T> :: operator*(Vector<T>  & rhs ){
        
        Vector<T> result(this->nrows);  T sum =0;   

        for(int i(0); i < nrows; ++i){
           sum=0;
           for(int j=this->RowPtr[i]; j<this->RowPtr[i+1]; j++) {
                int k =  this->ColIndices[j];
                sum += this->Values[j]*rhs[k];
              }
            result[i] = sum;
        } 

        return result;
        }

template<typename T>
Vector<T> const  CSRMatric<T> :: operator*(Vector<T> const & rhs ) const {
        
        Vector<T> result(this->nrows);  T sum =0;   

        for(int i(0); i < nrows; ++i){
           sum=0;
           for(int j=this->RowPtr[i]; j<this->RowPtr[i+1]; j++) {
                int k =  this->ColIndices[j];
                sum += this->Values[j]*rhs[k];
              }
            result[i] = sum;
        } 

        return  result;
        }

template<typename T>
CSRMatric<T> CSRMatric<T> :: operator * (CSRMatric<T> & rhs ){
           
       CSRMatric<T> result;
       CSRMatric<T> result2;

       int ncols=0;
       for(int j=0; j< rhs.get_nnzs(); j++){
          if(rhs.ColIndices[j] >ncols ){
             ncols = rhs.ColIndices[j];
           }
       }

       size_t nnz_estimate(0);
       for(int irow(0)
          ;    irow<this->get_nrows()
          ;    irow++)
        {
           size_t row_nnz_estimate(0);
           for(int   icsr=this->RowPtr[irow+0]
               ;     icsr<this->RowPtr[irow+1]
               ;     icsr++)
            {
             int jcols = this->ColIndices[icsr];
             row_nnz_estimate += rhs.RowPtr[jcols+1]
                               - rhs.RowPtr[jcols+0]; 
            }
           nnz_estimate += row_nnz_estimate;
        }
        result.init(this->get_nrows(), nnz_estimate,0);

        Vector<T> X(ncols+1);
        int irow(0); int nnz(0); 
         T A_ik; T B_kj;


       for(int  i=0
           ;    i<this->get_nrows()
           ;    i++)
        {
     
         for(size_t iA = this->RowPtr[i+0]
            ;       iA < this->RowPtr[i+1]
            ;       iA++)
          {
            int jj = this->ColIndices[iA]; A_ik = this->Values[iA];

           for(size_t jA=rhs.RowPtr[jj+0]
              ;       jA<rhs.RowPtr[jj+1]
              ;       jA++)
            {          
                B_kj = rhs.Values[jA];
                X[ rhs.ColIndices[jA] ] += A_ik*B_kj ;
            }  
          }
       
         for(int it=0
            ;    it<X.Size()
            ;    it++)
           {
             if( X[it] != 0){
                 result.ColIndices[nnz]=it;
                 result.Values[nnz] = X[it];
                 ++nnz; 
               }
           }  
       
      	   X.clear();            
     	   ++irow;
       	   result.RowPtr[irow] = nnz;
        }

        result2.init(this->get_nrows(), nnz,0);
        
       for(int  i=0
           ;    i<result.get_nrows()
           ;    i++)
        {
          result2.RowPtr[i] = result.RowPtr[i]; 
         for(size_t iA = result.RowPtr[i+0]
            ;       iA < result.RowPtr[i+1]
            ;       iA++)
          {
            result2.ColIndices[iA] = result.ColIndices[iA];
            result2.Values[iA] = result.Values[iA]; 
          }

        }
       result2.RowPtr[result.get_nrows()] = result.RowPtr[result.get_nrows()];
       
      return result2;
}

//Transpsoe of CSR Matrix
template<typename T> 
CSRMatric<T> CSRMatric<T> :: Transpose(){

       CSRMatric<T> transpose(this->nrows, this->nnz,0);

       for (int i=0; i <= this->nrows; i++)
           transpose.RowPtr[i] = 0;
 
       for (int i=0; i < this->nnz; i++)
         	transpose.RowPtr[this->ColIndices[i]+1]++;

 
        
       for (int i=0; i < this->nrows; i++)
            transpose.RowPtr[i+1] += transpose.RowPtr[i];
   

       int ROW = this->nrows;

       for (int i=0; i < ROW ; i++){
          for (int j=this->RowPtr[i]; j < this->RowPtr[i+1]; j++){

    	    int index = this->ColIndices[j];

       	    transpose.ColIndices[transpose.RowPtr[index]] = i;

            transpose.Values[ transpose.RowPtr[index] ] = this->Values[ j ]; 

       	    transpose.RowPtr[index]++;

            }
       }     


 

       for(int i = ROW; i > 0; i--){
         transpose.RowPtr[i] = transpose.RowPtr[i-1];
       }
       transpose.RowPtr[0] = 0;




      return transpose;

        }

template<typename T>
size_t const & CSRMatric<T> ::  get_nrows() const{
        return this-> nrows; 
        }


template<typename T>
int const  & CSRMatric<T> ::  get_nnzs() const {
        return this-> nnz; 
        }



#endif
