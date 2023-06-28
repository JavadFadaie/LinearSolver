#ifndef CSRMATRIC_HPP
#define CSRMATRIC_HPP
#include "Vector.hpp"
//#include "COOBuffer.hpp"

template<typename T>
class CSRMatric{

   public: 

   CSRMatric();
   virtual ~CSRMatric();

   CSRMatric( int _nrows, size_t _nnz, T _initial );
   void init( int _nrows, size_t _nnz, T _initial );
  
   void Display();
   
   void Display() const;
   // Constructor based on the COO Buffer
  // void COOBuffer(std::vector< COO_T > _buffer);
    
   // = operator   
   CSRMatric<T> & operator = (CSRMatric<T>  & rhs );

   //matrix vector multiplication
   Vector<T>  operator * ( Vector<T>  & rhs );  
 
   Vector<T> const operator * ( Vector<T> const & rhs ) const;  
 
   //matrix matrix multiplciation
   CSRMatric<T>  operator * ( CSRMatric<T> & rhs );

  
   //CSR Matrix Transpose 
   CSRMatric<T>  Transpose();

    
   size_t const & get_nrows() const ;
   int const  & get_nnzs() const ; 

  
   Vector<int> ColIndices;
   Vector<size_t> RowPtr;
   Vector<T> Values; 

   protected:

   int nnz;
   size_t nrows; 
      

 };

#include "csrmatric.cpp"
#endif
