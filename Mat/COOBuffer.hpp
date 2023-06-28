#ifndef COOBUFFER_HPP
#define COOBUFFER_HPP

#include "Vector.hpp"
#include <algorithm>
 #include "csrmatric.hpp"


template <typename T>
class COOBuffer{
  public : 
  
  COOBuffer(){
       nnz = 0;
      }

  ~COOBuffer(){
       delete P;
      }

  COOBuffer(int _nnz){
       nnz = _nnz;
      }

  void setEntry(int row, int colIndices, T values){
        buffer.push_back(COO_T(row, colIndices, values));
      }

  void finalize(){ 

	std::sort
	(buffer.begin(),buffer.end(), 
		[](const COO_T & a, const COO_T & b) -> bool
		{ return a.row < b.row || (a.row == b.row && a.colIndices < b.colIndices); } 
	 );
	nnz = buffer.size();

	   // Constructing the CSR format based on the COO format
   P = new CSRMatric<T> ( buffer[buffer.size()-1].row + 1 , nnz, 0);
   
  for(int j=0;j< nnz;j++){
     P->ColIndices[j] = buffer[j].colIndices;
     P->Values[j]     = buffer[j].values;
     }

      int nz(0);
  for(int i(0); i< P->get_nrows(); i++){
       for(int j=nz;j<buffer.size();j++){
             if( i == buffer[j].row ) {
                 P->RowPtr[i]++;
                 nz++;
                 }               
           }
       }
 
  for(int i(0); i< P->get_nrows(); i++)
       P->RowPtr[i+1] += P->RowPtr[i] ;

  for(int i = P->get_nrows(); i > 0; i--)
       P->RowPtr[i] = P->RowPtr[i-1];

    
    P->RowPtr[0] = 0;  
  }


  CSRMatric<T> & CSRFormat(){
      return  * P; 
      }

   
    protected : 

 struct COO_T{
  COO_T(int _row, int _colIndices, T _values) : row(_row) , colIndices(_colIndices), values(_values){}

     int row;
     int colIndices;
     T values;
     };

     std::vector<COO_T> buffer;
     CSRMatric<T> * P; 
     int nnz;

 };
#endif
