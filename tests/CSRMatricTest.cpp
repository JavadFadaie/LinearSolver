#include <iostream>
#include "../Mat/csrmatric.hpp"
#include "../Mat/Vector.hpp"
#include <gtest/gtest.h>

int const nrows = 6;
size_t const nnz = 6;

class CSRMatrixInit{

    public: 
    CSRMatrixInit( double  & val){
      A.init(nrows,nnz,0);
  
      for(int i(0); i<nrows; i++){
        A.RowPtr[i] = i;
       }
       
       A.RowPtr[nrows] = nrows;
       
      for(int i(0); i<nnz; i++){
        A.ColIndices[i] = i;
        A.Values[i] = val;
        }
     } 
   
     CSRMatric<double> & GetMatric(){
       return A;  
     }
  
    CSRMatric<double> A;
};

TEST(CSRMatricTests, Construct)
{   

   CSRMatric<double> A(nrows, nnz, 0); 
   
   CSRMatric<double> B;
   B.init(nrows, nnz,0);
 
   ASSERT_EQ(nrows,A.get_nrows());
   ASSERT_EQ(nnz, A.get_nnzs());

   ASSERT_EQ(nrows, B.get_nrows());
   ASSERT_EQ(nnz, B.get_nnzs());
}

TEST(CSRMatricTests, EqualOperator)
{
    double val = 1;
    CSRMatrixInit Amatric(val);
    CSRMatric<double> A (Amatric.GetMatric());
    CSRMatric<double> B (nrows, nnz,0);
      
        B = A;
    
    for(int i(0); i<nrows+1; i++){    
      ASSERT_EQ(B.RowPtr[i], A.RowPtr[i]);
     }
     
    for(int i(0); i<nnz; i++){    
      ASSERT_EQ(B.ColIndices[i], A.ColIndices[i]);
     }
    
    for(int i(0); i<nnz; i++){    
      ASSERT_EQ(B.Values[i], A.Values[i]);
     }   
}

TEST(CSRMatricTests, SpGEMM)
{
  double A_val = 1;
  CSRMatrixInit AMatric(A_val);
  CSRMatric<double> & A (AMatric.GetMatric());
  
  double B_val = 2;
  CSRMatrixInit BMatric(B_val);
  CSRMatric<double> & B (BMatric.GetMatric());
  
  CSRMatric<double>  C ( B*A );
  
   for(int i(0); i<nrows+1; i++){    
      ASSERT_EQ( C.RowPtr[i], B.RowPtr[i]);
     }
     
    for(int i(0); i<nnz; i++){    
      ASSERT_EQ( C.ColIndices[i], B.ColIndices[i]);
     }
    
    for(int i(0); i<nnz; i++){    
      ASSERT_EQ( C.Values[i], B.Values[i]);
     }   
     
}

TEST(CSRMatricTests, MatricVector)
{
  double A_val = 1;
  CSRMatrixInit AMatric(A_val);
  CSRMatric<double>  A (AMatric.GetMatric());
 
  Vector<double> V(nrows,2);
  
  Vector<double> const & RV((A*V));
  
   for(int i(0); i<nrows; i++){    
      ASSERT_EQ(RV[i],2);
     }
     
}

int CSR_Matric(int argc, char **argv){
   testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();

}
