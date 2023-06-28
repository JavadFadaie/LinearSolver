#include<iostream>
#include "../Mat/Vector.hpp"
#include <gtest/gtest.h>

int const n  = 6;

TEST(VectorTests, Construct)
{   

   Vector<int> vec(n);
   
   ASSERT_EQ(n,vec.Size());

   Vector<int> const vec2(n);
   ASSERT_EQ(n,vec2.Size());

}

TEST(VectorTests, Assign )
{
   Vector<int>  vec(n,n);
   ASSERT_EQ(n,vec[0]);
   
   
   Vector<int> const vec2(n,n);
   ASSERT_EQ(n,vec[0]);

}

TEST(VectorTests, set)
{
  Vector<int> vec(n);
  vec.set(n);
  
  ASSERT_EQ(vec[0],n);

}

TEST(VectorTests, equalOperator)
{
  Vector<int> vec(n,n);
  Vector<int> vec1(n);
  
  vec1 = vec;
  
  ASSERT_EQ(vec1.Size(),n);
  ASSERT_EQ(vec1[0],n);

}

TEST(VectorTests, minusOperator)
{
  Vector<int> vec(n,n);
  Vector<int> vec2(n,n+1);
  Vector<int> vec3(n);
 
  vec3 = vec2-vec;
  
  ASSERT_EQ(vec3[0],1);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
