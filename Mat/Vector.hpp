#ifndef  VECTOR_H
#define  VECTOR_H

#include <vector>
#include <iostream>

template <typename T>
class Vector {

   public:
    Vector();
    
    Vector(int const _size);

    Vector(int const _size, T const _initial);
     
    void init(int const _size, T const _initial);

    void set(T const _value);
     //for Disply the vector content  
    void Display();

     //clear the vector entries    
    void clear();

    //Accessing the individual elem
  
     T const &  operator[](const int i) const ;

    T &  operator[](const int i);

    //= operator
     Vector<T> & operator=( Vector<T> const & rhs);

     Vector<T> const & operator=( Vector<T> const & rhs) const;

    //+ operator 
     Vector<T>  operator + ( Vector<T> const & rhs);

    //- operator 
     Vector<T>  operator - ( Vector<T> const & rhs);
  
     Vector<T>  operator - ( Vector<T> const & rhs) const;
  
    //* operator 
     T  operator * ( Vector<T> const & rhs);
  
      //Max norm
     T  MaxNorm ();
  
      //norm2
     T  Norm2 ();
  
     bool contain(T const element);

   //getting the size
     size_t Size() const;
  

     void resize( int const _size );

   protected:
   // T * vec;
     std::vector<T> vec;
    int size;
 };


#include "Vector.cpp"

#endif
