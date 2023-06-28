#ifndef  VECTOR_CPP
#define  VECTOR_CPP
#include "Vector.hpp"
#include <cmath>  

template <typename T>
Vector<T> :: Vector(){
     }

template <typename T>
Vector<T> :: Vector(int const _size){
     size = _size;
     vec.resize(size);

     }

template <typename T>
Vector<T> :: Vector(int const _size, T const _initial){
     size = _size;
     vec.resize(size);
     
     for(int i=0; i < size; i++)
         vec[i] = _initial; 
  
     }

template <typename T>
void Vector<T> :: init(int const _size, T const _initial){
     size = _size;
     vec.resize(size);
     
     for(int i=0; i < size; i++)
         vec[i] = _initial; 
  
     }

template<typename T>
void Vector<T> :: set(T const _value){

     for(int i(0); i<this->size; i++){
        this->vec[i]= _value;
     } 

}

template<typename T>
void Vector<T> ::  Display(){
      
      std::cout << std::endl; std::cout << "Vector Content" << std::endl;
     for(size_t i(0); i < this->size; i++ ){
         std::cout << i << " : " << this->vec[i] << std::endl;
        }
}


template<typename T>
void Vector<T> :: clear(){

     for(int i=0; i<this->size; i++){
         this->vec[i]=0;
      }
}

//= operator
// in arg we have to add const because the vlue in receive might 
template <typename T>
Vector<T> & Vector<T> :: operator = ( Vector<T> const & rhs){
   
    size = rhs.Size();

    for(int i=0; i<size; i++){
       this->vec[i] = rhs.vec[i];
     }


    return *this;
    }

template <typename T>
Vector<T> const & Vector<T> :: operator = ( Vector<T> const & rhs) const {
   
    size = rhs.Size();

    for(int i=0; i<size; i++){
       this->vec[i] = rhs.vec[i];
     }


    return *this;
    }

// + operator
template <typename T>
Vector<T>  Vector<T> :: operator+(Vector<T> const & rhs){
    if(this->size != rhs.Size())
          std::cout<< "Run time error" << std::endl;
          
        
    Vector result(rhs.Size());

    for(int i=0; i< this->size; i++)
       result[i] = this->vec[i] + rhs.vec[i];
       
  
    return result;
    }

//- operator
template <typename T>
Vector<T>  Vector<T> :: operator-( Vector<T> const & rhs){
    if(this->size != rhs.Size())
          std::cout<< "+ Run time error" << std::endl;
        
    Vector result(rhs.Size());

    for(int i=0; i< this->size; i++)
       result[i] = this->vec[i] - rhs.vec[i];
       
  
    return result;
    }

template <typename T>
Vector<T>  Vector<T> :: operator-( Vector<T> const & rhs) const {
    if(this->size != rhs.Size())
          std::cout<< "+ Run time error" << std::endl;
        
    Vector result(rhs.Size());

    for(int i=0; i< this->size; i++)
       result[i] = this->vec[i] - rhs.vec[i];
       
  
    return result;
    }

//dot product
template <typename T>
T  Vector<T> :: operator * ( Vector<T> const & rhs){
    if(this->size != rhs.Size())
        std::cout<< "Run time error" << std::endl;
          
        
    T result = 0;
    for(int i=0; i< this->size; i++)
       result += (this->vec[i]) * (rhs.vec[i]);
  
 return result;
    }

//MaxNorm of vector
template <typename T>
T  Vector<T> :: MaxNorm (){
        
       T result(0);
      
    for(int i=0; i< this->size; i++)
       if( std::abs(this->vec[i]) >= result ){
           result = std::abs( this->vec[i]);
          }
       
  
    return result;
    }

//Norm2 of vector
template <typename T>
T  Vector<T> :: Norm2 (){
        
       T result(0);
       //  result(0);

    for(int i=0; i< this->size; i++){
         result += (this->vec[i])*(this->vec[i]);
      }
       
  
    return sqrt(result);
    }

// Accessing operator
template <typename T>
T const & Vector<T> :: operator[](const int i) const {
    return this-> vec[i];
    }


template <typename T>
T & Vector<T> :: operator[](const int i){
    return vec[i];
    }

// Size of the vector
template <typename T>
size_t Vector<T> :: Size() const {
    return this->size;
    }

template <typename T>
void Vector<T> :: resize( int _size ) {
    
     size = _size;
     vec.resize(size);

    }
#endif
