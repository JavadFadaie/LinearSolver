#ifndef RUGECOARSEN_HPP
#define RUGECOARSEN_HPP

#include <iostream>
#include "MeasureList.hpp"
#include <memory>
#include "Transpose.hpp"

int UNDINITIALIZED = -2;
int UNDECIDED = 0;


int C_PT = 1;
int F_PT = -1;
int F_ST_PT = -3;


template<typename T>
class RugeCoarsen{

  public :

       RugeCoarsen(CSRMatric<T> const & StrengthMatric, int const & nThread );

       void apply();

       void Second_Phase_Ruge();

       Vector<int> & CF_Marker();
  
    
  protected :
    static
    void Coarsening_Function(void *arg);

    static
    void Ruge_Stuben
           ( CSRMatric<T> const & S
           , CSRMatric<T> const & S_T
           , Vector<int> & CF_Marker
           , Vector<size_t> & measure_array   
           );

    static
    void Ruge_second_pass(void *arg);
     
    CSRMatric<T> const & _S;
    int const & _nThread;
    std::unique_ptr< Transpose <T> > trans_S;
    Vector<int> _CF_Marker;
    Vector<size_t> _measure_array;
    int nRows;
};



////////////////////////////////////////////////////////////////////////////////
//                                 implementation
////////////////////////////////////////////////////////////////////////////////



template<typename T>
RugeCoarsen<T> 
        :: RugeCoarsen( CSRMatric<T> const & StrengthMatric, int const  & nThread )
:_S( StrengthMatric )
,_nThread(nThread)
,trans_S(nullptr)
,_CF_Marker(0)
,_measure_array(0)
,nRows(0)
{
 
    trans_S = std::unique_ptr< Transpose <T> >
                    (new Transpose <T>( StrengthMatric
                                       , StrengthMatric.get_nrows()
                                       , _nThread ) );
    trans_S->apply();        
 
       
}


template<typename T>
void RugeCoarsen<T>
            ::apply(){
    
   nRows = _S.get_nrows();
   _CF_Marker.init(nRows, UNDINITIALIZED );
   _measure_array.init(nRows, UNDINITIALIZED );
    Coarsening_Function(this);
     }

template<typename T>
void RugeCoarsen<T> 
            :: Second_Phase_Ruge(){
   
   Ruge_second_pass(this);
}

template<typename T>
Vector<int> &
 RugeCoarsen<T> :: CF_Marker(){

   return _CF_Marker;
       
}

template<typename T>
void RugeCoarsen<T> :: Coarsening_Function(void *arg){
     
    RugeCoarsen<T> & rugeCoarsen(* reinterpret_cast<RugeCoarsen<T>*> (arg));

    CSRMatric<T>   const & S             (rugeCoarsen._S);
    CSRMatric<T>   const & S_T           (rugeCoarsen.trans_S->Get_Transpose());
    Vector<int>     & CF_Marker     (rugeCoarsen._CF_Marker);
    Vector<size_t>  & measure_array (rugeCoarsen._measure_array);
 
    Ruge_Stuben
           (  S
           ,  S_T
           ,  CF_Marker
           ,  measure_array   
           );

     }

template<typename T>
void RugeCoarsen<T> :: Ruge_Stuben
           ( CSRMatric<T> const & S
           , CSRMatric<T> const & S_T
           , Vector<int> & CF_Marker
           , Vector<size_t> & measure_array   
           ){


     size_t num_left=0; 
     int graph_index=0;  int Coarse_index=0;
     int nabor = 0;      int nabor_two=0;
     int num_variables=S.get_nrows();

     MeasureList<int> measureList(num_variables);
  
     for(int i=0
        ;    i<num_variables
        ;    ++i)
        {

         measure_array[i] = S_T.RowPtr[i+1]
                          - S_T.RowPtr[i+0];

         size_t measure = measure_array[i];

         if( ( S.RowPtr[i+1]-S.RowPtr[i+0]) == 0)
          {
           CF_Marker[i] = F_ST_PT;
           measure_array[i]=0;
          }
          else
          {
            CF_Marker[i] = UNDECIDED;
            num_left++;
          }
        }
        
      for(int i(0)
         ;    i<num_variables
         ;    ++i)
       {
         size_t const measure(measure_array[i] );
 
         if(CF_Marker[i] != F_ST_PT)
         {
           if(measure >0){
               measureList.insert(i,measure);
            }
           else
            {
             if(measure < 0) {
               std::cout <<"measure 0 detected " 
                         << measure<< std::endl;
              }
              
              CF_Marker[i] = F_PT;
              for( size_t j=S.RowPtr[i+0]
                 ;        j<S.RowPtr[i+1]
                 ;      ++j)
               {
                 int const neighbor (S.ColIndices[j]);
                 if(CF_Marker[neighbor] !=F_ST_PT)
                 { 
                    if(neighbor < i)
                    {
                     size_t neighbor_measure(measure_array[neighbor]);
                     if(neighbor_measure > 0){
                       measureList.erase(neighbor,neighbor_measure);
                      }
                      measure_array[neighbor]=++neighbor_measure;
                      measureList.insert(neighbor,neighbor_measure); 
                    }
                    else
                    {
                      ++(measure_array[neighbor]);
                    }
                 }
               }
               --num_left;
             } 
          }
       }  

      while( num_left > 0 )
       {
          
           int const coarseIndex(measureList.popIndexWithMaxMeasure());
             
           CF_Marker[coarseIndex] = C_PT;
           measure_array[coarseIndex]=0;
           --num_left;

           for(size_t j=S_T.RowPtr[coarseIndex+0]
               ;      j<S_T.RowPtr[coarseIndex+1]
               ;      ++j)   
              {
                 
                int const  neighbor (S_T.ColIndices[j]);
                if( CF_Marker[neighbor] == UNDECIDED )
                 {
                  CF_Marker[neighbor] = F_PT;
                  measureList.erase(neighbor,measure_array[neighbor]);
                  --num_left;  

                  for(size_t k=S.RowPtr[neighbor+0]
                      ;      k<S.RowPtr[neighbor+1]
                      ;      k++)
                   {
                            
                      int const neighbor_two (S_T.ColIndices[k]);
                      if( CF_Marker[neighbor_two] == UNDECIDED)
                       {
                          size_t measure(measure_array[neighbor_two]);
                          measureList.erase(neighbor_two,measure);
                          measure_array[neighbor_two] = ++measure;
                          measureList.insert(neighbor_two,measure);
                       }
                    }
                 }    
              }
             for(size_t j=S.RowPtr[coarseIndex+0]
               ;        j<S.RowPtr[coarseIndex+1]
               ;        ++j)   
               {
                 int const neighbor (S.ColIndices[j]) ;
                  if(  CF_Marker[neighbor] == UNDECIDED )
                   {
                     size_t measure(measure_array[neighbor]);
                     measureList.erase(neighbor,measure);
                     measure_array[neighbor]=--measure;
                     
                     if(measure >0){
                        measureList.insert(neighbor,measure);
                      }
                      else
                      {
                       CF_Marker[neighbor]=F_PT;
                       --num_left;

                       for(size_t k = S.RowPtr[neighbor+0]
                          ;       k < S.RowPtr[neighbor+1]
                          ;      ++k)
                        {
                          size_t const neighbor_two(S.ColIndices[k]);
                           if(CF_Marker[neighbor_two] == UNDECIDED)
                           {
                            size_t neighbor_measure(measure_array[neighbor_two]);
                            measureList.erase(neighbor_two,neighbor_measure);
                            measure_array[neighbor_two] = ++neighbor_measure;
                            measureList.insert(neighbor_two,neighbor_measure);
                           }
                         }             
                       } 
                    }
                 }
            }//end while


   }

template<typename T>
void
RugeCoarsen<T>::Ruge_second_pass(void *arg)
{
 RugeCoarsen<T> & rugeCoarsen ( * reinterpret_cast<RugeCoarsen<T>* >(arg) );

   CSRMatric<T>   const & S        (rugeCoarsen._S);
   Vector<int>     & CF_Marker     (rugeCoarsen._CF_Marker);
   
   Vector<int>     graph_array(S.get_nrows());
   
   int nrows=S.get_nrows(); 
   int index=0; int nabor=0;
   int index_two; int set_empty;

    for(int   i=0
       ;      i<nrows
       ;      i++){

       if( CF_Marker[index] ==F_PT ){

          for(size_t  ji=S.RowPtr[index+0]
             ;        ji<S.RowPtr[index+1]
             ;         ++ji){
              nabor = S.ColIndices[ji];
              if(CF_Marker[nabor] == C_PT){
                  graph_array[nabor] = index;  
               }
            }


           for(size_t   ji=S.RowPtr[index+0]
              ;         ji<S.RowPtr[index+1]
              ;         ++ji ){
              
              nabor = S.ColIndices[ji];
              set_empty=1;
              if( CF_Marker[nabor] == F_PT ){
                 set_empty = 0;
                 for(size_t  jj=S.RowPtr[nabor+0]
                    ;        jj<S.RowPtr[nabor+1]
                    ;        ++jj){

                    index_two = S.ColIndices[jj];
                    if( graph_array[index_two] == index){
                       set_empty = 1;
                       }
                    }
                    if(set_empty == 0){
                       CF_Marker[nabor] = C_PT;
                     }
                 } 
               }
            }

          for(size_t ji=S.RowPtr[index+0]
             ;       ji<S.RowPtr[index+1]
             ;       ++ji)
           {
               nabor = S.ColIndices[ji];
               if( CF_Marker[nabor] == C_PT ){
                   graph_array[nabor] = -1;
                } 
              }  
           } 

     }



 
       


#endif
