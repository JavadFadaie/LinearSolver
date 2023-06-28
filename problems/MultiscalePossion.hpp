#ifndef MULTISCALEPOISSION_HPP
#define MULTISCALEPOISSION_HPP

#include "../Mat/COOBuffer.hpp"
#include "../Mat/csrmatric.hpp"
#include "problem.hpp"
 
template<typename T>
class MultiscalePoission : public Problem<T> {
   
public :

 MultiscalePoission(int _n1d, double _kI)
   :Problem<T>(_n1d)
   ,n1d(_n1d)
   ,kI(_kI)     
   {}  
 
private :    
 void
 setMatric( COOBuffer<T> & A  ) const override {
  
   double x, y;
   double h = 1.0/(n1d+1 );
   
   //double kI = 1.0;
   double kO = 1.0;
   double kapa = (2*kI*kO)/(kI+kO);

   double left,right,up, down, center;
   std::cout <<"This is kapa  " <<kapa << std::endl;

  for(int k = 0; k < n1d*n1d ; ++k){

     	 x = (k%n1d + 1)*h;
         y = (k%(n1d*n1d)/n1d + 1)*h;
  	 
       //   std::cout<< "(" << x <<", " << y << ")" << std::endl;
 
         // The Halo region distritization 
         if( k%n1d == 0 ){
          std::cout << std::endl;
          std::cout << std::endl;
         }

        std::cout<< "(" << x <<", " << y << ")" << std::endl;

         if( (y >= 0.25-h) & (y<= 0.75+h) & (x >= 0.25-h) & (x<=0.75+h )  )         
          {
           
              if( (y > 0.25) & (y< 0.75) & (x > 0.25) & (x<0.75)  )         
              {
                  //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kI;
                    right = kI;
                    up = kI;
                    down = kI;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
		      }
              else
              {
               
               // std::cout<< "(" << x <<", " << y << ")" << "  ";

                     //the halo horizental below line
                if( (y == 0.25-h) &  ( x > 0.25-h ) & ( x < 0.75+h )  )
                 {
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kO;
                    right = kO;
                    down = kO;
                    up = kapa;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
		     
                 }

                    //the halo vertical left below line
                if( (x == 0.25-h) &  ( ( y > 0.25-h ) & ( y < 0.75+h ) )  )
                 {
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kO;
                    right = kapa;
                    down = kO;
                    up = kO;
                    center = left + right + up + down;
                    
                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
                 }

                    //the halo vertical right below line
                if( (x == 0.75+h) &  ( ( y > 0.25-h ) & ( y < 0.75+h ) )  )
                 {
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kapa;
                    right = kO;
                    down = kO;
                    up = kO;
                    center = left + right + up + down;
                  
                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
                                    
                 }

                     //the halo horizental upper line
                if( (y == 0.75+h) &  ( x > 0.25-h ) & ( x < 0.75+h )  )
                 {
                   // std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kO;
                    right = kO;
                    down = kapa;
                    up = kO;
                    center = left + right + up + down;
                
                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
                   
                } 


                       // The Boundary region distritization

                    // the horizental below line
                if( (y == 0.25) & ( ( x > 0.25 ) & ( x < 0.75 ) ) )
                { 
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kI;
                    right = kI;
                    down = kapa;
                    up = kI;
                    center = left + right + up + down;
                
                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
                     
                }

                      // the vertical left line
                if( (x == 0.25) & ( ( y > 0.25 ) & ( y < 0.75 ) ) )
                { 
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kapa;
                    right = kI;
                    down = kI;
                    up = kI;
                    center = left + right + up + down;
                
                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
                      
                }

                      // the vertical right line
                if( (x == 0.75) & ( ( y > 0.25 ) & ( y < 0.75 ) ) )
                { 
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kI;
                    right = kapa;
                    down = kI;
                    up = kI;
                    center = left + right + up + down;
                
                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
                
                }

                    // the horizental upper line
                if( (y == 0.75) & ( ( x > 0.25 ) & ( x < 0.75 ) ) )
                { 
                    //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kI;
                    right = kI;
                    down = kI;
                    up = kapa;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
        
                }
                
                        //the below left corner 
               if( (x == 0.25) &  ( y == 0.25 ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kapa;
                    right = kI;
                    down = kapa;
                    up = kI;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
           
               }
                      //the below right corner
               if( (x == 0.75) &  ( y == 0.25 ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kI;
                    right = kapa;
                    down = kapa;
                    up = kI;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
               }

                        //the upper left corner
               if( (x == 0.25) &  ( y == 0.75 ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kapa;
                    right = kI;
                    down = kI;
                    up = kapa;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
               }

                         //the upper right corner
               if( (x == 0.75) &  ( y == 0.75 ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kI;
                    right = kapa;
                    down = kI;
                    up = kapa;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);
               }




                         //The corner

               if( (x == 0.25-h) &  ( y == 0.25-h ) )
               {
                    // std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kO;
                    right = kO;
                    down = kO;
                    up = kO;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);   
               }
               if( (x == 0.75+h) &  ( y == 0.25-h ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kO;
                    right = kO;
                    down = kO;
                    up = kO;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);   
               }
               if( (x == 0.25-h) &  ( y == 0.75+h ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                    left = kO;
                    right = kO;
                    down = kO;
                    up = kO;
                    center = left + right + up + down;

                    A.setEntry(k, k, -center);
    				A.setEntry(k, k - 1, left);
			     	A.setEntry(k, k + 1 , right);
					A.setEntry(k, k - n1d, down);
					A.setEntry(k, k + n1d , up);   
               }
               if( (x == 0.75+h) &  ( y == 0.75+h ) )
               {
                     //std::cout<< "(" << x <<", " << y << ")" << "  ";
                   left = kO;
                   right = kO;
                   down = kO;
                   up = kO;
                   center = left + right + up + down;

                   A.setEntry(k, k, -center);
    			   A.setEntry(k, k - 1, left);
			       A.setEntry(k, k + 1 , right);
				   A.setEntry(k, k - n1d, down);
				   A.setEntry(k, k + n1d , up);   
               }
    

              } 
          }
          else
          {
                 // This is Outer Box
              //std::cout<< "(" << x <<", " << y << ")" << "  ";
           
                A.setEntry(k, k, -4.0);
    			if(k%n1d > 0){
			  	A.setEntry(k, k - 1, kO);
				}
				if(k%n1d < n1d - 1) {
			  	A.setEntry(k, k + 1 , kO);
				}
		   	    if(k%(n1d*n1d) >= n1d) {
				A.setEntry(k, k - n1d, kO);
				}
				if(k%(n1d*n1d) < n1d*n1d - n1d) {
				A.setEntry(k, k + n1d , kO);
				} 
         }
       
       
     }


       A.finalize();       
     }


   
 void
  setRhsVector
    ( Vector<T> & _rhs
    , Vector<T> & _exactSol ) const override {
   
    double x, y; double phi;
    double sol;

    double h = 1.0/(n1d + 1);
    double h2inv = (h*h);
      std::cout << kI << "  " << 1/kI << std::endl;
	for(int k=0; k <  n1d*n1d; k++){  
    	 x = (k%n1d + 1)*h;
   	     y = (k%(n1d*n1d)/n1d + 1)*h;
      
      	phi = ( -2*(M_PI)*(M_PI) )*sin(M_PI*x)*sin(M_PI*y);
      	_rhs[k] = h2inv*phi;
    
             // The inner box region
           sol = sin(M_PI*x)*sin(M_PI*y);
      	 
       if( (y > 0.25) & (y< 0.75) & (x > 0.25) & (x< 0.75 )   )         
        {
           _exactSol[k] = (1/kI)*sol; 
     
        }else
        {
           _exactSol[k] = sol; 
        }
      }

     }

  private :
      int  n1d;
      double kI;
 
 }; 


#endif
