#ifndef ITERATIVESOLVERTEST_HPP
#define ITERATIVESOLVERTEST_HPP

#include "../problems/csrstencil.hpp"
#include "../iterativeSolver/GuassSiedelSolver.hpp"

class IterativeSolverTest{

   public : 

    IterativeSolverTest()
    {}

    double generateProblemAndSolve(int const n1d){

          
      pProblem =  std::unique_ptr< PoissonBVP2D<double> >
                        ( new PoissonBVP2D<double> (n1d) ); 
 
      pProblem->generateDiscretization();
      pProblem->randomApproxVectorInit();  
 
      Solver::GSiedel<double> gSiedel
                    ( pProblem->getMatric() 
             	    , pProblem->getApproxVector()
                    , pProblem->getRhsVector() ); 
       
      gSiedel.setMaxIteration(100000);
      gSiedel.setRelativeTolorance(0.00001);
      gSiedel.solve();

      double Error( gSiedel.ErrorInfNorm
                      ( pProblem->getExactSol() 
                      , pProblem->getApproxVector() ) 
                   ); 

       return  Error;
    }

    void GSiedelTest(){
 
        double OldError(1);
        bool TestPass(true);
     for(int i(2); i<6;i++){
    
       int _n1 = pow(2,i) - 1;

       double Error (generateProblemAndSolve(_n1) );
        
       double ratio = Error/OldError;
       OldError = Error;
        
       if( i != 2 ){ 
 
            if( (ratio > 0.25+0.2) && (ratio < 0.25-0.2)){
                 std::cout << "GSiedel Iterative Solver failed" << std::endl;
              }
        } 
       
      }
      
      if(TestPass){
       std::cout << "GSiedel Iterative Solver Passed" << std::endl;
      }
      
    }
   private: 

    int _n1;

    std::unique_ptr< PoissonBVP2D<double> > pProblem;
};

#endif
