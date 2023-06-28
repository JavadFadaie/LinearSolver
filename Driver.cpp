
#include <iostream>
#include "JacobiSolver.hpp"
#include "GuassSiedelSolver.hpp"
#include "AMGSolver.hpp"
#include "problem.hpp"
#include "csrstencil.hpp"
#include "AnisotropicPoissonBVP2D.hpp"
#include "ConvectionDiffBVP2D.hpp"
#include "MultiscalePossion.hpp"
#include "RotatedAnisotropicBVP2D.hpp" 
#include <chrono>


int main(int argc, char* argv[] ){

   
   int n1d =  atoi(argv[1]);

   int nThread = atoi(argv[2]);

   int iproblem = atoi(argv[3]);
   


   
   Problem<double> * pProblem;
   switch (iproblem) {

     case 0: {
      std::cout << "Solving isotropic Poisson BVP in 2D" << std::endl;
      pProblem =  ( new PoissonBVP2D<double> (n1d) );
      break;
      }

     case 1: {
      std::cout << "Solving anisotropic Poisson BVP in 2D" << std::endl;
      pProblem =  ( new AnisotropicPoissonBVP2D <double>(n1d,0.01) );
      break;
    }
 
     case 2: {
      std::cout << "Solving Convection Diffusion BVP in 2D" << std::endl;
      pProblem =  ( new ConvectionDiffBVP2D<double>(n1d, 1) );
      break;
    }

     case 3: {
      std::cout << "Solving Multiscale Possion in 2D" << std::endl;
      pProblem =  ( new MultiscalePoission<double>(n1d, 1000) );
      break;
    }

      case 4: {
      std::cout << "Solving Rotated Anisotropic Possion in 2D" << std::endl;
      pProblem =  ( new RotatedAnisotropicBVP2D<double>(n1d, 1000) );
      break;
    }


    } 

     
     pProblem->generateDiscretization();
     pProblem->randomApproxVectorInit(); 
  
     
     std::cout << std::endl;
     std::cout << "AMG method " << std::endl;
     std::cout << std::endl;
     auto amgstart (std::chrono::high_resolution_clock::now());
     Solver::AMGSolver<double> amg( pProblem->getMatric()
                                  , pProblem->getApproxVector()
                                  , pProblem->getRhsVector()
                                  );
    
     amg.strengthThreshlod(0.25);
     amg.numPreSmoothing(1);
     amg.numPostSmoothing(1);
     amg.gridThresholdSize(0);
     amg.maxLevelAMG(10);
     amg.setMaxIteration(20);
     amg.setRelativeTolorance(1.e-7);
     amg.setNumThreads(nThread );
     amg.solve();
     auto amgstop (std::chrono::high_resolution_clock::now());
     
     std::cout <<"AMG method execution time " <<std::chrono::duration_cast<std::chrono::microseconds>
		                         (amgstop-amgstart).count() / 1.e6 << std::endl;
		                         
		                         
		                         
		                         
		                         
		                         
		                          
     std::cout << std::endl;
     std::cout << "Jacobi iterative method " << std::endl;
     std::cout << std::endl;
     auto jacobistart (std::chrono::high_resolution_clock::now());
     
      Solver:: JacobiSolver<double> jacobi( pProblem->getMatric()
                                      , pProblem->getApproxVector()
                                      , pProblem->getRhsVector()
                                     );
      jacobi.setRelativeTolorance(1.e-7);                               
      jacobi.solve();
      auto jacobistop (std::chrono::high_resolution_clock::now());
     
     std::cout <<"Jacobi method execution time " <<std::chrono::duration_cast<std::chrono::microseconds>
		                         (jacobistop-jacobistart).count() / 1.e6 << std::endl;
	 
  return 0;
  
 }




