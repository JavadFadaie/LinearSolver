
#include <iostream>
#include "../iterativeSolver/JacobiSolver.hpp"
#include "../iterativeSolver/GuassSiedelSolver.hpp"
#include "../iterativeSolver/AMGSolver.hpp"
#include "iterativeSolverTest.hpp"



int main( ){
 
   IterativeSolverTest test;
 
   test.GSiedelTest();
   
 
  return 0;
  
 }




