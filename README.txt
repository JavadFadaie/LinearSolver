README

With DriverAMG.cpp we try to solve several linear problem (Ax=b) with AMG (Algebraic Multi-Grid) linear Solver.

To compile redirect to build directory: cd/build
                                        cmake ..
                                        make
                                           
to run the program redirect to build directory type:./AMG 100  4  0 

100 : The problem size  
 4  : Number of thread used 
 0  : Type of the problem solved (2D Possion) 

