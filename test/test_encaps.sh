#! /bin/bash


echo "Testing zNdarray encap"	
echo "Making ERK"	
cd zNdarray/ERK_stepper
make
make DIM=2
make DIM=3
echo "Running ERK 1-d"	
mpirun -n 4 ./main.1d.exe test_1d.nml 
echo "Running ERK 2-d"
mpirun -n 4 ./main.2d.exe test_2d.nml
echo "Running ERK 3-d"
mpirun -n 4 ./main.3d.exe test_3d.nml
echo "Making RK "
cd ../RK_stepper
make
make DIM=2
make DIM=3  
echo "Running RK 1-d"	
mpirun -n 4 ./main.1d.exe test_1d.nml &&
echo "Running RK 2-d"
mpirun -n 4 ./main.2d.exe test_2d.nml 
echo "Running RK 3-d"
mpirun -n 4 ./main.3d.exe test_3d.nml 
echo "Making Running IMEX SDC "	
cd ../IMEX_sweeper
make
make DIM=2
make DIM=3 
echo "Running IMEX SDC 1-d"
mpirun -n 4 ./main.1d.exe test_1d.nml
echo "Running IMEX SDC 2-d"	
mpirun -n 4 ./main.2d.exe test_2d.nml 
echo "Running IMEX SDC 3-d"
mpirun -n 4 ./main.3d.exe test_3d.nml 
echo "Making Running EXP SDC "
cd ../EXP_sweeper &&  make && make DIM=2 && make DIM=3  
echo "Running EXP SDC 1-d"	
mpirun -n 4 ./main.1d.exe test_1d.nml 
echo "Running EXP SDC 2-d"	
mpirun -n 4 ./main.2d.exe test_2d.nml 
echo "Running EXP SDC 3-d"	
mpirun -n 4 ./main.3d.exe test_3d.nml 
cd ../..
echo "Testing Ndarray encap"
cd Ndarray/IMEX_sweeper
make
make DIM=2
make DIM=3 
echo "Running IMEX SDC 1-d"
mpirun -n 8 ./main.1d.exe test_1d.nml
echo "Running IMEX SDC 2-d"	
mpirun -n 4 ./main.2d.exe test_2d.nml 
echo "Running IMEX SDC 3-d"
mpirun -n 4 ./main.3d.exe test_3d.nml 

