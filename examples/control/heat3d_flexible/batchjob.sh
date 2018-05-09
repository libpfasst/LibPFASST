#!/bin/bash
#SBATCH --job-name=heat3d
#SBATCH --ntasks=32
#SBATCH --partition=HTC050
###SBATCH --exclusive
#SBATCH --output=htc052_n32.log


cd /nfs/datanumerik/bzfgoets/libpfasst/examples/control/heat3d_flexible
time mpirun -np 1 ./main_split.exe > opt32_1.log
time mpirun -np 2 ./main_split.exe > opt32_2.log
time mpirun -np 4 ./main_split.exe > opt32_4.log
time mpirun -np 8 ./main_split.exe > opt32_8.log
time mpirun -np 16 ./main_split.exe > opt32_16.log
time mpirun -np 32 ./main_split.exe > opt32_32.log

