#!/bin/bash
#SBATCH --job-name=heat3d
#SBATCH --ntasks=20
#SBATCH --partition=HTC050
#SBATCH --output=n10.log


cd /nfs/datanumerik/bzfgoets/libpfasst/examples/control/heat3d_flexible
time mpirun -np 1  ./main_split.exe > opt1.log
time mpirun -np 2  ./main_split.exe > opt2.log
time mpirun -np 5  ./main_split.exe > opt5.log
time mpirun -np 10 ./main_split.exe > opt10.log
time mpirun -np 20 ./main_split.exe > opt20.log

