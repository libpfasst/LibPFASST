#!/bin/bash
#SBATCH --job-name=heat3d
#SBATCH --ntasks=48
#SBATCH --partition=HTC050
###SBATCH --exclusive
#SBATCH --output=htc052_u1_n48_new_large.log
#SBATCH --open-mode=append
#SBATCH --mail-user=goetschel@zib.de
#SBATCH --mail-type=ALL


cd /nfs/datanumerik/bzfgoets/libpfasst/examples/control/heat3d_flexible
# time mpirun -np 1 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_1.log
# time mpirun -np 2 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_2.log
# time mpirun -np 4 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_4.log
# time mpirun -np 8 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_8.log
# time mpirun -np 16 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_16.log
# time mpirun -np 32 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_32.log
# time mpirun -np 48 ./main_split.exe probin.nml do_mixed=0 >  opt48_large_u1_48.log
# time mpirun -np 1 ./main_split.exe probin.nml do_mixed=1 > opt48_large_stop_mixed_u1_1.log
# time mpirun -np 2 ./main_split.exe probin.nml do_mixed=1 > opt48_large_stop_mixed_u1_2.log
# time mpirun -np 4 ./main_split.exe probin.nml do_mixed=1 > opt48_large_stop_mixed_u1_4.log
# time mpirun -np 8 ./main_split.exe probin.nml do_mixed=1 > opt48_large_stop_mixed_u1_8.log
# time mpirun -np 16 ./main_split.exe probin.nml do_mixed=1 > opt48_large_stop_mixed_u1_16.log
# time mpirun -np 32 ./main_split.exe probin.nml do_mixed=1 > opt48_large_stop_mixed_u1_32.log
time mpirun -np 48 ./main_split.exe probin.nml do_mixed=1 > opt48_new_large_stop_mixed_u1_48.log
time mpirun -np 1 ./main_split.exe probin.nml do_mixed=0 >  opt48_new_large_u1_1.log
time mpirun -np 1 ./main_split.exe probin.nml do_mixed=1 > opt48_new_large_stop_mixed_u1_1.log


