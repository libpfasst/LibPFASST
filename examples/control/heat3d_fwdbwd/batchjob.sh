#!/bin/bash
#SBATCH --job-name=heat3d
#SBATCH --ntasks=16
#SBATCH --partition=HTC040
###SBATCH --exclusive
#SBATCH --output=htc040_u1_n16.log
#SBATCH --open-mode=append
#SBATCH --mail-user=goetschel@zib.de
#SBATCH --mail-type=ALL


cd /nfs/datanumerik/bzfgoets/libpfasst/examples/control/heat3d_fwdbwd
# time mpirun -np 1 ./main_split.exe probin.nml do_both=0 >  opt16_u1_1.log
# time mpirun -np 4 ./main_split.exe probin.nml do_both=0 >  opt16_u1_4.log
# time mpirun -np 8 ./main_split.exe probin.nml do_both=0 >  opt16_u1_8.log
# time mpirun -np 16 ./main_split.exe probin.nml do_both=0 > opt16_u1_16.log
# time mpirun -np 1 ./main_split.exe probin.nml do_both=1 >  opt16_u1_both_1.log
# time mpirun -np 4 ./main_split.exe probin.nml do_both=1 >  opt16_u1_both_4.log
# time mpirun -np 8 ./main_split.exe probin.nml do_both=1 >  opt16_u1_both_8.log
time mpirun -np 16 ./main_split.exe probin.nml do_both=1 > opt16_u1_new_both_16.log
