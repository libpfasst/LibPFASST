#! /bin/bash
mpirun -n 1 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-4 magnus_order=1 outdir=\" ./test_p1_n64_o1 \"
mpirun -n 1 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-5 magnus_order=2 outdir=\"./test_p1_n64_o2\"
mpirun -n 1 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 outdir=\"./test_p1_n64_o3\"
mpirun -n 1 ./main.exe facke.nml nsteps=128 abs_res_tol=1e-4 magnus_order=1 outdir=\"./test_p1_n128_o1\"
mpirun -n 1 ./main.exe facke.nml nsteps=128 abs_res_tol=1e-6 magnus_order=2 outdir=\"./test_p1_n128_o2\"
mpirun -n 1 ./main.exe facke.nml nsteps=128 abs_res_tol=1e-8 magnus_order=3 outdir=\"./test_p1_n128_o3\"
mpirun -n 1 ./main.exe facke.nml nsteps=256 abs_res_tol=1e-5 magnus_order=1 outdir=\"./test_p1_n256_o1\"
mpirun -n 1 ./main.exe facke.nml nsteps=256 abs_res_tol=1e-7 magnus_order=2 outdir=\"./test_p1_n256_o2\"
mpirun -n 1 ./main.exe facke.nml nsteps=256 abs_res_tol=1e-10 magnus_order=3 outdir=\"./test_p1_n256_o3\"




