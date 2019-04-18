#! /bin/bash
mpirun -n 1 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p1_n64_o1\"
mpirun -n 1 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p1_n64_o2\"
mpirun -n 1 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p1_n64_o3\"
mpirun -n 1 ./main.exe facke.nml nsteps=128 abs_res_tol=1e-8 magnus_order=1 fbase=\"test_p1_n128_o1\"
mpirun -n 1 ./main.exe facke.nml nsteps=128 abs_res_tol=1e-8 magnus_order=2 fbase=\"test_p1_n128_o2\"
mpirun -n 1 ./main.exe facke.nml nsteps=128 abs_res_tol=1e-8 magnus_order=3 fbase=\"test_p1_n128_o3\"
mpirun -n 1 ./main.exe facke.nml nsteps=256 abs_res_tol=1e-10 magnus_order=1 fbase=\"test_p1_n256_o1\"
mpirun -n 1 ./main.exe facke.nml nsteps=256 abs_res_tol=1e-10 magnus_order=2 fbase=\"test_p1_n256_o2\"
mpirun -n 1 ./main.exe facke.nml nsteps=256 abs_res_tol=1e-10 magnus_order=3 fbase=\"test_p1_n256_o3\"

mpirun -n 2 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p2_n64_o1\"
mpirun -n 4 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p4_n64_o1\"
mpirun -n 8 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p8_n64_o1\"
mpirun -n 16 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p16_n64_o1\"
mpirun -n 32 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p32_n64_o1\"
mpirun -n 64 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=1 fbase=\"test_p64_n64_o1\"
#

mpirun -n 2 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p2_n64_o2\"
mpirun -n 4 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p4_n64_o2\"
mpirun -n 8 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p8_n64_o2\"
mpirun -n 16 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p16_n64_o2\"
mpirun -n 32 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p32_n64_o2\"
mpirun -n 64 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=2 fbase=\"test_p64_n64_o2\"
#

mpirun -n 2 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p2_n64_o3\"
mpirun -n 4 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p4_n64_o3\"
mpirun -n 8 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p8_n64_o3\"
mpirun -n 16 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p16_n64_o3\"
mpirun -n 32 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p32_n64_o3\"
mpirun -n 64 ./main.exe facke.nml nsteps=64 abs_res_tol=1e-6 magnus_order=3 fbase=\"test_p64_n64_o3\"



