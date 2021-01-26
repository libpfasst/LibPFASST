add_executable( EX5_KS_control src/main.f90 src/level.f90 src/sweeper.f90 src/probin.f90 src/hooks.f90 src/pf_optimization_1d.f90 )
target_link_libraries( EX5_KS_control PUBLIC pfasst )

