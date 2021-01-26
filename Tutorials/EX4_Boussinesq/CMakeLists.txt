add_executable( EX4_Boussinesq src/main.f90 src/level.f90 src/sweeper.f90 src/probin.f90 src/hooks.f90 )
target_link_libraries( EX4_Boussinesq PUBLIC pfasst )

