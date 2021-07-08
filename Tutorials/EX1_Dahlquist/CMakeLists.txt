add_executable( EX1_Dahlquist src/main.f90 src/level.f90 src/sweeper.f90 src/probin.f90 )
target_link_libraries( EX1_Dahlquist PUBLIC pfasst )
