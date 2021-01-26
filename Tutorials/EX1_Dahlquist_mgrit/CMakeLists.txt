add_executable( EX1_Dahlquist_mgrit src/hooks.f90 src/stepper.f90 src/main.f90 src/level.f90 src/probin.f90 )
target_link_libraries( EX1_Dahlquist_mgrit PUBLIC pfasst )

