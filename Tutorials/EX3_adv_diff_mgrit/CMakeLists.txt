add_executable( EX3_adv_diff_mgrit src/main.f90 src/level.f90 src/stepper.f90 src/probin.f90 src/hooks.f90 src/sweeper.f90 )
target_link_libraries( EX3_adv_diff_mgrit PUBLIC pfasst )
target_compile_options( EX3_adv_diff_mgrit  PRIVATE $<BUILD_INTERFACE:-cpp> )

