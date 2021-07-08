add_executable( EX3_adv_diff src/main.f90 src/level.f90 src/sweeper.f90 src/probin.f90 src/hooks.f90 )
target_link_libraries( EX3_adv_diff PUBLIC pfasst )
target_compile_options( EX3_adv_diff  PRIVATE $<BUILD_INTERFACE:-cpp> )

