add_executable( imex src/main.f90 src/sweeper.f90 src/level.f90 src/hooks.f90 
                     src/probin.f90 src/phi_mod.f90 )
target_link_libraries( imex PUBLIC pfasst )
target_compile_options( imex PRIVATE $<BUILD_INTERFACE:-cpp> )
