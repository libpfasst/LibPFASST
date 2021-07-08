add_executable( magpicard src/probin.f90 src/sweeper.f90 src/main.f90 src/hooks.f90
                          src/utils.f90 )
target_link_libraries( magpicard PUBLIC pfasst BLAS::BLAS )
target_compile_options( magpicard PRIVATE $<BUILD_INTERFACE:-cpp> )
