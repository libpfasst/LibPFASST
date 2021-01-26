add_executable( imk src/probin.f90 src/sweeper.f90 src/main.f90 src/hooks.f90
                    src/factory.f90 src/utils.f90 )
target_link_libraries( imk PUBLIC pfasst BLAS::BLAS )
target_compile_options( imk PRIVATE $<BUILD_INTERFACE:-cpp> )
