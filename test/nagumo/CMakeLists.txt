add_executable( nagumo src/feval.f90 src/solutions.f90 src/probin.f90 
                       src/pf_optimization.f90 src/main_split.f90 src/level.f90 
                       src/numpy.c src/hooks.f90 )
target_link_libraries( nagumo PUBLIC pfasst )
target_compile_options( nagumo PRIVATE $<BUILD_INTERFACE:-cpp> )
