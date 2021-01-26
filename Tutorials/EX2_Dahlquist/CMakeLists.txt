add_executable( EX2_Dahlquist src/main.f90 src/level.f90 src/sweeper.f90 src/probin.f90 src/hooks.f90 src/encap.f90 )
target_link_libraries( EX2_Dahlquist PUBLIC pfasst )
target_compile_options( EX2_Dahlquist  PRIVATE $<BUILD_INTERFACE:-cpp> )

