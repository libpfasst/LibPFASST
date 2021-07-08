cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( LibPFASST VERSION 0.0.0 LANGUAGES Fortran C )

# DL robust FindBLAS/LAPACK
include( FetchContent )
FetchContent_Declare( linalg-cmake-modules 
  GIT_REPOSITORY https://github.com/wavefunction91/linalg-cmake-modules.git
  GIT_TAG        main
)
FetchContent_GetProperties( linalg-cmake-modules )
if( NOT linalg-cmake-modules_POPULATED )
  FetchContent_Populate( linalg-cmake-modules )
  list( APPEND CMAKE_MODULE_PATH ${linalg-cmake-modules_SOURCE_DIR} )
endif()


add_subdirectory( src )
add_subdirectory( Tutorials )
add_subdirectory( test )

