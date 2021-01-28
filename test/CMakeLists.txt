if( NOT TARGET BLAS::BLAS )
  find_package( BLAS REQUIRED )
endif()

add_subdirectory( magpicard )
add_subdirectory( imk )
#add_subdirectory( nagumo )
add_subdirectory( EXP_adv_diff_fft/1d )
