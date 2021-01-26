
set( LIBPFASST_SRC 
     pfasst.f90 
     pf_comm.f90 
     pf_dtype.f90 
     pf_fft.f90 
     pf_hooks.f90 
     pf_interpolate.f90 
     pf_parallel.f90 
     pf_parareal.f90 
     pf_parallel_oc.f90 
     pf_pfasst.f90 
     pf_quadrature.f90 
     pf_restrict.f90 
     pf_results.f90 
     pf_rkstepper.f90 
     pf_erkstepper.f90 
     pf_solutions.f90 
     pf_stop.f90 
     pf_timer.f90 
     pf_utils.f90 
     pf_ham_sys_encap.f90 
     pf_ndarray_encap.f90 
     pf_ndarray_oc_encap.f90 
     pf_ndsysarray_encap.f90 
     pf_zndsysarray_encap.f90 
     pf_zndarray_encap.f90 
     pf_exp_sweeper.f90 
     pf_fexp_sweeper.f90 
     pf_imex_sweeper.f90 
     pf_imexR_sweeper.f90 
     pf_imexQ_oc_sweeper.f90 
     pf_imk_sweeper.f90 
     pf_magpicard_sweeper.f90 
     pf_misdcQ_sweeper.f90 
     pf_misdcQ_oc_sweeper.f90 
     pf_verlet_sweeper.f90 
     pf_fexp_sweeper.f90 
     pf_mgrit.f90 
#    pf_amisdc.f90 
#    pf_erkstepper.f90 
#    pf_amisdcQ.f90 
)

list( APPEND LIBPFASST_SRC dfftpack.f pf_fftpack.f90 ) # If no FFTW

find_package( MPI REQUIRED )
list( APPEND LIBPFASST_SRC pf_mpi.f90 )   # If MPI



add_library( pfasst ${LIBPFASST_SRC} )
target_link_libraries( pfasst PUBLIC MPI::MPI_C )
target_include_directories( pfasst PUBLIC ${CMAKE_BINARY_DIR}/src )
target_compile_options( pfasst PRIVATE $<BUILD_INTERFACE:-cpp> )
