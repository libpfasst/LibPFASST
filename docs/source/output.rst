
Output
========================

There are two main ways to produce output automatically from LibPFASST, namely through the setting of hooks and the output of the results module.

hooks
-----
The first output possibility is through the use of what are called hooks, which are user defined subroutines that can be called from many different places in the source code.  These are typically set in the users main routine.  For example in "LibPFASST/Tutorials/EX2_Dahlquist/src/main.f90"  a hook is added by the statement.

.. code-block:: fortran
		
    !> add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

This will force a call to  the user-supplied subroutine "echo_error" (supplied locally in the file "hooks.f90") after every PFASST iteration.      The second parameter specifies the level on which the hooks should be called with the value -1 indicating all levels.  See the file
"LibPFASST/src/pf_hooks.f90" for a list of hook locations.

results
-------
LibPFASST contains the module "pf_mod_results" which is used to collect information about the run during a call to LibPFASST.
The destination of the output can be changed from the default directory "dat/outdir" to an arbitrary path by specifying the variable "outdir" in the PF namelist.  As with all input flags, this can be done on the command line with the syntax outdir=\\\\"my_outdir\\\\".  The collection and output of the results can be controlled by the parameters save_residuals, save_delta_q0, save_errors, and save_timings.  The first will determine if the residuals are saved internally and
then output to the directory "my_outdir/residuals" after the completion of the iterations.  The second does essentially the change in intial conditions at each processor.  The third controls the behavior for errors. The last is for the output of
timings of various parts of the code and is an integer which can be set to 0,1,2, or 3.

* 0: no timings
* 1: just the total time
* 2: time all parts of the code
* 3: time all parts of the code and echo timers to the standard out

Note that for the library supplied data encapsulations there is a subroutine called "eprint" which dumps (some of) the content in the data type to the screen (typically for debugging purposes).  This routine can be modified for an additional avenue for output from any place in the code.

Depending on the above flag, several files will be created after the time integration and put in the directory /dat/outdir.  

* pfasst_params.json :  A description of all the run parameters
* runtimes :  A directory containing timing files for each processor
* residuals :  A directory containing residual information for each   processor
* delta_q0 :  A directory containing information for each  processor  concerning the changes to the initial condition at each processor
* errors :  A directory containing solution error files for each processor

The first 5 of these files are generated automatically  by the library unless the user turns of the option for each with flags in the input files.  The error statistics must be set by the user (typically with a "hook" since the library has no information on the error for a particular application.

At the moment, the easiest way to understand the output file formats is probably to consult the code in src/pf_results.f90, e.g.


.. code-block:: fortran
		
  
    do klevel=1,this%nlevs
       do kblock = 1, this%nblocks
          nstep=(kblock-1)*this%nprocs+this%rank+1
          do kiter = 0 , this%niters
             do ksweep = 1, this%max_nsweeps
                if (this%save_residuals) write(rstream,101 ) klevel,nstep,kblock,kiter,ksweep,this%residuals(klevel, kblock,kiter+1, ksweep)
                if (this%save_errors) write(estream,101) klevel,nstep,kblock,kiter,ksweep,this%errors(klevel,kblock,kiter+1,  ksweep)
                if (this%save_delta_q0) write(qstream,101) klevel,nstep,kblock,kiter,ksweep,this%delta_q0(klevel, kblock,kiter+1, ksweep)
             end do
          end do
       enddo
    enddo
    101 format(I3,I10, I10,I5, I4, e22.14)


The output files are set up to allow easy input into Python post-processing routines using  json.load or numpy.loadtxt routines (or equivalents).  


