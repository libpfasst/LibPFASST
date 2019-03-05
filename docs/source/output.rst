
Output
========================

There are two main ways to produce output automatically from LibPFASST, namely through the setting of hooks and the output of the results module.

hooks
-----
The first output possibility is through the use of the what are called hooks, which are user defined subroutines that can be called from many different places in the source code.  These are typically set in the users main routine.  For example in "LibPFASST/Tutorials/EX1_Dahlquist/src/main.f90"  a hook is added by the statement.

.. code-block:: fortran
		
    !> add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

This will force a call to  the user-supplied subroutine "echo_error" (supplied locally in the file "hooks.f90") after every PFASST iteration.      The second parameter specifies the level on which the hooks should be called with the value -1 indicating all levels.  See the file
"LibPFASST/src/pf_hooks.f90" for a list of hook locations.

results
-------
LibPFASST contains the module "pf_mod_results" which is used to collect information about the run during a call to library.
The destination of the output can be changed from default direct "dat/" to an arbitrary path by specifying the variable "outdir" in the PF namelist.  As with all input flags, this can be done on the command line with the syntax outdir=\\"my_outdir\\".  The collection and output of the results can be controlled by the parameters save_residuals, save_timings, and echo_timings.  The first will determine if the residuals are saved internally and
then output to the directory "my_outdir/residuals" after the completion of the iterations.  The second does essentially the same for total timing of various parts of the code.  The last, echo_timings prints timings to  stdout while the code is running.

Note that for the library supplied data encapsulations there is a subroutine called "eprint" which dumps (some of) the content in the data type to the screen (typically for debugging purposes).  This routine could be modified for an additional avenue for output from any place in the code.




