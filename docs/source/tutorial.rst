Tutorial
========

The following material will walk the user through a couple of examples to demonstrate how to set up an application using libpfasst.

Example 1
---------

Once libpfasst has been successfully built, move to the directory  libpfasst/Tutorials/EX1_Dahlquist.
This example solves the  scalar model problem or Dahlquist equation

.. math::

  y'  = \lambda y

  y(0) = 1

An implicit-explicit or IMEX  (also known as semi-implicit) splitting is used in this example, so the equation can be written 

.. math::
   y'  = \lambda_1 y + \lambda_2 y

Typing `make` in the directory should compile the example creating an exectuable called `main.exe`.  In the same directory, there are a few parameter files with the extension '.nml'.  You can run the example using one of these files, as in

    `$ ./main.exe sdc.nml`

Using your favorite editor, open the file sdc.nml.  There are two namelists here, 'PF_PARAMS' and 'PARAMS'.  The second of these is for local variables which in this simple example are just the values of lam1 and lam2, the simulation time and the number of steps.  
The 'PF_PARAMS' variables are discussed below in the section :doc:`Parameters <parameters>`.
Any parameter in either list can be overwritten by adding it to the command line, as in

    `$ ./main.exe sdc.nml lam1=3.0 niters=10`

The order of the command line parameters is not important except that they must come after the input file is specified.


The main program is in `src/main.f90`, and this format can easily be adapted to handle a different application.   
The routine to read the local parameters (located in  `src/probin.f90`)  is called first

.. code-block:: fortran
		
    !> read problem parameters
    call probin_init(probin_fname)

Most all of the changes needed to run a different example are contained in the block

.. code-block:: fortran
		
   !> loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)

       !>  Allocate the shape array for level (here just one dimension)
       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = 1

       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)
       call sweeper_setup(pf%levels(l)%ulevel%sweeper, pf%levels(l)%shape)

       !>  Set the size of the send/receive buffer
       pf%levels(l)%mpibuflen  = 1
    end do

The first line allocates the space for the `user_level`. In `src/feval.f90`, the local definition of the user level is given:

.. code-block:: fortran
		
  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: my_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type my_level_t

Here the routines to restrict and interpolate the solution in space are defined.  Since this is a scalar problem these routines are identity maps which can be seen by examining these routines at the bottom of   `src/feval.f90`.

In the same file, the definition of the sweeper being used is found

.. code-block:: fortran
		
  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imexQ_t) :: my_sweeper_t

   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves 

  end type my_sweeper_t

The IMEX sweeper needs 
functions to evaluate each term in the IMEX splitting and
a routine to solve an implicit equation equivalent to a
backward-Euler step.  These routines are in src/feval.f90 and are called
'f_eval' and 'f_comp'.
  

		
In EX1, the data structure of the solution is one of the included types in libpfasst called `ndarray`.
In the second line of the loop, the choice of the ndarray data structure is evident.

.. code-block:: fortran
		
       allocate(ndarray_factory::pf%levels(l)%ulevel%factory)


Here the solution is a one-dimensional array of length 1, and the  two lines

.. code-block:: fortran
		
       !>  Allocate the shape array for level (here just one dimension)
       allocate(pf%levels(l)%shape(1))
       pf%levels(l)%shape(1) = 1


allocate an integer array of the dimension of the ndarray, and then set the extent of each dimension.  Finally, the last line in the loop sets the total size of the problem on each level (again here it is a single variable).  This is used to control the size of the MPI buffers used to send and recieve the data.

The subroutine call that actually does the PFASST algorithm is 

.. code-block:: fortran
		
    !> do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)

Here we are passing the initial condition and solution at the final time, and these are type `ndarray` declared at the beginning of the subroutine

.. code-block:: fortran
		
    type(ndarray)     :: y_0      !<  the initial condition
    type(ndarray)     :: y_end    !<  the solution at the final time

The initial condition is set by a call to the subroutine `initial` which is in `src/feval.f90`

.. code-block:: fortran
		
    !> compute initial condition
    call initial(y_0)
   

The file `src/hooks.f90` contains some output routines, and in the main routine, a hook is set to call one of these routines after every iteration.  See the section on hooks for more info

.. code-block:: fortran
		
    !> add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_ITERATION, echo_error)

Example 2
---------

A second example is provided in the directory `Tutorials/EX2_Dahlquist` which solves exactly the same problem as example EX1.  The main difference in EX2 is that the ndarray encapsulation type is not used.  Instead, a user defined encapsulation type is created to demonstrate what a user needs to provide for a new application that does not use ndarray.  The data type here is called `scalar_encap` which is evident in the lines in `src/main.f90`

.. code-block:: fortran
		
    type(scalar_encap) :: y_0      !<  the initial condition
    type(scalar_encap) :: y_end    !<  the solution at the final time

If you compare the src/main.f90 files for EX1 and EX2, you will see they are very similar.    The only important change is in the allocation of the encapsulation factory
    
.. code-block:: fortran
		
       allocate(scalar_factory::pf%levels(l)%ulevel%factory)

       
Examining the difference between the files `src/probin.f90` and `src/hooks.f90` from the two examples will show that they too are nearly identical.  The big difference here is contained in the file `src/encap.f90`.  It is here that the user defined data type `scalar_encap` is defined.  Recall in EX1, we used the intrinsic data type `ndarray`, hence the reader might want to compare the file `src/encap.f90` with `libpfasst/src/pf_ndarray.f90`.

To create a new encapsulation of a data type for use by libpfasst, the user must provide the analog of the routines in `src/encap.f90`.
First, routines to create and destroy the data type (both a single copy and an array) must be provided.
The rhs of the following block defines the names of these user supplied routines

.. code-block:: fortran
		
  !>  Type to create and destroy the arrays
  type, extends(pf_factory_t) :: scalar_factory
   contains
     procedure :: create_single  => scalar_create_single
     procedure :: create_array  => scalar_create_array
     procedure :: destroy_single => scalar_destroy_single
     procedure :: destroy_array => scalar_destroy_array
  end type scalar_factory

Second, the data type must be defined, and here it is simply a scalar

.. code-block:: fortran
		
  !>  Type to extend the abstract encap and set procedure pointers		
  type, extends(pf_encap_t) :: scalar_encap
     real(pfdp) :: y   !  The scalar value

Finally, seven routines that act on the data type must be provided.

.. code-block:: fortran
		
     procedure :: setval => scalar_setval
     procedure :: copy => scalar_copy
     procedure :: norm => scalar_norm
     procedure :: pack => scalar_pack
     procedure :: unpack => scalar_unpack
     procedure :: axpy => scalar_axpy
     procedure :: eprint => scalar_eprint

The rhs of these statements defines the names of the user supplied subroutines.  The reader can see the actual subroutines in the remainder of `src/encap.f90` (they are all quite trivial).

Advection Diffusion Example
---------------------------

For a more complicated example,  see the Example in test/adv_diff_fft  included in libpfasst
for a simple PDE application of libpfasst.

This example solves a 1d linear advection diffusion equation

.. math::

  u_t  = - v u_x + \nu u_{xx}.

This PDE is set up to be solved in a way that is very similar to the EX1_Dahlquist example.  
This right hand side of the equation will be split into stiff terms handled implicitly
(:math:`\nu u_{xx}`) and non-stiff terms handled explicitly (:math:`-v u_x`),
and an IMEX SDC substepper will be used to evolve the equation through time.

The evaluation of the two rhs terms and the
solution of  implicit equation (feval and fcomp) are both
done in spectral space using the FFT.

The routines to spatially interpolate and
restrict solutions are in `src/feval.f90`.  It is assumed that the refinement ratio in space
is always 2.  Restriction is pointwise, and interpolation is done in spectral space again with the FFT.


Steps to build your own example
-------------------------------

The difficulting of using libpfasst to build a time parallelization scheme for a new application depends in large part on how different the application is to existing examples.  The key components of any example are

1. The data encapsulation type
2. The function evaluations
3. The sweeper type
4. The interpolation and restriction operators

In the simplest scenario, only the function evaluations would need to be changed.  This would be the case for example if the
included `ndarray` data type is used along with one of the sweepers included in `libpfasst` and no interpolation and restriction is
used.  

