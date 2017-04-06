LIBPFASST communicators
=======================

LIBPFASST includes several (time) communicators:

* MPI
* fake

The MPI communicator is the most actively developed and maintained.
The fake communicator is used for convergence testing (convergence
tests across many 'fake' cores can be performed on one core).


The number time steps taken |N| should be an integer multiple of the
number of PFASST processors |P| used.  If more time steps are taken
then PFASST processors used, the default behaviour of LIBPFASST is to
operate in |BLOCK| mode.

In |BLOCK| mode, |P| time steps are taken at a time.  Once these |P|
steps are completed, the final solution on the last time processor is
broadcast to all time processors (and is received as a new initial
condition) and the PFASST algorithm is started over (including the
predictor stage) again.  This process is repeated :math:`N/P` times
until all |N| time steps are complete.  This means that the first
processor may be idle for a time while it waits for the last processor
to finish (due to the burn in time associated with the predictor step)
after each block of |P| time steps are computed.

Alternatively, in |RING| mode, when a processor is finished iterating
on the time step that it is operating on, it moves itself to one time
step beyond the last processor.

At the beginning of each PFASST iteration (before receive requests are
posted) each processor: (i) determines if it has converged by checking
residuals, and (ii) receives a status message (this is blocking so
that the status send/recv stage is done in serial).  This status
message contains two pieces of information: whether the previous
processor is finished iterating (in which case the previous processor
will move itself to the end of the line), and how many of of the
previous processors are going to move to the end of the line.

If all residual conditions are met AND the previous processor is done
iterating, the current processor will decide that it is done iterating
too.  It will then send this information to the next processor (and
bump up the number of processors that are going to move by one).

If a processor is done iterating, then it will move to the end of the
line.  To do this, it does a blocking probe for a message with the
special "nmoved" tag and with source ``MPI_ANY_SOURCE``.  This special
message is sent by the first processor that is still iterating *and*
whose previous processor moved, **or** by the last processor in the
event that all processors finished at the same iteration.  In this
case, even though the last processor is done iterating, it will do one
final iteration before moving so that normal communication can
proceed.

In any case, all processors rotate their first/last markers by the
number of processors that moved.

Finally, if a processor steps beyond the number of requested time
steps it stop iterating.  In this case, other processors move their
last marker accordingly so that they do not include stopped processors
in their communications.


MPI
---

The PFASST MPI communicator can be used in conjunction with spatial
MPI communications.  The user is responsible for splitting
``MPI_COMM_WORLD`` appropriately (this can be done, eg, using a simple
colouring algorithm).

To support the |RING| mode, the MPI communicator maintains the
``first`` and ``last`` integers to keep track of the MPI ranks
corresponding to the first and last processors in the current block.

At the beginning of each PFASST iteration the MPI communicator posts
receive requests on all PFASST levels.


fake
----

The PFASST fake communicator allows us to study the convergence
behaviour of PFASST in |BLOCK| mode using only a single core.


.. |N| replace:: :math:`N`
.. |P| replace:: :math:`P`
.. |BLOCK| replace:: ``PF_WINDOW_BLOCK``
.. |RING| replace:: ``PF_WINDOW_RING``
