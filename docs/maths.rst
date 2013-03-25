Maths
=====

The PFASST algorithm is, in some ways, similar to the parareal method
for the temporal parallelization of ODEs and PDEs.  These methods can
be roughly described in terms of two numerical approximation methods,
denoted here by :math:`\mathcal{G}` and :math:`\mathcal{F}`.  Both
:math:`\mathcal{G}` and :math:`\mathcal{F}` propagate an initial value
:math:`U_n \approx u(t_n)` by approximating the solution to

.. math::

  u(t) = u_0 + \int_0^t f(\tau,u(\tau)) \,d\tau

from :math:`t_n` to :math:`t_{n+1}` (that is, :math:`\dot{u} =
f(u,t)`).  For the methods to be efficient, it must be the case that
the :math:`\mathcal{G}` propagator is computationally less expensive
than the :math:`\mathcal{F}` propagator, and hence
:math:`\mathcal{G}` is typically a low-order method.  Since the
overall accuracy of the methods is limited by the accuracy of the
:math:`\mathcal{F}` propagator, :math:`\mathcal{F}` is typically
higher-order and in addition may use a smaller time step than
:math:`\mathcal{G}`.  For these reasons, :math:`\mathcal{G}` is
referred to as the coarse propagator and :math:`\mathcal{F}` the fine
propagator.

For PDEs, it may be the case that the coarse and fine propagators
differ only in their level of refinement (in the spatial and/or
temporal directions).

The parareal method begins by computing a first approximation
:math:`U_{n+1}^0` for :math:`n = 0 \ldots N-1` where :math:`N` is the
number of time steps, often performed with the coarse propagator:

.. math::

   U_{n+1}^0 = \mathcal{G}(t_{n+1}, t_{n}, U_n^0).

The parareal method then proceeds iteratively, alternating between the
parallel computation of :math:`\mathcal{F}(t_{n+1},t_n,U_n^k)` and an
update of the initial conditions at each processor of the form

.. math::

  U_{n+1}^{k+1} = \mathcal{G}(t_{n+1}, t_n, U_n^{k+1})
                   + \mathcal{F}(t_{n+1}, t_n, U_n^k)
                   - \mathcal{G}(t_{n+1}, t_n, U_n^{k})

for :math:`n = 0 \ldots N-1`.  That is, the fine propagator is used
to refine the solution in each time slice in parallel, while the
coarse propagator is used to propagate the refinements performed by
the fine propagator through time to later processors.

The PFASST method operaters in a similar manner to the parareal
method, but it employs the SDC time-integration technique in a novel
way to improve the parallel efficiency of the method.  However, the
general structure and roles of the fine and coarse propagators are the
same as the parareal method.

Furthermore, the PFASST algorithm adds corrections to the coarse
propagator in a manner reminiscent of FAS corrections for non-linear
multi-grid methods.


References
----------

