# Script to do a plot of residual, delta_q0 and total iterations  per step single run
import numpy as np
import matplotlib.pyplot as plt
from pp.pfasst_run import PFASST_run as pf_run
import os


#  Define some things needed for run
Nprocs=4                    #  Number of processors
Nsteps=16
basename='test'  #  Base name output directory for runs
nmlfile='multi_level.nml'            #  Input file
#  Set up run
this_run=pf_run(Nprocs,Nsteps,basename,nmlfile)

#  Set some things that are needed for run
#this_run.Nlevs=2
#this_run.Dim=1
#this_run.Niters=20

#  Do the run
this_run.do_run(do_run=True)

#  Plot the total iterations by step
this_run.iter_stats.plot_iter_per_pstep(fignum=1,savefig=False)

#  Plot the residual for each iterations by step
this_run.iter_stats.plot_resid_per_iter_per_pstep(fignum=2,savefig=False)

#  Plot the initial condition jump for each iterations by step
this_run.iter_stats.plot_delq0_per_iter_per_pstep(fignum=3,savefig=False)

#  Plot the errors for each iterations by step
this_run.iter_stats.plot_error_per_iter_per_pstep(fignum=4,savefig=False)

#  Plot the timer info for each processor
this_run.run_times.plot_time_per_proc(fignum=5,savefig=True)



plt.show()



