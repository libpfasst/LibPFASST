#  This is the class to do a parareal run and get all the statistics
#  If Nlev=1, than this is just a serial stepper, so there will be no residuals or delta_q0 jumps, and the number of iters should be zero Niters=0
class Parareal_run:
      def __init__(self, Nprocs,Nsteps,Name,nmlfile):
            import numpy as np
            self.Nprocs=Nprocs
            self.Nsteps=Nsteps
            self.basename = Name
            self.nmlfile=nmlfile

            #  stuff that will be set from outside
            self.Dim=1
            self.Nlevs=2
            self.rk_order=4
            self.Niters=1
            self.Nstep_rk_coarse=1
            self.Nstep_rk_fine=32
            self.alpha=1
            #  Place to tack on command line options
            self.vstr=" "
            self.cmdopt=" "
            #  internal stuff
            self.param_dict={}
            
      def do_run(self,execute):
            import os            
            import numpy as np
            import json
            from time_per_proc import Time_per_proc as timings
            from stats_per_step import Stats_per_step as stats
            from errors_per_step import Errors_per_step as errors

            exe='mpirun -n '+str(self.Nprocs)+ ' ./main.'+str(self.Dim)+'d.exe '

            #  Build output directory
            self.outdir=self.basename+'P'+'{:04d}'.format(self.Nprocs)                
            outstr='\\\"'+self.basename+'\\\"'

            #  Build command line options
            self.vstr=self.nmlfile
            self.vstr+=' outdir='+outstr
            self.vstr+=' Nsteps='+str(self.Nsteps)


            #  add on external command line options
            myCmd =exe+self.vstr +self.cmdopt+' > out'
            if (execute):
                  print(myCmd)
                  os.system(myCmd)  #  Run the code

            #  Now get some stats from the run
            datdir='dat/'+self.outdir+'/'   #  Base output directory

            #  Get the run parameters in dictionary form
            with open(datdir+'pfasst_params.json', 'r') as f:
                  self.param_dict = json.load(f)

            #  Load the timer module
            if (self.param_dict['save_timings'] > 0):
                  self.run_times=timings(self.param_dict)

            #  Load the run stats, there are no residuals nor delta_q0 for one level
            if (self.Nlevs > 1):
                  self.iter_stats=stats(self.param_dict)

                        
