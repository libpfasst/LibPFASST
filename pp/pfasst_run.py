#  This runs a single instane of PFASST and then loads all the output stats and timers
class PFASST_run:
      method="PFASST"
      def __init__(self, Nprocs,Nsteps,Name,nmlfile):
            import numpy as np
            self.Nsteps=Nsteps    #  Number of steps
            self.Nprocs=Nprocs    #  Number of processors
            self.basename = Name  #  Base output name for run
            self.nmlfile=nmlfile  #  Input file to read
            self.outdir="foo"     #  Will be padded with number of processors  
            self.Nlevs=1
            self.Niters=1
            self.Dim=1
            self.param_dict={}
            self.cmdopt=" "

      def do_run(self,do_run):
            import os            
            import numpy as np
            import json
            from pp.time_per_proc import Time_per_proc as timings
            from pp.stats_per_step import Stats_per_step as stats
            
#            exe='mpirun -n '+str(self.Nprocs)+ ' ./main.'+str(self.Dim)+'d.exe ' 
            exe='mpirun -n '+str(self.Nprocs)+ ' ./main.exe ' 
            

            #  Build output directory            
            self.outdir=self.basename+'P'+'{:04d}'.format(self.Nprocs)    
            outstr='\\\"'+self.basename+'\\\"'

            #  Build command line options
            vstr=self.nmlfile
            vstr=vstr+' outdir='+outstr
            vstr=vstr+' Nsteps='+str(self.Nsteps)

            #  add on external command line options
            myCmd =exe+vstr +self.cmdopt+' > out'
            if(do_run):
                  print(myCmd)
                  os.system(myCmd)  #  Run the code

            # Set up output directories
            datdir='dat/'+self.outdir+'/'   #  Base output directory
            procdir='/Proc_'+str(self.Nprocs-1).zfill(3)+'/'

            #  Load the dictionary of parameters
            with open(datdir+'pfasst_params.json', 'r') as f:
                  self.param_dict = json.load(f)

            #  Load the timer module
            if (self.param_dict['save_timings'] > 0):
                  self.run_times=timings(self.param_dict)

            #  Load the stats module
            self.iter_stats=stats(self.param_dict)
                  

