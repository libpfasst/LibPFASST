#  This holds all the times per processor for pfasst or parareal run
class Time_per_proc:

      
      def __init__(self, param_dict):
            # On initialization, dat_dir should point to the directory where all the data lives (with trailing /)
            import numpy as np
            import json

            self.param_dict = param_dict
            if(self.param_dict['save_timings']==0):
                  print('Warning: No timing data for this run')
                  self.no_data=True
            else:
                  self.no_data=False
                  # Store the directory path
                  self.datdir='dat/'+param_dict['outdir']                        
                  
                  # Rip the problem params from the dictionary
                  self.Nprocs = self.param_dict['nproc']
                  self.Nlevs = self.param_dict['nlevels']
                  #  Set up a dictionary arrays to hold the data per timer
                  self.timers={}
                  self.timers['total']={'times': np.zeros(self.Nprocs),'do_plot': True,'color':'black'}
                  self.timers['predictor']={'times': np.zeros(self.Nprocs),'do_plot': True,'color':'orange'}
                  self.timers['block']={'times': np.zeros(self.Nprocs),'do_plot': False,'color':'green'}
                  self.timers['iteration']={'times': np.zeros(self.Nprocs),'do_plot': False,'color':'blue'}
                  self.timers['sweep']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'purple'}
                  self.timers['feval']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'red'}
                  self.timers['fcomp']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'lawngreen'}
                  self.timers['residual']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'yellow'}
                  self.timers['interp']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'aqua'}
                  self.timers['restrict']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'peru'}
                  self.timers['send']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'lightcoral'}
                  self.timers['receive']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'grey'}
                  self.timers['wait']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'darkblue'}
                  self.timers['hooks']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': True,'color':'pink'}
                  self.timers['aux']={'times': np.zeros([self.Nprocs,self.Nlevs]),'do_plot': False,'color':'slategrey'}
                  self.timer_per_proc=[]

                  #  Loop over the directory for each processor and load timers into dictionary
                  for kProcs in  range(self.Nprocs):
                        procdir=self.datdir+'/Proc_'+str(kProcs).zfill(3)+'/runtime.json'                  
                        with open(procdir, 'r') as f:
                              time_dict = json.load(f)    #  Loads the json file of timers
                              self.timer_per_proc.append(time_dict)
                              
                              for key in self.timers:
                                    if (self.timers[key]['times'].ndim > 1):
                                          self.timers[key]['times'][kProcs,:]=np.array(time_dict[key])
                                    else:
                                          self.timers[key]['times'][kProcs]=time_dict[key]
      #  Plot the different timers                        
      def plot_time_per_proc(self,fignum,savefig):                  
            import numpy as np
            import matplotlib.pyplot as plt
            if (self.no_data):
                  print('Warning: no timing data for this run')
            else:
                  
                  this_procs=np.arange(self.Nprocs)
                  plt.figure(fignum)
                  legstr1=[]
                  legstr2=[]
                  plt.subplot(1,2,1)                             
                  for key in self.timers:
                        if (self.timers[key]['do_plot']):
                              legstr1.append(key)
                              if (self.timers[key]['times'].ndim > 1):
                                    plt.plot(this_procs,np.sum(self.timers[key]['times'],axis=1),':*',color=self.timers[key]['color'])
                              else:
                                    plt.plot(this_procs,self.timers[key]['times'],':*',color=self.timers[key]['color'])
                  plt.legend(legstr1)    
                  plt.xlabel('Processor')
                  plt.ylabel('Time')
                                    
                  #  Plot the timers with multiple levels
                  legstr2=[]
                  plt.subplot(1,2,2)
                  linestyles = ['-' ,'--','-.', ':']
                  for pk in range(self.Nlevs):
                        for key in self.timers:
                              if (self.timers[key]['do_plot']):
                                    if (self.timers[key]['times'].ndim > 1):
                                          pkn=self.Nlevs-pk-1
                                          lst="$"+str(pkn+1)+"$"
                                          plt.semilogy(this_procs,self.timers[key]['times'][:,pkn],color=self.timers[key]['color'],linestyle=linestyles[pk],marker=lst,markersize=8)
#                                          plt.plot(this_procs,self.timers[key]['times'][:,pkn],color=self.timers[key]['color'],linestyle=linestyles[pk],marker=lst,markersize=8)
                                          if (pk == 0):
                                                legstr2.append(key)

                  plt.legend(legstr2)    
                  plt.xlabel('Processor')
                  plt.ylabel('Time')

                  if(savefig):
                        plt.savefig(self.datdir+'/time_per_proc.png',dpi=300)                        
