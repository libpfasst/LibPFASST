#  This holds all the stats per iteration per step of a pfasst or parareal run
#  In particular, total iteration counts, residuals, delta q0 (initial condition jumps), and errors
#  It is a work in progress
class Stats_per_step:
      def __init__(self, param_dict):
            import numpy as np
            import json
            self.param_dict=param_dict
            self.Nsteps = param_dict['nsteps']
            self.Nlevs = param_dict['nlevels']
            self.Niters = param_dict['niters']
            self.Nprocs = param_dict['nproc']
            self.max_Nsweeps = max(param_dict['nsweeps'])

            #  Make some space to hold stuff
            self.iter_per_step=np.zeros([self.Nsteps])
            self.resid_per_step=np.zeros([self.Nsteps])
            self.delq0_per_step=np.zeros([self.Nsteps])
            self.error_per_step=np.zeros([self.Nsteps])
            self.resid_per_iter=np.zeros([self.Nsteps,self.Niters+1])
            self.delq0_per_iter=np.zeros([self.Nsteps,self.Niters+1])
            self.error_per_iter=np.zeros([self.Nsteps,self.Niters+1])
            self.Nblocks=int(int(self.Nsteps)/int(self.Nprocs))

            #  Loop over processors output files and strip data into numpy arrays
            dpath='dat/'+param_dict['outdir']
            for kProcs in  range(self.Nprocs):
                  pdir='/Proc_'+str(kProcs).zfill(3)+'/'
                  iname=dpath+pdir+'iter.dat'
                  iarray=np.loadtxt(iname)
                  iterarray=iarray.reshape(self.Nblocks,3)

                  #  We are going to reshape the resid and q0 and error arrays                  
                  if (param_dict['save_residuals']):
                        rname=dpath+pdir+'residual.dat'
                        rarray=np.loadtxt(rname)
                        resarray=rarray[:,5].reshape(self.Nlevs,self.Nblocks,self.Niters+1,self.max_Nsweeps)
                  if (param_dict['save_delta_q0']):
                        qname=dpath+pdir+'delta_q0.dat'
                        qarray=np.loadtxt(qname)
                        q0array =qarray[:,5].reshape(self.Nlevs,self.Nblocks,self.Niters+1,self.max_Nsweeps)
                  if (param_dict['save_errors']>0):
                        ename=dpath+pdir+'error.dat'
                        earray=np.loadtxt(ename)
                        errarray =earray[:,5].reshape(self.Nlevs,self.Nblocks,self.Niters+1,self.max_Nsweeps)

                  for ks in  range(self.Nblocks):
                        thisStep=int(ks*self.Nprocs+kProcs+1)  # same as int(resarray[ks,0])
                        thisNiter=int(iterarray[ks,2])
                        self.iter_per_step[thisStep-1]=thisNiter

                        if (param_dict['method']=="PFASST"):                              
                              # Load the residual from the finest level  (first sweep) xxx
                              if (param_dict['save_residuals']):
                                    thisResid=resarray[self.Nlevs-1,ks,0:thisNiter,0]
                                    self.resid_per_iter[thisStep-1,0:thisNiter]=thisResid
                                    self.resid_per_step[thisStep-1]=self.resid_per_iter[thisStep-1,thisNiter-1]

                              # Load the delta_q0 from the first level  (first sweep) xxx
                              if (param_dict['save_delta_q0']):
                                    thisDelq=q0array[0,ks,0:thisNiter,0]
                                    self.delq0_per_iter[thisStep-1,0:thisNiter]=thisDelq
                                    self.delq0_per_step[thisStep-1]=self.delq0_per_iter[thisStep-1,thisNiter-1]
                                    
                              #  Load the error from the finest level
                              if (param_dict['save_errors']>0):
                                    thisError=errarray[self.Nlevs-1,ks,0:thisNiter,0]
                                    self.error_per_iter[thisStep-1,0:thisNiter]=thisError
                                    self.error_per_step[thisStep-1]=self.error_per_iter[thisStep-1,thisNiter-1]
                        else:   #  parareal run
                              # Load the residual from the finest level  (first sweep) xxx
                              if (param_dict['save_residuals']):
                                    thisResid=resarray[self.Nlevs-1,ks,0:thisNiter+1,0]
                                    self.resid_per_iter[thisStep-1,0:thisNiter+1]=thisResid
                                    self.resid_per_step[thisStep-1]=self.resid_per_iter[thisStep-1,thisNiter]
                                    
                              # Load the delta_q0 from the first level  (first sweep) xxx
                              if (param_dict['save_delta_q0']):                              
                                    thisDelq=q0array[0,ks,0:thisNiter+1,0]
                                    self.delq0_per_iter[thisStep-1,0:thisNiter+1]=thisDelq
                                    self.delq0_per_step[thisStep-1]=self.delq0_per_iter[thisStep-1,thisNiter]

                              # Load the error from the finest level  (first sweep) xxx
                              if (param_dict['save_errors']>0):                              
                                    thisError=errarray[self.Nlevs-1,ks,0:thisNiter+1,0]
                                    self.error_per_iter[thisStep-1,0:thisNiter+1]=thisError
                                    self.error_per_step[thisStep-1]=self.error_per_iter[thisStep-1,thisNiter]
                                    

      def plot_iter_per_pstep(self,fignum,savefig):                  
            import numpy as np
            import matplotlib.pyplot as plt
            plt.figure(fignum)
            
            for j in  range(self.Nblocks):
                  j0=j*self.Nprocs
                  this_iters=self.iter_per_step[j0:j0+self.Nprocs]
                  this_procs=np.arange(self.Nprocs)+j0+1
                  plt.plot(this_procs,this_iters,'-*')
            plt.xlabel('PStep')
            plt.ylabel('Total Iterations')
            plt.title('Total Iterations for Convergence Per Step')
            if (savefig):
                  plt.savefig('total_iter_per_step.png',dpi=300)
                  
      def plot_resid_per_iter_per_pstep(self,fignum,savefig):                  
            import numpy as np
            import matplotlib.pyplot as plt
            plt.figure(fignum)
            
            for j in  range(self.Nblocks):
                  j0=j*self.Nprocs
                  this_procs=np.arange(self.Nprocs)+j0+1
                  Niters=int(self.iter_per_step[j0+self.Nprocs-1])
                  for k in  range(Niters+1):
                        this_resids=self.resid_per_iter[j0:j0+self.Nprocs,k]
                        if (np.max(this_resids) > 0):
                              plt.semilogy(this_procs[this_resids>0],this_resids[this_resids>0],'-*')
            plt.xlabel('PStep')
            plt.ylabel('Residuals')
            plt.title('Residual Per Step Per Iteration')
            if (savefig):
                  plt.savefig('residual_iter_per_step.png',dpi=300)
                  
      def plot_delq0_per_iter_per_pstep(self,fignum,savefig):                  
            import numpy as np
            import matplotlib.pyplot as plt
            plt.figure(fignum)
            
            for j in  range(self.Nblocks):
                  j0=j*self.Nprocs
                  this_procs=np.arange(self.Nprocs)+j0+1
                  Niters=int(self.iter_per_step[j0+self.Nprocs-1])
                  for k in  range(Niters+1):
                        this_delta=self.delq0_per_iter[j0:j0+self.Nprocs,k]
                        if (np.max(this_delta) > 0):
                              plt.semilogy(this_procs[this_delta>0],this_delta[this_delta>0],'-*')
            plt.xlabel('PStep')
            plt.ylabel('Delta q0')
            plt.title('Initial Condition Jump Per Step Per Iteration')            
            if (savefig):
                  plt.savefig('delta_q0_iter_per_step.png',dpi=300)                        
                              
      def plot_error_per_iter_per_pstep(self,fignum,savefig):                  
            import numpy as np
            import matplotlib.pyplot as plt
            plt.figure(fignum)
            
            for j in  range(self.Nblocks):
                  j0=j*self.Nprocs
                  this_procs=np.arange(self.Nprocs)+j0+1
                  Niters=int(self.iter_per_step[j0+self.Nprocs-1])
                  for k in  range(Niters+1):
                        this_errors=self.error_per_iter[j0:j0+self.Nprocs,k]
                        if (np.max(this_errors) > 0):
                              plt.semilogy(this_procs[this_errors>0],this_errors[this_errors>0],'-*')
            plt.xlabel('PStep')
            plt.ylabel('Error')
            plt.title('Solution Error Step Per Iteration')            
            if (savefig):
                  plt.savefig('error_iter_per_step.png',dpi=300)                        
