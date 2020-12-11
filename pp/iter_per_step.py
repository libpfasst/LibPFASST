#  This holds all the iteration per step of a pfasst or parareal run
class Stats_per_step:
      Nsteps=1              #  Total number of steps
      Niters=1              #  Max iters from run
      Nblocks=1
      iters_per_step=[]
      resid_per_step=[]
      delq0_per_step=[]
      def __init__(self, param_dict,Nsteps):
            import numpy as np
            import json
            self.Nsteps = 1
            self.Nproc=1
            self.Nlev = param_dict['nlevels']
            self.Niters = param_dict['niters']
            self.Nproc = param_dict['nproc']

            self.iters_per_step=np.zeros([self.Nsteps])
            self.resid_per_step=np.zeros([self.Nsteps,self.Niters])
            self.delq0_per_step=np.zeros([self.Nsteps,self.Niters])
            Nblocks=int(self.Nsteps/self.Nproc)            
            for kProcs in  range(self.Nproc):
                iname='dat/'+param_dict['outdir']+'/residuals/Proc_'+str(kProcs).zfill(3)+'/Lev_'+str(self.Nlev).zfill(1)+'_iter.dat'
                rname='dat/'+param_dict['outdir']+'/residuals/Proc_'+str(kProcs).zfill(3)+'/Lev_'+str(self.Nlev).zfill(1)+'.dat'
                qname='dat/'+param_dict['outdir']+'/delta_q0/Proc_'+str(kProcs).zfill(3)+'/Lev_'+str(self.Nlev).zfill(1)+'.dat'
                k0=0

                iterarray=np.loadtxt(iname)
                resarray=np.loadtxt(rname)
                delqarray=np.loadtxt(qname)
                for ks in  range(Nblocks):
                    thisStep=int(ks*self.Nproc+kProcs+1)  # same as int(resarray[ks,0])
                    if (Nblocks == 1):
                        thisNiter=int(iterarray[2])
                        thisResid=resarray[k0:k0+thisNiter-1,4]
                        thisDelq=resarray[k0:k0+thisNiter-1,4]
                    else:
                        thisNiter=int(iterarray[ks,2])
                        thisResid=resarray[k0:k0+thisNiter-1,4]
                        thisDelq=resarray[k0:k0+thisNiter-1,4]
                    k0=k0+int(thisNiter)                        
                    self.iters_per_step[thisStep-1]=thisNiter
                    self.resid_per_step[thisStep-1,0:int(thisNiter)-1]=thisResid
                    self.delq0_per_step[thisStep-1,0:int(thisNiter)-1]=thisDelq

 
