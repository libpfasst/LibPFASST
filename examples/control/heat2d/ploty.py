import numpy as np
import matplotlib.pyplot as plt
nnodes = 5;
nvars = 64;
nsteps = 4;
levelstring = 'l01m'

pltrng = range(2*nnodes,3*nnodes)

dat = np.zeros([nnodes*nsteps,nvars,nvars]);
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+ levelstring + str(m).zfill(2) + '.npy'
    fn='Dat/pfasst_V/numpy/pcompr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/ydr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    dat[counter,:,:] = np.load(fn)
    counter = counter + 1
#for i in pltrng:
  ##plt.figure(figsize=(16,12))
  #xx = np.linspace(0,1,nvars)
  #yy = np.linspace(0,1,nvars)
  #x, y = np.meshgrid(xx,yy)
  #plt.pcolormesh(x, y, dat[i,:,:])
  #plt.xlabel('x')
  #plt.ylabel('y')
  #plt.title('(y-yd) comp at t=t_' + str(i)+', max = '+str(np.max(dat[i,:,:])))
  #plt.colorbar()
##plt.legend()

plt.figure()
plt.plot(dat[:,16,16], marker='o')
plt.title('comp, max = '+str(np.max(np.abs(dat[:,16,16]))))

dat = np.zeros([nnodes*nsteps,nvars,nvars]);
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+ levelstring + str(m).zfill(2) + '.npy'
    fn='Dat/pfasst_V/numpy/pexr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/ydr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    dat[counter,:,:] = np.load(fn)
    counter = counter + 1
#for i in pltrng:
  #plt.figure()
  #xx = np.linspace(0,1,nvars)
  #yy = np.linspace(0,1,nvars)
  #x, y = np.meshgrid(xx,yy)
  #plt.pcolormesh(x, y, dat[i,:,:])
  #plt.xlabel('x')
  #plt.ylabel('y')
  #plt.title('error in y at t=t_' + str(i)+', max = '+str(np.max(dat[i,:,:])))
  #plt.colorbar()
plt.figure()
plt.plot(dat[:,16,16], marker='o')
plt.title('ex, max = '+str(np.max(np.abs(dat[:,16,16]))))
  
#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    ###fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+ levelstring + str(m).zfill(2) + '.npy'
    #fn='Dat/pfasst_V/numpy/prhsr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    ###fn='Dat/pfasst_V/numpy/ydr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    #counter = counter + 1
#for i in pltrng:
  #plt.figure()
  #xx = np.linspace(0,1,nvars)
  #yy = np.linspace(0,1,nvars)
  #x, y = np.meshgrid(xx,yy)
  #plt.pcolormesh(x, y, dat[i,:,:])
  #plt.xlabel('x')
  #plt.ylabel('y')
  #plt.title('F(i,1) for adjoint at t=t_' + str(i)+', max = '+str(np.max(dat[i,:,:])))
  #plt.colorbar()
  
#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    ###fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+ levelstring + str(m).zfill(2) + '.npy'
    #fn='Dat/pfasst_V/numpy/prhsdiffr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    ###fn='Dat/pfasst_V/numpy/ydr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    #counter = counter + 1
#for i in pltrng:
  #plt.figure()
  #xx = np.linspace(0,1,nvars)
  #yy = np.linspace(0,1,nvars)
  #x, y = np.meshgrid(xx,yy)
  #plt.pcolormesh(x, y, dat[i,:,:])
  #plt.xlabel('x')
  #plt.ylabel('y')
  #plt.title('error in F(i,1) at t=t_' + str(i)+', max = '+str(np.max(dat[i,:,:])))
  #plt.colorbar()


#datyd = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l01m' + str(m).zfill(2) + '.npy'
    #fn='Dat/pfasst_V/numpy/yexr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #datyd[counter,:,:] = np.load(fn)
    #datyd[counter,:,:] -= dat[counter,:,:]
    ##datyd[counter,:,:] *= (-1)
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
##plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()
#plt.figure(figsize=(16,12))
#xx = np.linspace(0,1,nvars)
#yy = np.linspace(0,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, datyd[nnodes-1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('yex-y at t=1')
#plt.colorbar()



#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+ levelstring + str(m).zfill(2) + '.npy'
    #fn='Dat/pfasst_V/numpy/pcompr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
    
#for i in pltrng:
  #plt.figure()
  #xx = np.linspace(0,1,nvars)
  #yy = np.linspace(0,1,nvars)
  #x, y = np.meshgrid(xx,yy)
  #plt.pcolormesh(x, y, dat[i,:,:])
  #plt.xlabel('x')
  #plt.ylabel('y')
  #plt.title('adjoint at t_' + str(i)+', max = '+str(np.max(dat[i,:,:])))
  #plt.colorbar()
##plt.savefig('adjoint_exp.png')

#datpex = np.zeros([nnodes*nsteps,nvars,nvars]);
#plt.figure(figsize=(16,12))
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/pexr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/y_dr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #datpex[counter,:,:] = np.load(fn)
    #datpex[counter,:,:] -= dat[counter,:,:]
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
##plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()

#xx = np.linspace(0,1,nvars)
#yy = np.linspace(0,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, datpex[1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('diff to exact adjoint pex-p at t=0')
#plt.colorbar()


#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/ur'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
##plt.figure(figsize=(16,12))
##xx = np.linspace(-1,1,nvars)
##yy = np.linspace(-1,1,nvars)
##x, y = np.meshgrid(xx,yy)
##plt.pcolormesh(x, y, dat[1,:,:])
##plt.xlabel('x')
##plt.ylabel('y')
##plt.title('u at T=0')
##plt.colorbar()


#datuex = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/uexactr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #datuex[counter,:,:] = np.load(fn);
    #datuex[counter,:,:] -= dat[counter,:,:];
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
#plt.figure(figsize=(16,12))
#xx = np.linspace(0,1,nvars)
#yy = np.linspace(0,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, datuex[nnodes-1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('diff to exact control, uex-u at t=1')
#plt.colorbar()

#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/udiffr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    #counter = counter + 1
#plt.figure(figsize=(16,12))
#xx = np.linspace(-1,1,nvars)
#yy = np.linspace(-1,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, dat[nnodes-1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('udiff at t=1')
#plt.colorbar()


plt.show()