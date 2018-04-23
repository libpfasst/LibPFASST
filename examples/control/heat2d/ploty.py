import numpy as np
import matplotlib.pyplot as plt
nnodes = 5;
nvars = 128;
nsteps = 1;
dat = np.zeros([nnodes*nsteps,nvars,nvars]);
plt.figure(figsize=(16,12))
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l02m' + str(m).zfill(2) + '.npy'
    #fn='Dat/pfasst_V/numpy/ydr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    dat[counter,:,:] = np.load(fn)
    #plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    counter = counter + 1
#plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()

xx = np.linspace(0,1,nvars)
yy = np.linspace(0,1,nvars)
x, y = np.meshgrid(xx,yy)
plt.pcolormesh(x, y, dat[nnodes-1,:,:])
plt.xlabel('x')
plt.ylabel('y')
plt.title('state at t=0')
plt.colorbar()
#plt.legend()

datyd = np.zeros([nnodes*nsteps,nvars,nvars]);
plt.figure(figsize=(16,12))
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l01m' + str(m).zfill(2) + '.npy'
    fn='Dat/pfasst_V/numpy/ydr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    datyd[counter,:,:] = np.load(fn)
    #plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    counter = counter + 1
#plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()

xx = np.linspace(0,1,nvars)
yy = np.linspace(0,1,nvars)
x, y = np.meshgrid(xx,yy)
plt.pcolormesh(x, y, datyd[nnodes-1,:,:])
plt.xlabel('x')
plt.ylabel('y')
plt.title('y_d at t=0')
plt.colorbar()


##dat = np.zeros([nnodes*nsteps,nvars,nvars]);
##plt.figure(figsize=(16,12))
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ###fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l01m' + str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/y_dr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##dat[counter,:,:] = np.load(fn)
    ###plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    ##counter = counter + 1
###plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
###plt.colorbar()

##xx = np.linspace(0,1,nvars)
##yy = np.linspace(0,1,nvars)
##x, y = np.meshgrid(xx,yy)
##plt.pcolormesh(x, y, dat[nnodes-1,:,:])
##plt.xlabel('x')
##plt.ylabel('y')
##plt.title('desired state at T=0.5')
##plt.colorbar()

dat = np.zeros([nnodes*nsteps,nvars,nvars]);
plt.figure(figsize=(16,12))
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+'l02m' + str(m).zfill(2) + '.npy'
    #fn='Dat/pfasst_V/numpy/y_dr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    dat[counter,:,:] = np.load(fn)
    #plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    counter = counter + 1
#plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()

xx = np.linspace(0,1,nvars)
yy = np.linspace(0,1,nvars)
x, y = np.meshgrid(xx,yy)
plt.pcolormesh(x, y, dat[1,:,:])
plt.xlabel('x')
plt.ylabel('y')
plt.title('adjoint at t=0')
plt.colorbar()

#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#plt.figure(figsize=(16,12))
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/pexr'+str(step).zfill(2)+'m' + str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/y_dr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
##plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()

#xx = np.linspace(0,1,nvars)
#yy = np.linspace(0,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, dat[1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('exact adjoint at t=0')
#plt.colorbar()


#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#plt.figure(figsize=(16,12))
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/uexactr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
##plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()

#xx = np.linspace(0,1,nvars)
#yy = np.linspace(0,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, dat[1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('uexact at T=0')
#plt.colorbar()

dat = np.zeros([nnodes*nsteps,nvars,nvars]);
plt.figure(figsize=(16,12))
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    fn='Dat/pfasst_V/numpy/ur'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    dat[counter,:,:] = np.load(fn)
    #plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    counter = counter + 1
#plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()

xx = np.linspace(-1,1,nvars)
yy = np.linspace(-1,1,nvars)
x, y = np.meshgrid(xx,yy)
plt.pcolormesh(x, y, dat[1,:,:])
plt.xlabel('x')
plt.ylabel('y')
plt.title('u at T=0')
plt.colorbar()


#dat = np.zeros([nnodes*nsteps,nvars,nvars]);
#plt.figure(figsize=(16,12))
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/u0r'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:,:] = np.load(fn)
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    #counter = counter + 1
##plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()

#xx = np.linspace(0,1,nvars)
#yy = np.linspace(0,1,nvars)
#x, y = np.meshgrid(xx,yy)
#plt.pcolormesh(x, y, dat[nnodes-1,:,:])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('u0 at T=0.5')
#plt.colorbar()

##dat = np.zeros([nnodes*nsteps,nvars,nvars]);
##counter = 0
##xx = np.linspace(0,1,nvars)
##yy = np.linspace(0,1,nvars)
##x, y = np.meshgrid(xx,yy)

##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ###fn='Dat/pfasst_V/numpy/ur'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l01m' + str(m).zfill(2) + '.npy'
    ##dat[counter,:,:] = np.load(fn)
    ##plt.figure()
    ##plt.pcolormesh(x, y, dat[counter,:,:])
    ##plt.xlabel('x')
    ##plt.ylabel('y')
    ##plt.title('y at node ' + str(m))
    ##plt.colorbar()

    ###plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    ##counter = counter + 1
###plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
###plt.colorbar()

##plt.figure(figsize=(16,12))
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/uexactr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##dat[counter,:] = np.load(fn);
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    ##counter = counter + 1
###plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
###plt.colorbar()
##plt.xlabel('x')
###plt.ylabel('t')
##plt.title('exact control')
##plt.legend()


##plt.figure(figsize=(16,12))
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/udiffr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##dat[counter,:] = np.load(fn);
    ##plt.plot(dat[counter,:], label='interval ' + str(step) + ' node ' + str(m))
    ##counter = counter + 1
###plt.imshow(dat,extent=[0,1,0,1],aspect='auto',origin='lower',interpolation='none')
###plt.colorbar()
##plt.xlabel('x')
###plt.ylabel('t')
##plt.title('error in computed control')
##plt.legend()


plt.show()