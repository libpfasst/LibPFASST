import numpy as np
import matplotlib.pyplot as plt
nnodes = 9;
nvars = 128;
nsteps = 20;
dat = np.zeros([nnodes*nsteps,nvars]);

##counter = 0
##for step in range(1,nsteps+1):
  ##for m in range(1,nnodes+1):
    ###fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l01m' + str(m).zfill(2) + '.npy'
    ##fn='Dat/pfasst_V/numpy/ur00'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##dat[counter,:] = np.load(fn);
    ##counter = counter + 1
##plt.figure()
##plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none', vmin=-0.65, vmax=0.0)
##plt.colorbar()
##plt.xlabel('x')
##plt.ylabel('t')
##plt.title('control')
####plt.savefig('controlmisdcnew.png')

##plt.figure()
##plt.plot(dat[10,:])
##plt.savefig('dat10.png', dpi=plt.dpi)

##plt.figure()
##dat = np.zeros([nnodes*nsteps,nvars]);
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/uexactr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##dat[counter,:] = np.load(fn);
    ##counter = counter + 1
##plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none', vmin=-0.65, vmax=0.0)
##plt.colorbar()
##plt.xlabel('x')
##plt.ylabel('t')
###v = np.linspace(-.65, 0., 15, endpoint=True)
###plt.colorbar(tics=v)
##plt.title('exact control')
##plt.savefig('exactcontrolmisdc10finefix.png')


##plt.figure()
##dat = np.zeros([nnodes*nsteps,nvars]);
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/udiffr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##dat[counter,:] = np.load(fn);
    ##counter = counter + 1
##plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()
##plt.xlabel('x')
##plt.ylabel('t')
##plt.title('diff in control')
###plt.savefig('diffcontrolmisdcnew.png')


##plt.figure()
##dat = np.zeros([nnodes*nsteps,nvars]);
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l03m' + str(m).zfill(2) + '.npy'
    ###fn='Dat/pfasst_V/numpy/ur'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ###print fn
    ##dat[counter,:] = np.load(fn);
    ##counter = counter + 1
##plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()
##plt.xlabel('x')
##plt.ylabel('t')
##plt.title('state')
##plt.savefig('statemisdcnew.png')

###plt.figure()
###plt.plot(dat[10,:])
###plt.plot(dat[nnodes*nsteps/2-1,:])
###plt.plot(dat[nnodes*nsteps-1,:])

##plt.figure()
##dat = np.zeros([nnodes*nsteps,nvars]);
##counter = 0
##for step in range(0,nsteps):
  ##for m in range(1,nnodes+1):
    ##fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+'l03m' + str(m).zfill(2) + '.npy'
    ##dat[counter,:] = np.load(fn);
    ##counter = counter + 1
##plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
##plt.colorbar()
##plt.xlabel('x')
##plt.ylabel('t')
##plt.gca().set_yticks(np.arange(0., 5., 0.2))
##plt.title('adjoint')
##plt.savefig('adjointmisdcnew.png')

#plt.figure()
#counter = 0
#datyd = np.zeros([nnodes*nsteps,nvars]);
#for step in range(1,nsteps+1):
  #index = 0
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/y_dr00s'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    ##print fn, counter + index
    #datyd[counter+index,:] = np.load(fn);
    #index = index + 1
  #counter = counter + index
#plt.imshow(datyd,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()
#plt.xlabel('x')
#plt.ylabel('t')
#plt.title('desired state')
####plt.savefig('desiredmisdc10finefix.png')

plt.figure()
dat = np.zeros([nnodes*nsteps,nvars]);
counter = 0
for step in range(0,nsteps):
  index = 0
  for m in range(1,nnodes+1):
    fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l02'+'m'+ str(m).zfill(2) + '.npy'
    #print fn, counter + index
    dat[counter+index,:] = np.load(fn);
    index = index + 1
  counter = counter + index
plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('state')

#plt.figure()
#dat = dat - datyd
#plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()
#plt.xlabel('x')
#plt.ylabel('t')
#plt.title('diff to desired state')


#plt.figure()
#dat = np.zeros([nnodes*nsteps,nvars]);
#counter = 0
#for step in range(1,nsteps+1):
  #index = 0
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/ur00s'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter+index,:] = np.load(fn);
    #index = index + 1
  #counter = counter + index
#plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()
#plt.xlabel('x')
#plt.ylabel('t')
#plt.title('control')
####plt.savefig('controlmisdcnew.png')

#plt.figure()
#dat = np.zeros([nnodes*nsteps,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #index = 0
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+'l03'+'m'+ str(m).zfill(2) + '.npy'
    ##print fn, counter + index
    #dat[counter+index,:] = np.load(fn);
    #index = index + 1
  #counter = counter + index
#plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
#plt.colorbar()
#plt.xlabel('x')
#plt.ylabel('t')
#plt.title('adjoint')

plt.show()
