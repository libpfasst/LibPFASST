import numpy as np
import matplotlib.pyplot as plt
nnodes = 9;
nvars = 128;
nsteps = 20;
dat = np.zeros([nnodes*nsteps,nvars]);

counter = 0
for step in range(0,nsteps):
 for m in range(1,nnodes+1):
   #fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l01m' + str(m).zfill(2) + '.npy'
   fn='Dat/pfasst_V/numpy/ur'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
   dat[counter,:] = np.load(fn);
   counter = counter + 1
plt.figure()
plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')#, vmin=-0.65, vmax=1.0)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('control')
###plt.savefig('u_misdc_r50_dy_wp_cold_1em11_256_7.png')
##plt.savefig('u_tozero_imex_r20_dy_wp_cold_1em11_128_5-5-9_200iter.png')


#plt.figure()
#plt.plot(dat[:,16])
##plt.savefig('dat16_misdc_r20_dy_wp_cold_1em7_128_11.png')#, dpi=plt.dpi)

#plt.figure()
#dat = np.zeros([nnodes*nsteps,nvars]);
#counter = 0
#for step in range(0,nsteps):
  #for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/uexactr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    #dat[counter,:] = np.load(fn);
    #counter = counter + 1
#plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')#, vmin=-0.65, vmax=0.0)
#plt.colorbar()
#plt.xlabel('x')
#plt.ylabel('t')
###v = np.linspace(-.65, 0., 15, endpoint=True)
###plt.colorbar(tics=v)
##plt.title('exact control')
###plt.savefig('exactcontrolmisdc10finefix.png')


plt.figure()
dat = np.zeros([nnodes*nsteps,nvars]);
counter = 0
for step in range(0,nsteps):
  for m in range(1,nnodes+1):
    fn='Dat/pfasst_V/numpy/udiffr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
    dat[counter,:] = np.load(fn);
    counter = counter + 1
plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('diff to exact control (unregularized problem)')
#plt.savefig('udiff_misdc_r20_dy_wp_cold_1em11_128_9_200iter.png')
##plt.savefig('diffcontrolmisdcnew.png')


plt.figure()
dat = np.zeros([nnodes*nsteps,nvars]);
counter = 0
for step in range(0,nsteps):
 for m in range(1,nnodes+1):
   fn='Dat/pfasst_V/numpy/y_s'+str(step).zfill(2)+'l03m' + str(m).zfill(2) + '.npy'
   #fn='Dat/pfasst_V/numpy/ur'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
   #print fn
   dat[counter,:] = np.load(fn);
   counter = counter + 1
plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('state')
####plt.savefig('y_misdc_r50_dy_wp_cold_1em11_256_7.png')
#plt.savefig('y_tozero_imex_r20_dy_wp_cold_1em11_128_5-5-9_200iter.png')

#plt.figure()
#plt.plot(dat[:,32])
###plt.plot(dat[nnodes*nsteps/2-1,:])
###plt.plot(dat[nnodes*nsteps-1,:])

plt.figure()
counter = 0
datyd = np.zeros([nnodes*nsteps,nvars]);
for step in range(0,nsteps):
 for m in range(1,nnodes+1):
   fn='Dat/pfasst_V/numpy/y_dr'+str(step).zfill(2)+'m'+ str(m).zfill(2) + '.npy'
   datyd[counter,:] = np.load(fn);
   counter = counter + 1
#dat = dat - datyd
plt.imshow(datyd,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
#plt.title('diff to desired state')
#plt.savefig('ydiff_tozero_imex_r20_dy_wp_cold_1em11_128_5-5-9_200iter.png')
#####plt.savefig('desiredmisdc10finefix.png')

plt.figure()
dat = np.zeros([nnodes*nsteps,nvars]);
counter = 0
#counter = counter - 1
for step in range(0,nsteps):
 for m in range(1,nnodes+1):
   fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+'l03m' + str(m).zfill(2) + '.npy'
   dat[counter,:] = np.load(fn);
   counter = counter + 1
plt.imshow(dat,extent=[0,20,0,5],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
#plt.gca().set_yticks(np.arange(0., 5., 0.2))
#plt.gca().invert_yaxis()
plt.title('adjoint')
#plt.savefig('p_tozero_imex_r20_dy_wp_cold_1em11_128_5-5-9_200iter.png')

#plt.savefig('p_misdc_r50_dy_wp_cold_1em11_256_7.png')

#plt.figure()
#dat = np.ones([2*nnodes+1,nvars]);
#counter = 0
#step = 10;
#for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+'l03m' + str(m).zfill(2) + '.npy'
    #dat[counter,:] = np.load(fn);
    #counter = counter + 1
##plt.plot(dat[:,10])  #plot for fixed x = idx*dx = 10 * 20/128 over time at single time step
##
##plt.figure()
##dat = np.zeros([nnodes,nvars]);
##counter = 0
#counter = counter+1
#step = 11;
#for m in range(1,nnodes+1):
    #fn='Dat/pfasst_V/numpy/p_s'+str(step).zfill(2)+'l03m' + str(m).zfill(2) + '.npy'
    #dat[counter,:] = np.load(fn);
    #counter = counter + 1
#masked = np.ma.masked_where(dat==1,dat)
#plt.plot(masked[:,10], linewidth=3.0)
#plt.savefig('p_full_ts10and11.svg',  format='svg', dpi=1200)


plt.show()
