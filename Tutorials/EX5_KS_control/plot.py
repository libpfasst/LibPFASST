import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

nnodes = 5;
nvars = 1024;
nsteps = 100;

path = './npy/'

dat = np.zeros([nnodes*nsteps,nvars]);

x = np.linspace(0,32*np.pi,nvars)

#plt.figure()
#dat = np.zeros([nvars])
#fn = path+'uexact.npy'
#dat = np.load(fn)
#plt.plot(x,dat)

plt.figure()
dat = np.zeros([nsteps+1,nvars]);
counter = 0
for step in range(0,nsteps+1):
#  for m in range(1,nnodes+1):
   #fn=path+'ys'+str(step).zfill(4)+'m' + str(m).zfill(2) + '.npy'
   fn=path+'ys'+str(step).zfill(4) + '.npy'
   dat[counter,:] = np.load(fn);
   counter = counter + 1
plt.imshow(dat,extent=[0,101,0,60],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('state')

plt.figure()
dat = np.zeros([nsteps*nnodes,nvars]);
counter = 0
for step in range(1,nsteps+1):
 for m in range(1,nnodes+1):
   fn=path+'ytargets'+str(step).zfill(4)+'m' + str(m).zfill(2) + '.npy'
   #fn=path+'ys'+str(step).zfill(4) + '.npy'
   dat[counter,:] = np.load(fn);
   counter = counter + 1
plt.imshow(dat,extent=[0,101,0,60],aspect='auto',origin='lower',interpolation='none')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('target state')


#plt.figure()
#plt.plot(dat[:,15])
#plt.title('state at node 15')

#plt.figure()
#plt.plot(dat[:,16])
#plt.title('state at node 16')

#plt.figure()
#plt.plot(dat[:,-1])
#plt.title('state at node 1024')

plt.figure()
#plt.plot(dat[-1,:])

datuex = np.zeros([nvars])
fn = path+'uexact.npy'
datuex = np.load(fn)
plt.plot(datuex, label='exact')


#plt.figure()
datu = np.zeros([nvars])
fn = path+'uk0001.npy'
datu = np.load(fn)
plt.plot(datu, label='initial')


#plt.figure()
datu = np.zeros([nvars])
fn = path+'ufinal.npy'
datu = np.load(fn)
plt.plot(datu, label='final')
plt.legend(loc='upper right')
plt.title('controls')

plt.figure()
plt.plot(datu-datuex)
plt.title('diff computed - exact control')


plt.figure()
fn = path+'gradientk0001.npy'
datu = np.load(fn)
#plt.plot(datu,label='it 1')
fn = path+'gradientk0005.npy'
datu = np.load(fn)
plt.plot(datu,label='it 5')
fn = path+'gradientk0020.npy'
datu = np.load(fn)
plt.plot(datu,label='it 20')
fn = path+'gradientk0060.npy'
datu = np.load(fn)
plt.plot(datu,label='it 60')
plt.legend(loc='upper left')
plt.title('gradients')


plt.title('gradients')

plt.show()
