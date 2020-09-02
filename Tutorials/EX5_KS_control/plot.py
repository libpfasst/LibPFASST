import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

nnodes = 5;
nvars = 1024;
nsteps = 500;

path = './npy/'

dat = np.zeros([nnodes*nsteps,nvars]);

x = np.linspace(0,32*np.pi,nvars)

#plt.figure()
#dat = np.zeros([nvars])
#fn = path+'uexact.npy'
#dat = np.load(fn)
#plt.plot(x,dat)

#plt.figure()
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
plt.plot(dat[-1,:])
dat = np.zeros([nvars])
fn = path+'uexact.npy'
dat = np.load(fn)
plt.plot(dat)

#plt.figure()
#plt.plot(dat[:,15*16])
#plt.title('state at node 15')

#plt.figure()
#plt.plot(dat[:,16*16])
#plt.title('state at node 16')

#plt.figure()
#plt.plot(dat[:,32*16])
#plt.title('state at node 32')


plt.show()
