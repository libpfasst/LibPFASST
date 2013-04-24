import h5py
import sys
import getMesh
import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import pylab

filename='test.h5'
file = h5py.File(filename,'r')
Nx_dset = file.get('input/discretization/Nx')

# corresponding dataset for some reasons is an array in an array; [0][0] retrieves integer
Nx           = Nx_dset[()][0][0]
xleft_dset   = file.get('input/problemdefinition/x_left')
x_left       = xleft_dset[()][0][0]
xright_dset  = file.get('input/problemdefinition/x_right')
x_right      = xright_dset[()][0][0]
xnodes  = numpy.linspace(x_left, x_right, Nx+1)
xcenter = [0 for x in range(Nx)]
for ii in range(0,Nx):
    xcenter[ii] = 0.5*( xnodes[ii] + xnodes[ii+1])

group1 = file.get('output')
if group1==None:
    print "No group 'output' exists in selected HDF5 file. Maybe not yet used to create data? Now exiting."
    sys.exit()

sol = group1.get('solution')
if sol==None:
    print "Could not find dataset 'solution' in selected HDF5 file. Now exiting."
    sys.exit()
sol = sol[...]
sol = sol.reshape(2800, Nx, 3)

pylab.plot(xcenter, sol[-1,:,0])
pylab.show()
