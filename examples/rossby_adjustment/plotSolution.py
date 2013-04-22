import h5py
import sys
import getMesh
import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

filename='test.h5'
file = h5py.File(filename,'r')
group1 = file.get('output')
Nx_dset = file.get('input/discretization/Nx')
# corresponding dataset for some reasons is an array in an array; [0][0] retrieves integer
Nx           = Nx_dset[()][0][0]
Ny_dset      = file.get('input/discretization/Ny')
Ny           = Ny_dset[()][0][0]
xleft_dset   = file.get('input/problemdefinition/x_left')
x_left       = xleft_dset[()][0][0]
xright_dset  = file.get('input/problemdefinition/x_right')
x_right      = xright_dset[()][0][0]
yup_dset     = file.get('input/problemdefinition/y_up')
y_up         = yup_dset[()][0][0]
ydown_dset   = file.get('input/problemdefinition/y_down')
y_down       = ydown_dset[()][0][0]

XX, YY, dx, dy = getMesh.getMesh(Nx, Ny, x_left, x_right, y_up, y_down)

if group1==None:
    print "No group 'output' exists in selected HDF5 file. Maybe not yet used to create data? Now exiting."
    sys.exit()

sol = group1.get('solution')
if sol==None:
    print "Could not find dataset 'solution' in selected HDF5 file. Now exiting."
    sys.exit()
sol = sol[...]
sol = sol.reshape(Ny, Nx)

fig = plt.figure()
ax  = Axes3D(fig)
N  = sol/sol.max()
surf = ax.plot_surface(XX, YY, sol, rstride=1, cstride=1, cmap='hot')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
