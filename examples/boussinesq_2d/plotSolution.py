import h5py
import sys
import getMesh
import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, animation, rcParams
import matplotlib.pyplot as plt
import pylab

filename='test.h5'
file = h5py.File(filename,'r')

nr_fields = (file.get('input/problemdefinition/nr_fields'))[()][0][0]
Nx      = (file.get('input/discretization/Nx'))[()][0][0]
Ny      = (file.get('input/discretization/Ny'))[()][0][0]
xleft   = (file.get('input/problemdefinition/x_left'))[()][0][0]
xright  = (file.get('input/problemdefinition/x_right'))[()][0][0]
yup     = (file.get('input/problemdefinition/y_up'))[()][0][0]
ydown   = (file.get('input/problemdefinition/y_down'))[()][0][0]

XX, YY, dx, dy = getMesh.getMesh(Nx,Ny,xleft,xright,yup,ydown)

group1 = file.get('output')
if group1==None:
    print "No group 'output' exists in selected HDF5 file. Maybe not yet used to create data? Now exiting."
    sys.exit()

sol = group1.get('solution')
if sol==None:
    print "Could not find dataset 'solution' in selected HDF5 file. Now exiting."
    sys.exit()
Ndumps = sol.shape[0]
sol    = sol[...]

print sol[-1,:,:,2].min()
print sol[-1,:,:,2].max()
#fig = [None for x in range(0,nr_fields)]
#ax  = [None for x in range(0,nr_fields)]
#for ii in range(0,nr_fields):
#    fig[ii] = plt.figure()
#    ax[ii] = fig[ii].gca(projection='3d')
#    ax[ii].plot_surface(XX,YY,sol[-1,:,:,0].T,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0)
#    ax[ii].set_title(ii+1)
    #ax[ii].view_init(30,60)

fig = plt.figure()
N = int((1e-1-5e-4)/5e-4)
contour_lines = numpy.linspace(5e-4,1e-1,N+1)
#contours = numpy.empty(2*N)
#for ii in range(0,N):
#    contours[ii]   = -contour_lines[N-ii-1]
#    contours[N+ii] = contour_lines[ii]

contours = [-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
CS = plt.contour(XX,YY,1e3*sol[-1,:,:,2].T, contours, colors='k')
#plt.clabel(CS, inline=1, fontsize=10)

plt.show()



