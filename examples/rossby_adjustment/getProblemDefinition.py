import math
import getMesh
import numpy
from collections import defaultdict
def get():
    Nx     = 600
    Ny     = 1
    xleft  = -15.0
    xright = 15.0
    ydown  = 0.0
    yup    = 1.0
    domain_height = abs(yup - ydown)
    domain_width  = abs(xright - xleft)
    Nsteps        = 2800
    Tend          = 10.0
    Niter         = 5
    grav          = 1.0
    coriolisPar   = 0.0

   # define function that provides initial value
    def q0(x,y):
#            return math.sin(2.0*math.pi*x)*math.sin(2.0*math.pi*y)
        if x < 0:
            return 1.1
        else:
            return 0.9
        
    # if this is an actual 2D problem, get a 2D mesh and evaluate q0 accordingly
    if Ny>1:
        XX, YY, dx, dy = getMesh.getMesh(Nx,Ny,xleft,xright,yup,ydown)
    
    # evaluate q0 on generate mesh --> adapt for IV with multiple components
        Q = 0.0*XX
        for ii in range(0,Nx):
            for jj in range(0,Ny):
                Q[jj,ii] = q0(XX[jj,ii],YY[jj,ii])

    else:
        xnodes  = numpy.linspace(xleft, xright, Nx+1)
        xcenter = [0 for x in range(Nx)]
        for ii in range(0,Nx):
            xcenter[ii] = 0.5*( xnodes[ii] + xnodes[ii+1])
# Note that the order here is reversed compared to how it is read in the FORTRAN code. Apparently, h5py also
# assumes a C-type ordering, so the "wrong" sorting here gets written in a order that is correct when read by
# FORTRAN HDF5
            Q = numpy.empty([ Ny, Nx, 3])
            for kk in range(0,3):
                for ii in range(0,Nx):
                    if kk==0:
                        Q[0,ii,kk] = q0(xcenter[ii],0.0)
                    else:
                        Q[0,ii,kk] = 0.0
    
    dx = (xcenter[1] - xcenter[0])
    dy = 0.0

    # Initialize 
    problem = {}
    problem['input'] = {}
    problem['input']['problemdefinition'] = defaultdict(list)
    problem['input']['problemdefinition']['nr_fields']     = ['i', 3, (1,1)]
    problem['input']['problemdefinition']['sound_speed']   = ['d', 0.0, (1,1)]
    problem['input']['problemdefinition']['stabFreq']      = ['d', 0.0, (1,1)]
    problem['input']['problemdefinition']['grav']          = ['d', grav, (1,1)]
    problem['input']['problemdefinition']['coriolisPar']   = ['d', coriolisPar, (1,1)]
    problem['input']['problemdefinition']['domain_height'] = ['d', domain_height, (1,1)]
    problem['input']['problemdefinition']['domain_width']  = ['d', domain_width, (1,1)] 
    problem['input']['problemdefinition']['BC']            = ['i', 2, (1,1)]
    problem['input']['problemdefinition']['x_left']        = ['d', xleft, (1,1)]
    problem['input']['problemdefinition']['x_right']       = ['d', xright, (1,1)]
    problem['input']['problemdefinition']['y_up']          = ['d', yup, (1,1)]
    problem['input']['problemdefinition']['y_down']        = ['d', ydown, (1,1)]
    problem['input']['problemdefinition']['q_initial']     = ['d', Q, (Ny,Nx, 3)]
    problem['input']['problemdefinition']['Uadv']          = ['d', 0.0, (1, 1)]
    problem['input']['problemdefinition']['Vadv']          = ['d', 0.0, (1, 1)]

    problem['input']['discretization'] = {}
    problem['input']['discretization']['order_coarse_advection'] = ['i', 5, (1,1)]
    problem['input']['discretization']['order_coarse_sound']     = ['i', 4, (1,1)]
    problem['input']['discretization']['order_fine_advection']   = ['i', 5, (1,1)]
    problem['input']['discretization']['order_fine_sound']       = ['i', 4, (1,1)]
    problem['input']['discretization']['Nx']                     = ['i', Nx, (1,1)]
    problem['input']['discretization']['Ny']                     = ['i', Ny, (1,1)]
    problem['input']['discretization']['Nx_coarse']              = ['i', Nx, (1,1)]
    problem['input']['discretization']['Ny_coarse']              = ['i', Ny, (1,1)]
    problem['input']['discretization']['dx']                     = ['d', dx, (1,1)]
    problem['input']['discretization']['dy']                     = ['d', dy, (1,1)]
    problem['input']['discretization']['dx_coarse']              = ['d', dx, (1,1)]
    problem['input']['discretization']['dy_coarse']              = ['d', dy, (1,1)]

    problem['input']['integration'] = {}
    problem['input']['integration']['global_tend']             = ['d', Tend, (1,1)]
    problem['input']['integration']['global_tstart']           = ['d', 0.0, (1,1)]
    problem['input']['integration']['nu_coarse']               = ['d', 0.02, (1,1)]
    problem['input']['integration']['nu']                      = ['d', 0.02, (1,1)]
    problem['input']['integration']['Nsteps_fine_total']       = ['i', Nsteps, (1,1)]
    problem['input']['integration']['Nsteps_coarse_total']     = ['i', Nsteps, (1,1)]
    problem['input']['integration']['Nsteps_sound_coarse']     = ['i', -1, (1,1)]
    problem['input']['integration']['Nsteps_sound_fine']       = ['i', -1, (1,1)]

    problem['input']['parareal'] = {}
    problem['input']['parareal']['maxit']    = ['i', Niter, (1,1)]
    problem['input']['parareal']['tolit']    = ['d', 1e-10, (1,1)]
    problem['input']['parareal']['nthreads'] = ['i', 16, (1,1)]

    Nc = problem['input']['integration']['Nsteps_coarse_total'][1]
    Nt = problem['input']['parareal']['nthreads'][1]
    Nc_per_slice = Nc/Nt
    problem['input']['integration']['Nsteps_coarse_per_slice'] = ['i', Nc_per_slice, (1,1)]

    problem['input']['topology'] = {}
    problem['input']['topology']['nprocs_x'] = ['i', 1, (1,1)]
    problem['input']['topology']['nprocs_y'] = ['i', 1, (1,1)]
    return problem
