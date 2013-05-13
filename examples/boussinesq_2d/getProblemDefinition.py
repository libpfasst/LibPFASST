import math
import getMesh
import numpy
from collections import defaultdict
# Formulas for compution nondimensional parameter:
#
# Given : reference length L and time T such that U = L/T = 20 m/s
# stabFreq^2 = L*N^2/g with N = 0.01 1/s and g = 10 m/s^2
# grav       = T^2*g/L
# c_s        = T*c_s_dim/L with c_s_dim = 300m/s
# xright     = 300km/L
# yup        = 10km/L
# Tend       = 3000s/T
# dt         = 1s/T or 2s/T
# Nsteps     = Tend/dt
# dx=dz      = 1km/L , compute Nx, Ny accordingly
# a          = 5km/L

# Use: L = 10km and T = 500s ==> U = 10,000m/500s = 20 m/s
# ...and stabFreq = sqrt(0.1), grav = 250, c_s = 15
def get():
    nr_fields = 4
    Nx        = 300
    Ny        = 60
    xleft     = 0.0
    xright    = 30.0
    ydown     = 0.0
    yup       = 1.0
    Nsteps    = 6000
    Tend      = 6.0
    Niter     = 6
    grav      = 250.0
    stabFreq  = math.sqrt(0.1)
    c_s       = 15.0

   # define function that provides initial value
    x_c = 10
    a   = 0.5
    def q0(x,y):
            #return math.sin(2.0*math.pi*x)*math.sin(2.0*math.pi*y)
        #return math.exp(-(x-x_c)**2/3.0**2 - (y-0.5)**2/0.1**2)
        return 0.01*math.sin(math.pi*y/yup)/( 1 + (x - x_c )**2/a**2 )

    # if this is an actual 2D problem, get a 2D mesh and evaluate q0 accordingly
    XX, YY, dx, dy = getMesh.getMesh(Nx,Ny,xleft,xright,yup,ydown)
    
    # evaluate q0 on generate mesh --> adapt for IV with multiple components
    Q = numpy.empty([Nx,Ny,nr_fields])
    for jj in range(0,Ny):
        for ii in range(0,Nx):
            for kk in range(0,nr_fields):
                if kk==2:
                    Q[ii,jj,kk] = q0(XX[jj,ii],YY[jj,ii])
                else:
                    Q[ii,jj,kk] = 0.0

    domain_height = abs(yup - ydown)
    domain_width  = abs(xright - xleft)


    # Initialize 
    problem = {}
    problem['input'] = {}
    problem['input']['problemdefinition'] = defaultdict(list)
    problem['input']['problemdefinition']['nr_fields']     = ['i', nr_fields, (1,1)]
    problem['input']['problemdefinition']['sound_speed']   = ['d', c_s, (1,1)]
    problem['input']['problemdefinition']['stabFreq']      = ['d', stabFreq, (1,1)]
    problem['input']['problemdefinition']['grav']          = ['d', grav, (1,1)]
    problem['input']['problemdefinition']['coriolisPar']   = ['d', 0.0, (1,1)]
    problem['input']['problemdefinition']['domain_height'] = ['d', domain_height, (1,1)]
    problem['input']['problemdefinition']['domain_width']  = ['d', domain_width, (1,1)] 
    problem['input']['problemdefinition']['BC']            = ['i', 3, (1,1)]
    problem['input']['problemdefinition']['x_left']        = ['d', xleft, (1,1)]
    problem['input']['problemdefinition']['x_right']       = ['d', xright, (1,1)]
    problem['input']['problemdefinition']['y_up']          = ['d', yup, (1,1)]
    problem['input']['problemdefinition']['y_down']        = ['d', ydown, (1,1)]
    problem['input']['problemdefinition']['q_initial']     = ['d', Q, (Nx, Ny, nr_fields)]
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
    problem['input']['integration']['nu_coarse']               = ['d', 0.0, (1,1)]
    problem['input']['integration']['nu']                      = ['d', 0.0, (1,1)]
    problem['input']['integration']['Nsteps_fine_total']       = ['i', Nsteps, (1,1)]
    problem['input']['integration']['Nsteps_coarse_total']     = ['i', Nsteps, (1,1)]
    problem['input']['integration']['Nsteps_sound_coarse']     = ['i', -1, (1,1)]
    problem['input']['integration']['Nsteps_sound_fine']       = ['i', -1, (1,1)]

    problem['input']['parareal'] = {}
    problem['input']['parareal']['maxit']    = ['i', Niter, (1,1)]
    problem['input']['parareal']['tolit']    = ['d', 1e-10, (1,1)]
    problem['input']['parareal']['nthreads'] = ['i', 1, (1,1)]

    Nc = problem['input']['integration']['Nsteps_coarse_total'][1]
    Nt = problem['input']['parareal']['nthreads'][1]
    Nc_per_slice = Nc/Nt
    problem['input']['integration']['Nsteps_coarse_per_slice'] = ['i', Nc_per_slice, (1,1)]

    problem['input']['topology'] = {}
    problem['input']['topology']['nprocs_x'] = ['i', 1, (1,1)]
    problem['input']['topology']['nprocs_y'] = ['i', 1, (1,1)]
    return problem
