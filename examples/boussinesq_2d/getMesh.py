import numpy
def getMesh(Nx,Ny,xleft,xright,yup,ydown):
    xnodes  = numpy.linspace(xleft, xright, Nx+1)
    xcenter = [0 for x in range(Nx)]
    for ii in range(0,Nx):
        xcenter[ii] = 0.5*( xnodes[ii] + xnodes[ii+1])

        ynodes = numpy.linspace(ydown, yup, Ny+1)
        ycenter = [0 for x in range(Ny)]
        
        for jj in range(0,Ny):
            ycenter[jj] = 0.5*( ynodes[jj] + ynodes[jj+1])
            
            XX, YY = numpy.meshgrid(xcenter, ycenter)
            YY = numpy.flipud(YY)
    dx = xcenter[1]-xcenter[0]                    
    dy = ycenter[1]-ycenter[0]
    return XX, YY, dx, dy

