import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from subprocess import call
from pf.pfasst import PFASST, Experiment, Params

def init():
    lines[0].set_data([], [])
    lines[1].set_data([], [])
    lines[2].set_data([], [])

    return lines

def animate(i):
#    for j, line in enumerate(lines):
#        line.set_data(np.sin(q[i:i+10, j]), np.cos(q[i:i+10, j]))

    thisx = [0.75*np.sin(q[i, j])-0.02 for j in range(nparticles)]
    thisy = [0.75*np.cos(q[i, j]) for j in range(nparticles)]
    lines[0].set_data(thisx, thisy)

    histx.append([0.75*np.sin(q[i, j])-0.02 for j in [0]])#, nparticles-1]])
    histy.append([0.75*np.cos(q[i, j]) for j in [0]])#, nparticles-1]])
    lines[1].set_data(histx, histy)

#    histx.append([0.5*np.sin(q[i, j])-0.05 for j in [0]])#, nparticles-1]])
#    histy.append([0.5*np.cos(q[i, j])-0.05 for j in [0]])#, nparticles-1]])
    lines[2].set_data(histx, histy)
    return lines
#    return lines[0], lines[1], lines[2], lines[3], \
#            lines[4], lines[5], lines[6], lines[7], \
#            lines[8], lines[9], lines[10], time_text



if __name__ == '__main__':
    exe = '/home/bkrull/devel/pfasst-nwchem/libpfasst/tests/toda/main.exe'
    exp = Experiment()
    cmap = plt.cm.Dark2.colors*2

    nsteps = 4096
    nparticles = 11
    periodic = True
    tfinal = 25.0

    dt = tfinal/nsteps

    params = Params(tfinal=tfinal, nodes=[3], magnus=[2], tasks=8,
                    sweeps_pred=[2], 
                    particles=nparticles, periodic=periodic, tolerance=1e-10,
                    iterations=20, nsteps=nsteps, solutions=True)
    toda = PFASST(exe, params)

    results = toda.run()[0]

    q, p = toda.get_toda_solutions(results)
    
    fig, ax = plt.subplots(figsize=(3,3), dpi=300)
    lines = [ax.plot([],[], marker='o', color=cmap[0])[0], 
             ax.plot([],[], lw=1, color='red')[0],
             ax.plot([],[], lw=1, color='black')[0]]
    histx = []
    histy = []

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.axis('off')
    ax.grid();

    ani = animation.FuncAnimation(fig, animate, np.arange(0, nsteps), \
            interval=5, blit=True, init_func=init)

    ani.save('input.mp4', fps=60)
    #plt.show()

    #convert_to_gif = 'ffmpeg -i input.mp4 -r 10 -f image2pipe -vcodec ppm - | ' + 'convert -delay 5 -loop 0 - output.gif'

    #convert_to_gif = convert_to_gif.split()

    #shrink_gif = 'gifsicle output.gif --resize 300x300 > output300.gif'.split()

    #call(convert_to_gif)
    #call(shrink_gif, shell=True)
