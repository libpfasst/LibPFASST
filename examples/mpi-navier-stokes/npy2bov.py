
import sys
import numpy as np

for fname in sys.argv[1:]:
    z = np.load(fname)
    u = np.real(np.fft.ifftn(z[:,:,:,0]))
    v = np.real(np.fft.ifftn(z[:,:,:,1]))
    w = np.real(np.fft.ifftn(z[:,:,:,2]))

    basename = fname.split('.')[0]
    dat = basename + '.dat'
    bov = basename + '.bov'

    n = u.shape[0]

    print 'converting:', fname

    b = np.empty((3,) + u.shape)
    b[0] = u
    b[1] = v
    b[2] = w
    with open(dat, 'wb') as f:
        f.write(b)

    with open(bov, 'w') as f:
        f.write('DATA_FILE: %s\n' % dat)
        f.write('DATA_SIZE: %d %d %d\n' % (n, n, n))
        f.write('DATA_COMPONENTS: 3\n')


