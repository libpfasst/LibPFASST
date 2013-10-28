
import subprocess
import numpy as np
import tempfile
import shutil
from cStringIO import StringIO

import pf

with open('gpanim.gnu', 'r') as f:
    gp = f.read()

dname = tempfile.mkdtemp(dir='.')

for i, sol in enumerate(pf.io.read_final('test')):
    data = StringIO()
    np.savetxt(data, np.load(sol.fname)[::4])

    src = gp.format(
        terminal='pngcairo',
        output='%s/f%05d.png' % (dname, i),
        data=data.getvalue())

    p = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
    p.communicate(src)

subprocess.call(['ffmpeg -i %s/f%%05d.png -vcodec mpeg4 anim.avi' % dname], shell=True)
    
shutil.rmtree(dname)
    
    
