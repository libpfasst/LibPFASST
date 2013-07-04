"""JUQUEEN scheduler (IBM LoadLeveler)."""

from StringIO import StringIO as strio

import base

from fabric.api import *


class JQLL(base.Container):

    def submit(self, job, rwd=None, exe=None, inputs=None, queue=None,
               width=None, depth=None, pernode=None,
               walltime=None, stdout=None, stderr=None,
               verbose=False, dry_run=None, 
               cmds=None, **kwargs):

        #
        # import defaults from fabric env
        #

        if width is None:
            width = getattr(env, 'width', None)

        if depth is None:
            depth = getattr(env, 'depth', None)

        if pernode is None:
            pernode = getattr(env, 'pernode', 1)

        if queue is None:
            queue = getattr(env, 'queue', None)

        if walltime is None:
            walltime = getattr(env, 'walltime', None)

        if stdout is None:
            stdout = getattr(env, 'stdout', 'stdout')

        if stderr is None:
            stderr = getattr(env, 'stderr', 'stderr')

        verbose = verbose or getattr(env, 'verbose', False)

        #
        # sanity checks
        #

        if not rwd:
            raise ValueError("Invalid remote working directory (rwd): not set.")

        if not exe:
            raise ValueError("Invalid executable (exe): not set.")

        #
        # construct ll script
        #

        ll = []
        ll.append("# @ job_name = " + job.name)
        # if queue:
        #     ll.append("# @ -q " + queue)
        if width:
            ll.append("# @ bg_size = " + str(width))
        # if depth:
        #     ll.append("# @ -l mppdepth=" + str(depth))
        if walltime:
            ll.append("# @ wall_clock_limit = " + walltime)
        if stdout:
            ll.append("# @ output = {rwd}/" + stdout)
        if stderr:
            ll.append("# @ error = {rwd}/" + stderr)
        
        ll.append('# @ notification = always')
        ll.append('# @ notify_user = memmett@gmail.com')
        ll.append('# @ environment = COPY_ALL')

        ll.append('# @ job_type = bluegene')
        ll.append('# @ queue')

        ll.append("")
        ll.append("cd {rwd}")
        if cmds:
            ll.extend(cmds)

        # if depth:
        #     ll.append("export OMP_NUM_THREADS=" + str(depth))

        runjob = [ 'runjob' ]
        if pernode:
            runjob.append('--ranks-per-node ' + str(pernode))
        if width < 32:
            runjob.append('--np ' + str(width))
        runjob.append(' : {exe} {inputs}')
        ll.append(' '.join(runjob))
        ll.append('')

        ll = '\n'.join(ll).format(rwd=rwd, exe=exe, inputs=inputs)

        if not dry_run:
            run_script = rwd + '/ll.sh'
            put(strio(ll), run_script)
            run('llsubmit ' + run_script)

        else:
            print ll
