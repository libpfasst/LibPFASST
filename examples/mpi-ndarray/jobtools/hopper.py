"""Hopper scheduler (Cray PBS)."""

from StringIO import StringIO as strio

import base

from fabric.api import *


class HopperPBS(base.Container):

    def submit(self, job, rwd=None, exe=None, inputs=None, queue=None,
               width=None, depth=None, pernode=None,
               walltime=None, stdout=None, stderr=None,
               pbs_opts=None, aprun_opts=None, verbose=False, dry_run=None, 
               cmds=None, **kwargs):

        #
        # import defaults from fabric env
        #

        if width is None:
            width = getattr(env, 'width', None)

        if depth is None:
            depth = getattr(env, 'depth', None)

        if pernode is None:
            pernode = getattr(env, 'pernode', None)

        if queue is None:
            queue = getattr(env, 'queue', None)

        if pbs_opts is None:
            pbs_opts = getattr(env, 'pbs_opts', None)

        if aprun_opts is None:
            aprun_opts = getattr(env, 'aprun_opts', None)

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
        # construct pbs script
        #

        pbs = [ '#!/bin/sh' ]
        pbs.append("#PBS -N " + job.name)
        if queue:
            pbs.append("#PBS -q " + queue)
        if width:
            pbs.append("#PBS -l mppwidth=" + str(width))
        if depth:
            pbs.append("#PBS -l mppdepth=" + str(depth))
        if pernode:
            pbs.append("#PBS -l mppnppn=" + str(pernode))
        if walltime:
            pbs.append("#PBS -l walltime=" + walltime)
        if stdout:
            pbs.append("#PBS -o {rwd}/" + stdout)
        if stderr:
            pbs.append("#PBS -e {rwd}/" + stderr)
        if pbs_opts:
            pbs.extend(pbs_opts)
        pbs.append("#PBS -V")

        pbs.append("")
        pbs.append("cd {rwd}")
        if cmds:
            pbs.extend(cmds)
        if depth:
            pbs.append("export OMP_NUM_THREADS=" + str(depth))
        pbs.append("aprun {opts} -B {exe} {inputs}") # change this

        if aprun_opts:
            opts = ' '.join(aprun_opts)
        else:
            opts = ''

        pbs = '\n'.join(pbs).format(opts=opts, rwd=rwd, exe=exe, inputs=inputs)


        if not dry_run:
            # push to remote host
            run_script = rwd + '/pbs.sh'
            put(strio(pbs), run_script)
            
            # submit to queue
            run('qsub ' + run_script)

        else:
            print pbs


