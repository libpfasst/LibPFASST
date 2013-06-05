"""Hopper scheduler (Cray PBS)."""

import base

from fabric.api import *


class Serial(base.Container):

    def submit(self, job, rwd=None, exe=None, inputs=None,
               width=None, depth=None, stdout=None, stderr=None,
               verbose=False, dry_run=None, **kwargs):

        #
        # import defaults from fabric env
        #

        if width is None:
            width = getattr(env, 'width', None)

        if depth is None:
            depth = getattr(env, 'depth', None)

        if stdout is None:
            stdout = getattr(env, 'stdout', 'stdout')

        if stderr is None:
            stderr = getattr(env, 'stderr', 'stderr')

        if dry_run is None:
            dry_run = getattr(env, 'dry_run', False)

        verbose = verbose or getattr(env, 'verbose', False)

        mpirun = getattr(env, 'mpirun', 'mpirun')

        #
        # sanity checks
        #

        if not rwd:
            raise ValueError("Invalid remote working directory (rwd): not set.")

        if not exe:
            raise ValueError("Invalid executable (exe): not set.")

        #
        # construct run script
        #

        sh = [ '#!/bin/sh' ]
        sh.append("")
        sh.append("cd {rwd}")
        if depth:
            sh.append("export OMP_NUM_THREADS={depth}")
        sh.append("{mpirun} -n {width} {exe} {inputs} 2> {stderr} | tee {stdout}")

        sh = '\n'.join(sh).format(rwd=rwd, exe=exe, width=width, depth=depth, inputs=inputs,
                                  mpirun=mpirun, stderr=stderr, stdout=stdout)

        # push to remote host
        run_script = rwd + '/run.sh'

        # launch
        if dry_run:
            print sh
        else:
            import StringIO
            put(StringIO.StringIO(sh), run_script)
            run('sh ' + run_script)
