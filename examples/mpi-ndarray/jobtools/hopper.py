"""Hopper scheduler (Cray PBS)."""

from StringIO import StringIO as strio

from base import Container

from fabric.api import *


class HopperPBS(object):

    def submit(self, **kwargs):

        a = Container(**kwargs)

        if not a.rwd:
            raise ValueError("Invalid remote working directory (rwd): not set.")

        if not a.exe:
            raise ValueError("Invalid executable (exe): not set.")


        pbs = [ '#!/bin/sh' ]
        pbs.append("#PBS -N " + a.name)
        if a.queue:
            pbs.append("#PBS -q " + a.queue)
        if a.pernode:
            mppwidth = a.width * 24 / a.pernode
        else:
            mppwidth = a.width * a.get('depth', 1)
        pbs.append("#PBS -l mppwidth=" + str(mppwidth))
        if a.walltime:
            pbs.append("#PBS -l walltime=" + a.walltime)
        if a.stdout:
            pbs.append("#PBS -o {rwd}/" + a.stdout)
        else:
            pbs.append("#PBS -o {rwd}/stdout")
        if a.stderr:
            pbs.append("#PBS -e {rwd}/" + a.stderr)
        else:
            pbs.append("#PBS -e {rwd}/stderr")
        if a.pbs_opts:
            pbs.extend(a.pbs_opts)
        pbs.append("#PBS -V")

        pbs.append("")
        pbs.append("cd {rwd}")
        if a.depth:
            pbs.append("export OMP_NUM_THREADS=" + str(a.depth))

        if a.pbs_cmds:
            pbs.extend(a.pbs_cmds)
        pbs.append("aprun {opts} {exe} {inputs} {cmd_opts}")

        opts = []
        if a.width:
            opts.append("-n " + str(a.width))
        if a.depth:
            opts.append('-d ' + str(a.depth))
        if a.pernode:
            opts.append('-N ' + str(a.pernode))
            # opts.append('-S ' + str(a.pernode/2))
            opts.append('-cc numa_node')
        if a.specialized:
            opts.append('-r ' + str(a.specialized))
            if a.pernode is None:
                opts.append('-cc numa_node')
        if a.aprun_opts:
            opts.extend(a.aprun_opts)

        pbs = '\n'.join(pbs).format(opts=' '.join(opts), rwd=a.rwd, exe=a.exe, inputs=a.inputs, cmd_opts=a.cmd_opts or '')

        if not a.dry_run:
            run_script = a.rwd + '/pbs.sh'
            put(strio(pbs), run_script)
            run('qsub ' + run_script)
        else:
            print pbs


