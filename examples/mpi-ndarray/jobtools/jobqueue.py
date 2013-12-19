
import os
import shutil

import base
import schedulers

from fabric.api import *


class JobQueue(base.Container):

    def __init__(self, **jqargs):
        base.Container.__init__(self, **jqargs)
        self.jobs = []
        self.stage = 'staging'


    def add(self, job):
        """Add *job* to the queue.

        Job resource (specific to the queueing system in use) can be
        specified as keyword arguments.
        """
        self.jobs.append(job)


    def stage_all(self):
        for job in self.jobs:
            os.makedirs(os.path.join(self.stage, job.rwd))
            if job.has_params:
                job.write_params(os.path.join(self.stage, job.rwd, 'probin.nml'))


    def rsync_stage(self):
        run('mkdir -p {dst}'.format(dst=self.rwd))
        local('rsync -auz {src}/* {host}:{dst}/'.format(src=self.stage, host=env.host, dst=self.rwd))
        shutil.rmtree(self.stage)


    def submit_all(self, **kwargs):
        self.stage_all()
        self.rsync_stage()

        scheduler = self.scheduler or env.get('scheduler', None)
        if isinstance(scheduler, str):
            Scheduler = schedulers.by_name[scheduler]
            scheduler = Scheduler()

        for job in self.jobs:
            rwd = self.rwd or env.get('rwd', '')
            exe = self.exe or job.exe or env.get('exe', None)
            if not exe.startswith('/'):
                exe = rwd + '/' + exe

            if not job.rwd.startswith('/'):
                rwd = rwd + '/' + job.rwd
            else:
                rwd = job.rwd

            kwargs.update(env)
            kwargs.update(self)
            kwargs.update(job)
            kwargs.pop('rwd', None)
            kwargs.pop('exe', None)

            if job.has_params:
                scheduler.submit(exe=exe, rwd=rwd, inputs='probin.nml', **kwargs)
            else:
                scheduler.submit(exe=exe, rwd=rwd, inputs='', **kwargs)
