
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
            job.write_params(os.path.join(self.stage, job.rwd, 'probin.nml'))


    def rsync_stage(self):
        local('rsync -auz {src}/* {host}:{dst}/'.format(src=self.stage, host=env.host, dst=self.rwd))
        shutil.rmtree(self.stage)


    def submit_all(self):
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

            rwd = rwd + '/' + job.rwd

            kwargs = {}
            kwargs.update(job.attrs)
            kwargs.update(self.attrs)
            kwargs.pop('rwd', None)
            kwargs.pop('exe', None)

            # print some message...
            scheduler.submit(job, rwd=rwd, exe=exe, inputs='probin.nml', **kwargs)
