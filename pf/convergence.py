
from itertools import product

import numpy as np
import pylab
import pf.io


class Cache(object):
  '''Simple cache for solutions.'''
  def __init__(self):
    self.cache = {}
  def get(self, name):
    if name not in self.cache:
      self.cache[name] = np.load(name)
    return self.cache[name]


def step_iter_level_map(available):
  return { (x.step, x.iter, x.level): x.fname for x in available }


def find(available, step, iter, level):
  singleton = [ x.fname for x in available 
                if  x.step  == step
                and x.iter  == iter
                and x.level == level ]

  if len(singleton) == 1:
    return singleton[0]

  raise ValueError('more than one solution found')


def errors(reference, approximate, **kwargs):
  '''Compute errors of approximate solutions relative to the reference solution.

  :param reference:   list of available reference solutions
  :param approximate: list of available approximate solutions

  Basic usage:

  >>> reference = pf.io.read_avail('ref_output_dir')
  >>> approximate = pf.io.read_avail('app_output_dir')
  >>> errors, steps, iters, levels = pf.convergence.errors(reference, approximate)

  '''

  cache = Cache()

  steps  = sorted(list(set([ x.step for x in approximate ])))
  iters  = sorted(list(set([ x.iter for x in approximate ])))
  levels = sorted(list(set([ x.level for x in approximate ])))
  riter  = max([ x.iter for x in reference ])
  rlev   = max([ x.level for x in reference ])
  rmap   = step_iter_level_map(reference)
  amap   = step_iter_level_map(approximate)

  trat   = max([ x.step for x in reference ]) / max([ x.step for x in approximate ])

  errors = {}
  for step, aiter, alev in product(steps, iters, levels):
    ref = cache.get(rmap[step*trat, riter, rlev])
    app = cache.get(amap[step, aiter, alev])
    err = abs(ref - app).max()
    errors[step, aiter, alev] = err

  return errors, steps, iters, levels


def plot(errs, steps, iters, levels, **kwargs):
  '''Plot error vs time (for each iteration) of approximate solutions
  relative to the reference solution.

  :param errs:   dictionary of errors
  :param steps:  list of steps
  :param iters:  list of iterations
  :param level:  list of levels

  Basic usage:

  >>> reference   = pf.io.read_avail('ref_output_dir')
  >>> approximate = pf.io.read_avail('app_output_dir')
  >>> errs, steps, iters, levels = pf.convergence.errors(reference, approximate)
  >>> pf.convergence.plot(errs, steps, iters, levels)

  '''

  fig, ax = pylab.subplots(ncols=len(levels))

  if not isinstance(ax, list):
    ax = [ ax ]
  
  for l, level in enumerate(levels):
    for i, iter in enumerate(iters):

      x = steps
      y = [ errs[step,iter,level] for step in steps ]

      ax[l].semilogy(x, y, **kwargs)
      ax[l].set_title('level %d' % level)
      ax[l].set_xlabel('step/processor')
      ax[l].set_ylabel('max abs. error')

  return fig, ax



  

  
    
    
