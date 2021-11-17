"""Wrapper classes for running pfasst.exe executables

Classes
======
PFASST     : contains all setup/run information for a PFASST calculation
             PFASST.run() returns a trajectory (a dataframe of all states of calculation)
Experiment : mostly empty class with methods that take PFASST as input and
             perform a given 'experiment'
             Experiments return a Results object, a dataframe-derived class that has
             additional plot functionality
Params     : stores all PFASST parameters, mostly implemented to reduce number of
             PFASST class attributes
Results    : pd.DataFrame child that implements pfasst-related plots, and results
             storing
"""
import glob
import re
import attr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
from os import remove, mkdir
from subprocess import check_output, STDOUT, CalledProcessError
from scipy.io import FortranFile
from pf.io import read_all_timings, Timing

try:
    termtype = get_ipython().__class__.__name__
except:
    from tqdm import tqdm as tqdm
else:
    if 'ZMQ' in termtype:
        from tqdm import tqdm_notebook as tqdm
    else:
        from tqdm import tqdm as tqdm

@attr.s(slots=True)
class Params(object):
    """Class containing all sweeper-independent parameters necessary for
    running a PFASST calculation.

    This is the main class that is ultimately passed around from all three of
    the other classes found in this file. Its structure/usability is very
    important, it's still in the works. Class-level assignment via attr is used
    in order to make access attribute style rather than dictionary.

    Sets defaults first, then, if `nb` is set to False (i.e for use in a
    shell), it'll then invoke an ArgumentParser to get command-line options,
    including the possibility to read in a file. These read-in-from-cli options
    are used to overwrite default parameter values.
    """
    #pylint: disable=too-many-instance-attributes
    nb = attr.ib(default=True, repr=False)
    exe = attr.ib(default=None)
    levels = attr.ib(default=1)
    tfinal = attr.ib(default=10.0)
    iterations = attr.ib(default=30)
    nsteps = attr.ib(default=16)
    nodes = attr.ib(default=[2], validator=attr.validators.instance_of(list))
    sweeps = attr.ib(default=[1], validator=attr.validators.instance_of(list))
    sweeps_pred = attr.ib(default=[1],
                          validator=attr.validators.instance_of(list))
    tasks = attr.ib(default=1)
    base_dir = attr.ib(default='output', repr=False)
    verbose = attr.ib(default=False, repr=False)
    nersc = attr.ib(default=False)
    dt = attr.ib(default=None)
    timings = attr.ib(default=False)
    vcycle = attr.ib(default=False)
    tolerance = attr.ib(default=1e-12)
    qtype = attr.ib(default='lobatto')
    sdc = attr.ib(default=True)

    def _init_dt(self):
        if self.dt is None:
            self.dt = self.tfinal / self.nsteps
        else:
            self.nsteps = self.tfinal / self.dt

    def print_params(self):
        params = {
            'tfinal': self.tfinal,
            'nsteps': self.nsteps,
            'dt': self.dt,
            'nodes': self.nodes,
        }
        pprint(params, width=1)

    def pack(self):
        return attr.asdict(self)

    def unpack(self, params):
        for k, v in params.items():
            setattr(self, k, v)


@attr.s(slots=True)
class MagpicardParams(Params):
    """ Classes that derive from Params are implementations of sweeper/problem-specific
    parameters. When building a new derived Params class one should do the following:
    1. Define new params (e.g. magnus, exptol, solutions, ...)
    2. Define new 'base_string' which is the parameter-less string that will be
       written to disk with the appropriate params when a calculation is run
    3. Define new 'pkl' string which is a string that uniquely identifies a previous
       run whose results were stored on disk (this will likely utilize the new params)
    4. Define method make_list() that constructs a list of the params that will be
       injected into base_strings before writing to disk
    5. Define method get_pkl_path() which constructs the string to/from which data
       will be saved/loaded
    6. Define ref_param_list which is a dictionary with k, v pairs that indicate the
       parameters for a reference calculation
    7. Implement method __attrs_post_init_() that does some final initialization that
       is required (a quirk of the attrs library)
    """
    #pylint: disable=too-many-instance-attributes
    filename = attr.ib(default='mag.nml')
    inttype = attr.ib(default='mag')
    magnus = attr.ib(default=[1], validator=attr.validators.instance_of(list))
    exptol = attr.ib(default=['1.d-15'],
                     validator=attr.validators.instance_of(list))
    solutions = attr.ib(default=False)
    particles = attr.ib(default=11)
    periodic = attr.ib(default=True)
    exptol = attr.ib(default=['1.d-15'],
                     validator=attr.validators.instance_of(list))
    base_string = attr.ib(default=
                          '&pf_params\n\tnlevels = {}\n\tniters = {}\n\t'+ \
                          'nnodes = {}\n\tnsweeps_pred = {}\n\tnsweeps = {}\n\t'+ \
                          'qtype = {}\n\techo_timings = {}\n\tabs_res_tol = {}\n\t'+ \
                          'rel_res_tol = {}\n\tpipeline_pred = .true.\n\t' + \
                          'pfasst_pred = .true.\n\tvcycle = {}\n/\n\n'+ \
                          '&params\n\tfbase = {}\n\t'+ \
                          'magnus_order = {}\n\ttfin = {}\n\tnsteps = {}\n\t'+ \
                          'exptol = {}\n\tnparticles = {}\n\t'+ \
                          'save_solutions = {}\n\ttoda_periodic = {}\n\t'+ \
                          'use_sdc = {}\n/\n')
    param_list = attr.ib(default=None)
    pkl = attr.ib(default=None)
    ref_param_list = attr.ib(default={
        'nodes': [3],
        'magnus': [3],
        'qtype': 'gauss'
    })

    def __attrs_post_init__(self):
        self.param_list = self._make_list()
        self.pkl = self.base_dir + \
                  '/tfinal_{}-dt_{}-'+ \
                  'particles_{}-periodic_{}-qtype-{}-'+ \
                  'levels_{}-nodes_{}-magorder_{}-tasks_{}-sdc_{}-inttype_{}.pkl'
        self._init_dt()

    def _make_list(self):
        nodes = ' '.join(map(str, self.nodes))
        magnus = ' '.join(map(str, self.magnus))
        sweeps = ' '.join(map(str, self.sweeps))
        sweeps_pred = ' '.join(map(str, self.sweeps_pred))
        exptol = ' '.join(self.exptol)
        basedir = '"{}"'.format(self.base_dir)

        if self.solutions == True:
            solns = '.true.'
        else:
            solns = '.false.'

        if self.timings == True:
            timings = '.true.'
        else:
            timings = '.false.'

        if self.periodic == True:
            periodic = '.true.'
        else:
            periodic = '.false.'

        if self.sdc == True:
            sdc = '.true.'
        else:
            sdc = '.false.'

        if self.vcycle == True:
            vcycle = '.true.'
        else:
            vcycle = '.false.'

        if self.qtype == 'gauss':
            qtype = 5
        else:
            qtype = 1

        return [self.levels, self.iterations,
                nodes, sweeps_pred, sweeps, qtype, timings,
                self.tolerance, self.tolerance, vcycle, basedir,
                magnus, self.tfinal,
                self.nsteps, exptol, self.particles,
                solns, periodic, sdc]

    def get_pkl_path(self):
        nodes = '_'.join(map(str, self.nodes))
        magnus = '_'.join(map(str, self.magnus))
        pkl_path = self.pkl.format(self.tfinal, self.dt,
                                   self.particles, self.periodic,
                                   self.qtype, self.levels,
                                   nodes, magnus, self.tasks, self.sdc, self.inttype)
        return pkl_path


@attr.s(slots=True)
class IMKParams(Params):
    filename = attr.ib(default='imk.nml')
    inttype = attr.ib(default='imk')
    exptol = attr.ib(default=['1.d-15'],
                     validator=attr.validators.instance_of(list))
    solutions = attr.ib(default=False)
    particles = attr.ib(default=11)
    periodic = attr.ib(default=True)
    nterms = attr.ib(default=[5], validator=attr.validators.instance_of(list))
    rk = attr.ib(default=False)
    mkrk = attr.ib(default=False)
    exptol = attr.ib(default=['1.d-15'],
                     validator=attr.validators.instance_of(list))
    base_string = attr.ib(default=
                          '&pf_params\n\tnlevels = {}\n\tniters = {}\n\t'+ \
                          'nnodes = {}\n\t' + \
                          'nsweeps_pred = {}\n\tnsweeps = {}\n\t'+ \
                          'qtype = {}\n\techo_timings = {}\n\tabs_res_tol = {}\n\t'+ \
                          'rel_res_tol = {}\n\tpipeline_pred = .true.\n\t' + \
                          'pfasst_pred = .true.\n\tvcycle = {}\n/\n\n'+ \
                          '&params\n\tfbase = {}\n\t'+ \
                          'nterms = {}\n\ttfin = {}\n\tnsteps = {}\n\t'+ \
                          'exptol = {}\n\tnparticles = {}\n\t'+ \
                          'save_solutions = {}\n\ttoda_periodic = {}\n\t'+ \
                          'use_sdc = {}\n\trk = {}\n\tmkrk = {}\n/\n')
    param_list = attr.ib(default=None)
    pkl = attr.ib(default=None)
    ref_param_list = attr.ib(default={
        'inttype': 'imk',
        'rk': False,
        'mkrk': False,
        'nodes': [9],
        'nterms': [20],
        'qtype': 'lob',
    })

    def __attrs_post_init__(self):
        pieces = self.inttype.split('.')
        if len(pieces) > 1:
            if 'mk' in pieces[-1]:
                self.mkrk = True
            else:
                self.rk = True
        self.param_list = self._make_list()
        self.pkl = self.base_dir + \
                  '/tfinal_{}-dt_{}-'+ \
                  'particles_{}-periodic_{}-qtype-{}-'+ \
                  'levels_{}-nodes_{}-nterms_{}-tasks_{}-sdc_{}-inttype_{}.pkl'
        self._init_dt()

    def _make_list(self):
        nodes = ' '.join(map(str, self.nodes))
        nterms = ' '.join(map(str, self.nterms))
        sweeps = ' '.join(map(str, self.sweeps))
        sweeps_pred = ' '.join(map(str, self.sweeps_pred))
        exptol = ' '.join(self.exptol)
        basedir = '"{}"'.format(self.base_dir)

        if self.solutions == True:
            solns = '.true.'
        else:
            solns = '.false.'

        if self.timings == True:
            timings = '.true.'
        else:
            timings = '.false.'

        if self.periodic == True:
            periodic = '.true.'
        else:
            periodic = '.false.'

        if self.sdc == True:
            sdc = '.true.'
        else:
            sdc = '.false.'

        if self.vcycle == True:
            vcycle = '.true.'
        else:
            vcycle = '.false.'

        if self.qtype == 'gauss':
            qtype = 5
        else:
            qtype = 1

        if self.rk == True:
            rk = '.true.'
            mkrk = '.false.'
        else:
            rk = '.false.'
            if self.mkrk == True:
                mkrk = '.true.'
            else:
                mkrk = '.false.'

        return [self.levels, self.iterations,
                nodes, sweeps_pred, sweeps,qtype, timings,
                self.tolerance, self.tolerance, vcycle, basedir,
                 nterms, self.tfinal,
                self.nsteps, exptol, self.particles,
                solns, periodic, sdc, rk, mkrk]

    def get_pkl_path(self):
        nodes = '_'.join(map(str, self.nodes))
        nterms = '_'.join(map(str, self.nterms))
        pkl_path = self.pkl.format(self.tfinal, self.dt,
                                   self.particles, self.periodic,
                                   self.qtype, self.levels,
                                   nodes, nterms, self.tasks, self.sdc, self.inttype)
        return pkl_path


class PFASST(object):
    """A wrapper class for compiled libpfasst programs.
    1. Get all the parameters in place
    2. (Optional) Change parameters (for experiments)
    3. Run
       a. Create the pfasst input file contents as str
       b. Write this string to file
       c. Check for saved pickle, if yes -> load, and return Results
       d. (If necessary) Run calculation via check_output
       e. (If necessary) Get and save trajectory with timings from Output

    Parameters
    ==========
    params : A Params class object is  either passed in, or a default created. Contains all
    necessary parameters for running a PFASST calculation.
    exe : a (str) path to the project's root e.g. '/exe/bkrull/apps/pfasst/dev' is used to
    path to binary

    Exposed Methods
    ===============
    write_pfstring_to_file : writes the input file for executable in base_dir
    run                    : with params set as desired, invokes exe and returns trajectory
    compute_reference      : does a very very very small dt run

    Example
    ======
    >>> from pf.pfasst import PFASST
    >>> pf = PFASST(params)
    >>> pf.p.nsteps = 32
    >>> pf.p.tfinal = 5.0
    >>> results = pf.run()
    """

    def __init__(self, params=None, **kwargs):
        if params is None:
            self.p = Params()
        else:
            self.p = params

        self.exe = self.p.exe
        for k, v in kwargs.items():
            setattr(self.p, k, v)

        try:
            mkdir(self.p.base_dir)
        except OSError:
            pass

    def _create_pf_string(self):
        param_list = self.p._make_list()
        return self.p.base_string.format(*param_list)

    def _write_to_file(self, string):
        """creates the input file on disk in the base_dir for the exe to run"""
        with open(self.p.base_dir + '/' + self.p.filename, 'w') as f:
            f.write(string)

    def _build_command(self):
        if self.p.nersc:
            command = ['srun', '-n', str(self.p.tasks)]
            command.extend([self.p.exe, self.p.base_dir + '/' + self.p.filename])
        else:
            command = ['mpirun', '-np', str(self.p.tasks), self.p.exe, \
                        self.p.base_dir + '/' + self.p.filename]

        return command

    def _pre_run_setup(self, ref):
        pf_string = self._create_pf_string()
        self._write_to_file(pf_string)
        pkl_path = self.p.get_pkl_path()

        if ref:
            pkl_path = pkl_path+'_ref'
        return pkl_path

    def _cleanup(self):
        for fname in glob.iglob(self.p.base_dir + '/*_soln'):
            remove(fname)

        for fname in glob.iglob('fort.*'):
            remove(fname)

        remove('final_solution')

    def run(self, ref=False):
        pkl_path = self._pre_run_setup(ref)

        try:
            trajectory = pd.read_pickle(pkl_path)
            total_times = {}
        except:
            try:
                if self.p.verbose:
                    nodes = ' '.join(map(str, self.p.nodes))

                    print('---- running pfasst: tasks={}, nodes={}, dt={} ----'.format(self.p.tasks, nodes, self.p.dt))

                command = self._build_command()

                output = check_output(command, stderr=STDOUT)
            except CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
            else:
                trajectory, total_times = self._get_trajectory_from_output(
                    output, self.p.nsteps, ref=ref)
                trajectory.to_pickle(pkl_path)
                self._cleanup()

        return trajectory, total_times

    def _get_trajectory_from_output(self, output, nsteps, ref=False):
        """If one was to replace reading stdout with reading dat files from disk
        this is the function that you would want to redefine.
        """
        trajectory = pd.DataFrame(columns=[
            'time', 'rank', 'step', 'iter', 'level', 'residual', 'solution', 'eigval'
        ])
        prog_state = re.compile(
            "resid: time: (.*) rank: (.*) step: (.*) iter: (.*) level: (.*) resid: (.*)"
        )
        time_state = re.compile(
            r'processor\ +(.*) returned from pfasst CPU time:\ +(.*) seconds')

        idx = 0
        total_times = {}
        if not ref:
            for line in tqdm(output.split('\n'), leave=False,
                             desc='parsing nsteps{}'.format(nsteps)):
                match = prog_state.search(line)
                if match:
                    time = float(match.group(1))
                    rank, step, iteration, level = map(int, match.group(2, 3, 4, 5))
                    residual_value = float(match.group(6))
                    if iteration < 0: continue

                    if self.p.solutions:
                        if time >= 1000.0:
                            path_to_solution = self.p.base_dir+'/'+\
                                            "time_{:09.5f}-rank_{:03d}-step_{:05d}-iter_{:03d}-level_{:01d}_soln".format(
                                                time, rank, step, iteration, level)
                        elif time >= 100.0:
                            path_to_solution = self.p.base_dir+'/'+\
                                            "time_{:08.5f}-rank_{:03d}-step_{:05d}-iter_{:03d}-level_{:01d}_soln".format(
                                                time, rank, step, iteration, level)
                        elif time >= 10.0:
                            path_to_solution = self.p.base_dir+'/'+\
                                            "time_{:07.5f}-rank_{:03d}-step_{:05d}-iter_{:03d}-level_{:01d}_soln".format(
                                                time, rank, step, iteration, level)
                        else:
                            path_to_solution = self.p.base_dir+'/'+\
                                            "time_{:06.5f}-rank_{:03d}-step_{:05d}-iter_{:03d}-level_{:01d}_soln".format(
                                                time, rank, step, iteration, level)
                        solution = self._get_solution(path_to_solution)
                        eigval = np.linalg.eigvals(solution)
                    else:
                        solution = None
                        eigval = None

                    trajectory.loc[idx] = time, rank, step, iteration, \
                                level, residual_value, solution, eigval
                    idx += 1

                time_match = time_state.search(line)
                if time_match:
                    rank = int(time_match.group(1))
                    cpu_time = float(time_match.group(2))
                    total_times[rank] = cpu_time


        if self.p.timings:
            timings = read_all_timings(dname='.')
            self._merge_timings_into(trajectory, timings)

        trajectory.sort_values(['step', 'level', 'iter'], inplace=True)
        trajectory.reset_index(inplace=True)
        if not self.p.solutions:
            trajectory.solution = trajectory.solution.astype('object')

        sol = self._get_solution('final_solution')
        last_row = len(trajectory) - 1
        if ref:
            last_row = 0
        trajectory.set_value(last_row, 'solution', sol)

        return trajectory, total_times

    def _merge_timings_into(self, trajectory, timings):
        timers_of_interest = ['exp', 'feval', 'omega']
        df_timing = pd.DataFrame(columns=Timing._fields)

        for i, time in enumerate(timings):
            if time.iter > 0:
                df_timing.loc[i] = time

        df_timing.drop(['start', 'end'], axis=1, inplace=True)
        pivoted = df_timing.pivot_table(
            index=['rank', 'step', 'iter'],
            columns='timer',
            values='delta',
            aggfunc=np.sum)

        for rank in trajectory['rank'].unique():
            for timer in timers_of_interest:
                trajectory.loc[(trajectory['level'] == 1), timer] = pivoted[
                    timer].values
        return trajectory

    @staticmethod
    def _get_solution(path_to_solution):
        print(path_to_solution)
        f = FortranFile(path_to_solution)
        solution = f.read_record(np.complex_)
        f.close()

        dim = int(np.sqrt(solution.shape))
        solution = np.reshape(solution, (dim, dim))

        return solution

    def compute_reference(self):

        params = self.p.pack()

        self.p.tasks = 1
        self.p.levels = 1
        self.p.nsteps = 2**11
        self.p.dt = self.p.tfinal / self.p.nsteps

        for k, v in self.p.ref_param_list.iteritems():
            setattr(self.p, k, v)

        traj, _ = self.run(ref=True)
        last_row = len(traj) - 1
        final_solution = traj.loc[last_row, 'solution']

        self.p.unpack(params)


        return final_solution, traj

    @staticmethod
    def back_transform(l, nparticles):
        nparts = l.shape[0]
        alpha = np.zeros(nparts, dtype=np.complex_)
        q = np.zeros(nparts, dtype=np.complex_)

        for i in range(nparts-1):
            alpha[i] = l[i, i+1]

        alpha[nparts-1] = l[0, nparts-1]

        q[0] = 0

        for j in range(nparts-1):
            q[j+1] = -2.0*np.log(2.0*alpha[j]) + q[j]

        q = q - q[np.int(nparticles/2)]
        return q

    def get_toda_solutions(self, traj):
        solutions = []
        for step in traj.step.unique():
            max_iter = traj[traj.step == step].iter.max()
            sol = traj.loc[((traj.step == step) &
                            (traj.iter == max_iter)), 'solution']
            solutions.append(sol.values[0])

        solutions = np.asarray(solutions)

        q_traj = np.zeros((self.p.nsteps, self.p.particles), dtype=np.complex_)
        p_traj = np.zeros((self.p.nsteps, self.p.particles), dtype=np.complex_)
        for j in range(self.p.nsteps):
            q_traj[j, :] = self.back_transform(solutions[j, :, :], nparticles=self.p.particles)
            p_traj[j, :] = 2.0*np.diag(solutions[j])

        return q_traj, p_traj

    def plot_toda(self, traj, maxparticles=None):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

        q_traj, p_traj = self.get_toda_solutions(traj)
        if maxparticles:
            particles = maxparticles
        else:
            particles = self.p.particles

        for i in range(particles):
            ax1.plot(np.linspace(0, self.p.tfinal, num=self.p.nsteps), q_traj[:, i])
            ax2.plot(np.linspace(0, self.p.tfinal, num=self.p.nsteps), p_traj[:, i])

        ax1.set_title('Position')
        ax2.set_title('Momentum')

        ax1.set_xlabel('Time')
        ax2.set_xlabel('Time')

        return fig, (ax1, ax2)


class Results(pd.DataFrame):
    """DataFrame derived container for holding all results. Implements
    auxiliary manipulations of df entries generally necessary for analyses
    related to PFASST calculations. Also implements plot functions that are
    derived from DataFrame.plots

    Parameters
    ==========
    params : holds all the PFASST class parameters as a Params object
    trajectory : run data from a PFASST.run() evaluation, default None
    pkl_path : path to a pickle file to recover data
    """

    def __init__(self, params, trajectory=None, pkl_path=None):
        columns = [
            'dt', 'nsteps', 'nodes', 'iterations', 'tfinal',
            'final_solution', 'total_times', 'trajectory'
        ]

        super(Results, self).__init__(columns=columns)
        self.p = params
        if pkl_path is not None:
            self.load(pkl_path)
        elif trajectory is not None:
            self._create_result_row(trajectory)

    def _create_result_row(self, idx, traj, total_times):
        final_solution = traj.loc[len(traj)-1, 'solution']
        iterations = 0

        for step in traj.step.unique():
            iterations += traj[traj['step'] == step].iter.max()

        iterations = iterations / len(traj.step.unique()+1)

        self.loc[idx] = self.p.dt, self.p.nsteps, self.p.nodes,  \
                      iterations, self.p.tfinal, \
                      final_solution, total_times, traj
        return

    def get_final_block(self, idx=0):
        """returns a the final block of the results class"""
        traj = self.loc[idx, 'trajectory']
        last_steps = self.loc[idx, 'nsteps'] - self.p.tasks

        final_block = traj[(traj['step'] > last_steps) & (traj['iter'] > 0)]
        return final_block

    def save(self, pkl_path):
        self.to_pickle(pkl_path)

    def load(self, pkl_path):
        try:
            with open(pkl_path, 'rb') as p:
                traj = pd.read_pickle(p)
        except:
            raise
        else:
            print('recovering results from pickle!')
            self.loc[0] = traj.loc[0]
            return

    def plot_convergence(self, x, y, **kwargs):
        """really nothing more than a wrapper to the df.plot function. maybe
        time to deprecate it.
        """
        self.plot(x, y, **kwargs)

    def plot_residual_vs_iteration_for_each_cpu(self,
                                                idx=0,
                                                trajectory=None,
                                                legend=True):
        """plots residual vs iteration for each cpu (nicely named function),
        optional input: a potentially longer trajectory. if no trajectory is
        supplied then it defaults to grabbing the final block of steps.
        """
        fig, ax = plt.subplots()
        ax.set_title('Residuals of last {} steps'.format(self.p.tasks))
        ax.set_xlabel('Iteration Number')
        ax.set_ylabel('Residual')
        # ax.set_ylim(bottom=1e-20)

        if idx is None and trajectory is None:
            trajectory = self.get_final_block(idx=0)
        elif idx:
            trajectory = self.get_final_block(idx=idx)
        else:
            trajectory = self.get_final_block(idx=0)

        for cpu in range(self.p.tasks):
            try:
                trajectory[(trajectory['rank'] == cpu) & (trajectory[
                    'level'] == self.p.levels)].plot(
                                                    'iter',
                                                    'residual',
                                                    logy=True,
                                                    label='cpu{}'.format(cpu),
                                                    marker='o',
                                                    ax=ax)
            except:
                pass

        if not legend:
            ax.legend_.remove()

        return fig, ax

    def plot_residual_vs_cpu_for_each_iteration(self,
                                                idx=None,
                                                trajectory=None,
                                                legend=True):
        """plots residual vs cpu for each iteration (nicely named function),
        optional input: a potentially longer trajectory. if no trajectory is
        supplied then it defaults to grabbing the final block of steps.
        """
        fig, ax = plt.subplots()
        ax.set_title('Residuals of last {} steps'.format(self.p.tasks))
        ax.set_xlabel('CPU Number')
        ax.set_ylabel('Residual')
        # ax.set_ylim(bottom=1e-20)

        if idx is None and trajectory is None:
            trajectory = self.get_final_block(idx=0)
        elif idx:
            trajectory = self.get_final_block(idx=idx)
        else:
            trajectory = self.get_final_block(idx=0)

        for iteration in range(self.p.iterations):
            try:
                trajectory[(trajectory['iter'] == iteration + 1) & (trajectory[
                    'level'] == self.p.levels)].plot(
                            'rank',
                            'residual',
                            logy=True,
                            label='iter{}'.format(iteration + 1),
                            marker='o',
                            ax=ax)
            except:
                pass

        if not legend:
            ax.legend_.remove()

        return fig, ax


class Experiment(object):
    """A variety of different pfasst-related experiments to be performed
    on a PFASST object
    """

    def __init__(self):
        pass

    def convergence_exp(self, pf, steps=[3, 4, 5, 6]):
        """convergence experiment for testing residual vs nsteps. has default
        number of steps, but can be overridden for larger tests. returns a Results
        class object where each row in the object corresponds to the parameters
        of the calculation and the trajectory from the calculation.
        """
        results = Results(pf.p)

        ref_soln, ref_traj = pf.compute_reference()

        errors = []
        nsteps = [2**i for i in steps]
        for i, step in tqdm(enumerate(nsteps), total=len(nsteps),
                            desc='nodes{}'.format(pf.p.nodes)):
            pf.p.nsteps = step
            pf.p.dt = pf.p.tfinal / pf.p.nsteps
            trajectory, total_times = pf.run()
            results._create_result_row(i, trajectory, total_times)

            errors.append(
                self.compute_error(results.loc[i]['final_solution'], ref_soln))

        results['error'] = errors
        results.astype({'dt': np.float_, 'nsteps': np.int_, 'error': np.float_})

        slope, _ = self.get_linear_fit_loglog(results.dt, results.error)
        print('slope = {}'.format(slope))

        return results

    @staticmethod
    def get_linear_fit_loglog(x, y):
        slope, intercept = np.polyfit(np.log10(x.values), np.log10(y.values), 1)

        return slope, intercept

    @staticmethod
    def compute_error(soln, ref_soln):
        error = np.linalg.norm(
            abs(soln - ref_soln))  # / np.linalg.norm(ref_soln)

        return error.max()


