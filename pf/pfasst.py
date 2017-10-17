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
import argparse
import glob
import re
from pprint import pprint
from os import remove, mkdir
from subprocess import check_output, STDOUT, CalledProcessError
from scipy.io import FortranFile
from pf.io import read_all_timings, Timing
import attr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


@attr.s(slots=True)
class Params(object):
    """Class containing all parameters necessary for running an NWC/PFASST
    calculation.

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
    filename = attr.ib(default=None)
    levels = attr.ib(default=1)
    tfinal = attr.ib(default=1.0)
    iterations = attr.ib(default=10)
    nsteps = attr.ib(default=16)
    nodes = attr.ib(default=[2], validator=attr.validators.instance_of(list))
    magnus = attr.ib(default=[1], validator=attr.validators.instance_of(list))
    sweeps = attr.ib(default=[3], validator=attr.validators.instance_of(list))
    exptol = attr.ib(default=['1.d-15'], validator=attr.validators.instance_of(list))
    nprob = attr.ib(default=4)
    tasks = attr.ib(default=1)
    basis = attr.ib(default='')
    molecule = attr.ib(default='')
    exact_dir = attr.ib(default='')
    base_dir = attr.ib(default='output', repr=False)
    verbose = attr.ib(default=False, repr=False)
    nersc = attr.ib(default=False)
    dt = attr.ib(default=None)
    timings = attr.ib(default=False)
    plist = attr.ib(repr=False)

    @plist.default
    def get_options_from_cli(self):
        if not self.nb:
            parser = argparse.ArgumentParser(
                description='From nml file for PFASST, generate exact/ref solutions'
            )
            parser.add_argument('--filename', type=str)
            parser.add_argument('--tfinal', type=float)
            parser.add_argument('--nsteps', type=int)
            parser.add_argument('--levels', type=int)
            parser.add_argument('--iterations', type=int)
            parser.add_argument('--sweeps', type=int, nargs='*')
            parser.add_argument('--nodes', type=int, nargs='*')
            parser.add_argument('--magnus', type=int, nargs='*')
            parser.add_argument(
                '--nprob', type=int, help='Default problem: toda')
            parser.add_argument('--tasks', type=int, help='Number of MPI tasks')
            parser.add_argument('--basis', type=str)
            parser.add_argument('--molecule', type=str)
            parser.add_argument('--exact_dir', type=str)
            parser.add_argument('--base_dir', type=str)
            parser.add_argument('-v', '--verbose', type=bool)
            args = parser.parse_args()
            self.overwrite_with(args)
            return args.__dict__

        return None

    def overwrite_with(self, args):
        for k, v in args.__dict__.iteritems():
            if v is not None:
                setattr(self, k, v)

    def __attrs_post_init__(self):
        if self.dt is None:
            self.dt = self.tfinal / self.nsteps
        else:
            self.nsteps = self.tfinal / self.dt

    def asdict(self):
        return attr.asdict(self)

    def print_params(self):
        params = {
            'tfinal': self.tfinal,
            'nsteps': self.nsteps,
            'dt': self.dt,
            'nodes': self.nodes,
            'magnus': self.magnus
        }
        pprint(params, width=1)


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
    home : a (str) path to the project's root e.g. '/home/bkrull/apps/pfasst/dev' is used to
    path to binary

    Exposed Methods
    ===============
    write_pfstring_to_file : writes the input file for executable in base_dir
    run                    : with params set as desired, invokes exe and returns trajectory
    compute_reference      : does a very very very small dt run

    Example
    ======
    >>> from pf.pfasst import PFASST
    >>> pf = PFASST('/home/bkrull/pfasst/')
    >>> pf.p.nsteps = 32
    >>> pf.p.tfinal = 5.0
    >>> results = pf.run()
    """

    def __init__(self, home, params=None, **kwargs):
        self.home = home
        if params is None:
            self.p = Params()
        else:
            self.p = params

        self.exe = self.home + 'main.exe'

        for k, v in kwargs.iteritems():
            settatr(self.p, k, v)

        self.base_string = "&PF_PARAMS\n\tnlevels = {}\n\tniters = {}\n\tqtype = 1\n\techo_timings = .false.\n\t\
abs_res_tol = 0.d-12\n\trel_res_tol = 0.d-12\n\tPipeline_G = .true.\n\tPFASST_pred = .true.\n/\n\n\
&PARAMS\n\tnnodes = {}\n\tnsweeps_pred = 1\n\tnsweeps = {}\n\t\
magnus_order = {}\n\tTfin = {}\n\tnsteps = {}\n\texptol = {}\
\n\tnprob = {}\n\tbasis = {}\n\tmolecule = {}\n\texact_dir = {}\n/\n"

        if self.p.filename:
            with open(self.p.base_dir + '/' + self.p.filename, 'r') as f:
                content = f.read()
            self._override_default_params(content)

        if self.p.nprob == 1:
            self.p.filename = 'rabi.nml'
        elif self.p.nprob == 2:
            self.p.filename = 'tdrabi.nml'
        elif self.p.nprob == 3:
            mol_name = re.findall(r'([A-Z]+)', self.p.molecule, re.IGNORECASE)
            self.p.filename = ''.join(mol_name) + '_' + self.p.basis + '.nml'
        else:
            self.p.filename = 'toda.nml'

        try:
            mkdir(self.p.base_dir)
        except OSError:
            pass

        self.pkl = self.p.base_dir + '/nprob_{}-tfinal_{}-dt_{}-levels_{}'+ \
                   '-coarsenodes_{}-coarsemagnus_{}-tasks_{}.pkl'

    def _create_pf_string(self):
        self.p.magnus = self._make_sure_is_list(self.p.magnus)
        self.p.nodes = self._make_sure_is_list(self.p.nodes)
        self.p.sweeps = self._make_sure_is_list(self.p.sweeps)
        self.p.exptol = self._make_sure_is_list(self.p.exptol)
        nodes = ' '.join(map(str, self.p.nodes))
        magnus = ' '.join(map(str, self.p.magnus))
        sweeps = ' '.join(map(str, self.p.sweeps))
        exptol = ' '.join(self.p.exptol)

        self.pfstring = self.base_string.format(self.p.levels, self.p.iterations,\
                                                nodes, sweeps, magnus, self.p.tfinal, \
                                                self.p.nsteps, exptol, self.p.nprob, \
                                                "\'"+self.p.basis+"\'", \
                                                "\'"+self.p.molecule+"\'", \
                                                "\'"+self.p.exact_dir+"\'")

        return self.pfstring

    def _override_default_params(self, content):
        try:
            match = re.search(r'magnus_order = (.+)', content)
        except AttributeError:
            pass
        else:
            self.p.magnus = map(int, match.group(1).split())

        try:
            match = re.search(r'nsweeps = (.+)', content)
        except AttributeError:
            pass
        else:
            self.p.sweeps = map(int, match.group(1).split())

        try:
            match = re.search(r'nnodes = (.+)', content)
        except AttributeError:
            pass
        else:
            self.p.nodes = map(int, match.group(1).split())

        try:
            self.p.levels = int(
                re.search(r'nlevels\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.p.iterations = int(
                re.search(r'niters\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.p.tfinal = float(re.search(r'Tfin\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.p.nsteps = int(re.search(r'nsteps\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.p.dt = self.p.tfinal / self.p.nsteps
        except AttributeError:
            pass
        try:
            self.p.nprob = int(
                re.search(r'nprob\ =\ ([0-9])', content).group(1))
        except AttributeError:
            pass

        if self.p.nprob == 3:
            try:
                self.p.basis = re.search(r'basis\ =\ \'(.+)\'',
                                         content).group(1)
            except AttributeError:
                pass
            try:
                self.p.molecule = re.search(r'molecule\ =\ \'(.+)\'',
                                            content).group(1)
            except AttributeError:
                pass
            try:
                self.p.exact_dir = re.search(r'exact\_dir\ =\ \'(.+)\'',
                                             content).group(1)
            except AttributeError:
                pass
        else:
            self.p.basis = ''
            self.p.molecule = ''
            self.p.exact_dir = ''

    def write_to_file(self):
        """creates the input file on disk in the base_dir for the exe to run"""
        self.pfstring = self._create_pf_string()

        with open(self.p.base_dir + '/' + self.p.filename, 'w') as f:
            f.write(self.pfstring)

    def _build_command(self):
        if self.p.nersc:
            command = ['srun', '-n', str(self.p.tasks)]
            for option in self.p.nersc:
                command.extend(option)
            command.extend([self.exe, self.p.base_dir + '/' + self.p.filename])
        else:
            command = ['mpirun', '-np', str(self.p.tasks), self.exe, \
                        self.p.base_dir + '/' + self.p.filename]

        return command

    def _pre_run_setup(self):
        self.p.magnus = self._make_sure_is_list(self.p.magnus)
        self.p.nodes = self._make_sure_is_list(self.p.nodes)
        self.p.sweeps = self._make_sure_is_list(self.p.sweeps)
        self._create_pf_string()
        self.write_to_file()
        pkl_path = self.pkl.format(self.p.nprob, self.p.tfinal, self.p.dt,
                                   self.p.levels, self.p.nodes[0],
                                   self.p.magnus[0], self.p.tasks)

        return pkl_path

    def _cleanup(self):
        for file in glob.iglob(self.p.base_dir + '/*_soln'):
            remove(file)

        for file in glob.iglob('fort.*'):
            remove(file)

    def run(self, ref=False):
        pkl_path = self._pre_run_setup()

        try:
            trajectory = pd.read_pickle(pkl_path)
        except:
            try:
                if self.p.verbose:
                    nodes = ' '.join(map(str, self.p.nodes))
                    magnus = ' '.join(map(str, self.p.magnus))

                    print '---- running pfasst: tasks={}, nodes={}, magnus={}, dt={} ----'.format(
                        self.p.tasks, nodes, magnus, self.p.dt)

                command = self._build_command()
                output = check_output(command, stderr=STDOUT)
            except CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
            else:
                trajectory = self._get_trajectory_from_output(output, ref=ref)
                trajectory.to_pickle(pkl_path)
                self._cleanup()

        return trajectory

    def _get_trajectory_from_output(self, output, ref=False):
        columns = [
            'time', 'rank', 'step', 'iter', 'level', 'residual', 'solution']
        trajectory = pd.DataFrame(columns=columns)
        prog_state = re.compile(
            "resid: time: (.*) rank: (.*) step: (.*) iter: (.*) level: (.*) resid: (.*)"
        )

        idx = 0
        for line in output.split('\n'):
            match = prog_state.search(line)
            if match:
                time = float(match.group(1))
                rank, step, iteration, level = map(int, match.group(2, 3, 4, 5))
                residual_value = float(match.group(6))

                path_to_solution = self.p.base_dir+'/'+\
                                   "time_{:05.3f}-rank_{:03d}-step_{:04d}-iter_{:02d}-level_{:01d}_soln".format(
                                       time, rank, step, iteration, level)
                solution = self._get_solution(path_to_solution)
                trajectory.loc[idx] = time, rank, step, iteration, \
                              level, residual_value, solution
                idx += 1


        if not ref and self.p.timings:
            timings = read_all_timings(dname='.')
            self._merge_timings_into(trajectory, timings)

        trajectory.sort_values(['step', 'level', 'iter'], inplace=True)

        return trajectory

    def _merge_timings_into(self, trajectory, timings):
        timers_of_interest = ['exp', 'feval', 'omega']
        df_timing = pd.DataFrame(columns=Timing._fields)

        for i, time in enumerate(timings):
            if time.iter > 0:
                df_timing.loc[i] = time

        df_timing.drop(['start', 'end'], axis=1, inplace=True)
        pivoted = df_timing.pivot_table(index=['rank', 'step', 'iter'],
                                        columns='timer', values='delta', aggfunc=np.sum)

        for rank in trajectory['rank'].unique():
            for timer in timers_of_interest:
                    trajectory.loc[((trajectory['rank'] == rank-1) &
                                    (trajectory['level'] == 1)), timer] = pivoted[timer].values
        return trajectory

    @staticmethod
    def _get_solution(path_to_solution):
        f = FortranFile(path_to_solution)
        solution = f.read_record(np.complex_)
        f.close()

        dim = int(np.sqrt(solution.shape))
        solution = np.reshape(solution, (dim, dim))

        return solution

    def compute_reference(self):
        self.p.tasks = 1
        self.p.nsteps = 2**8
        self.p.nodes = 3
        self.p.magnus = 2
        self.p.dt = self.p.tfinal / self.p.nsteps

        traj = self.run(ref=True)

        return traj

    @staticmethod
    def _make_sure_is_list(thing):
        if type(thing) is not list:
            return [thing]
        else:
            return thing


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
            'dt', 'nsteps', 'nodes', 'magnus', 'iterations', 'tfinal',
            'final_solution', 'trajectory'
        ]

        super(Results, self).__init__(columns=columns)
        self.p = params
        if pkl_path is not None:
            self.load(pkl_path)
        elif trajectory is not None:
            self._create_result_row(trajectory)

    def _create_result_row(self, idx, trajectory):
        final_solution = self.get_final_solution(trajectory)

        self.loc[idx] = self.p.dt, self.p.nsteps, \
                      self.p.nodes, self.p.magnus, \
                      self.p.iterations, self.p.tfinal, \
                      final_solution, trajectory
        return

    def get_final_solution(self, trajectory):
        """gets the final solution for the finest level from a trajectory"""
        solution = trajectory[(trajectory['rank'] == self.p.tasks - 1) & \
                        (trajectory['iter'] == self.p.iterations) & \
                        (trajectory['time'] == self.p.tfinal) & \
                        (trajectory['level'] == self.p.levels) & \
                        (trajectory['step'] == self.p.nsteps)]['solution'].values

        return solution

    def get_final_block(self, idx=0):
        """returns a the final block of the results class"""
        traj = self.loc[idx]['trajectory']
        last_steps = self.loc[idx]['nsteps'] - self.p.tasks

        final_block = traj[(traj['step'] > last_steps)]
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
            print 'recovering results from pickle!'
            self.loc[0] = traj.loc[0]
            return

    def plot_convergence(self, x, y, **kwargs):
        """really nothing more than a wrapper to the df.plot function. maybe
        time to deprecate it.
        """
        self.plot(x, y, **kwargs)

    def plot_residual_vs_iteration_for_each_cpu(self, trajectory=None, legend=True):
        """plots residual vs iteration for each cpu (nicely named function),
        optional input: a potentially longer trajectory. if no trajectory is
        supplied then it defaults to grabbing the final block of steps.
        """
        fig, ax = plt.subplots()
        ax.set_title('Residuals of last {} steps'.format(self.p.tasks))
        ax.set_xlabel('Iteration Number')
        ax.set_ylabel('Log$_{10}$ Residual)')

        if trajectory is None:
            trajectory = self.get_final_block()

        for cpu in range(self.p.tasks):
            trajectory[(trajectory['rank'] == cpu) & (trajectory[
                'level'] == 1)].plot(
                    'iter',
                    'residual',
                    logy=True,
                    label='cpu{}'.format(cpu),
                    marker='o',
                    ax=ax)
        if legend:
            ax.legend()

    def plot_residual_vs_cpu_for_each_iteration(self, trajectory=None, legend=True):
        """plots residual vs cpu for each iteration (nicely named function),
        optional input: a potentially longer trajectory. if no trajectory is
        supplied then it defaults to grabbing the final block of steps.
        """
        fig, ax = plt.subplots()
        ax.set_title('Residuals of last {} steps'.format(self.p.tasks))
        ax.set_xlabel('CPU Number')
        ax.set_ylabel('Log$_{10}$ Residual)')

        if trajectory is None:
            trajectory = self.get_final_block()

        for iteration in range(self.p.iterations):
            trajectory[(trajectory['iter'] == iteration + 1) & (trajectory[
                'level'] == 1)].sort_values('rank').plot(
                    'rank',
                    'residual',
                    logy=True,
                    label='iter{}'.format(iteration + 1),
                    marker='o',
                    ax=ax)
        if legend:
            ax.legend()


class Experiment(object):
    """A variety of different pfasst-related experiments to be performed
    on a PFASST object
    """

    def __init__(self):
        pass

    def nodes_exp(self, pf, nodes=[2, 3], magnus=[1, 2]):
        """an experiment for testing the convergence of different number of sdc
        nodes. wraps around the convergence experiment.
        """
        results = []

        for i, node in enumerate(nodes):
            pf.nodes = node
            pf.magnus = magnus[i]
            r = self.convergence_exp(pf, steps=[4, 7])
            results.append(r)

        return results

    def convergence_exp(self, pf, steps=[3, 4, 5, 6]):
        """convergence experiment for testing residual vs nsteps. has default
        number of steps, but can be overridden for larger tests. returns a Results
        class object where each row in the object corresponds to the parameters
        of the calculation and the trajectory from the calculation.
        """
        results = Results(pf.p)

        nsteps = pf.p.nsteps
        nodes = pf.p.nodes
        magnus = pf.p.magnus
        tasks = pf.p.tasks

        ref_traj = pf.compute_reference()
        ref_soln = results.get_final_solution(ref_traj)

        pf.p.nsteps = nsteps
        pf.p.nodes = nodes
        pf.p.magnus = magnus
        pf.p.tasks = tasks

        errors = []
        nsteps = [2**i for i in steps]
        for i, step in enumerate(nsteps):
            pf.p.nsteps = step
            pf.p.dt = pf.p.tfinal / pf.p.nsteps
            trajectory = pf.run()
            results._create_result_row(i, trajectory)

            errors.append(
                self.compute_error(results.loc[i]['final_solution'], ref_soln))

        results['error'] = errors
        results.astype({'dt': np.float_, 'nsteps': np.int_, 'error': np.float_})

        slope, _ = self.get_linear_fit_loglog(results.dt, results.error)
        print 'slope = {}'.format(slope)

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
