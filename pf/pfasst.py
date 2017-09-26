import argparse
import glob
import re
from collections import namedtuple
from pprint import pprint
from os import remove
from subprocess import check_output
from scipy.io import FortranFile
import pandas as pd
import numpy as np

State = namedtuple(
    'State', ['time', 'rank', 'step', 'iter', 'level', 'residual', 'solution'])
Timing = namedtuple('Timing', ['task', 'step', 'iter', 'level', 'time'])

def create_pf_params():
    parser = argparse.ArgumentParser(
        description='From nml file for PFASST, generate exact/ref solutions')
    parser.add_argument('--pfasst_input', type=str, default=None)
    parser.add_argument('--tfinal', type=float, default=1.0)
    parser.add_argument('--dt', type=float, default=0.1)
    parser.add_argument('--nsteps', type=float, default=16)
    parser.add_argument('--levels', type=int, default=1)
    parser.add_argument('--iters', type=int, default=5)
    parser.add_argument('--sweeps', type=int, nargs='*', default=[3])
    parser.add_argument('--nodes', type=int, nargs='*', default=[2])
    parser.add_argument('--magnus', type=int, nargs='*', default=[1])
    parser.add_argument(
        '--nprob', type=int, default=4, help='Default problem: toda')
    parser.add_argument(
        '--tasks', type=int, default=1, help='Number of MPI tasks')
    parser.add_argument('--basis', type=str, default='')
    parser.add_argument('--molecule', type=str, default='')
    parser.add_argument('--exact_dir', type=str, default='')
    parser.add_argument('--base_dir', type=str, default='output')
    parser.add_argument('-v', '--verbose', type=bool, default=False)

    args = parser.parse_args()
    return args


def create_nb_pf_params():
    params = {
        'pfasst_input': None,
        'levels': 1,
        'sweeps': [3],
        'tfinal': 1.0,
        'iters': 5,
        'dt': 0.1,
        'nsteps': 16,
        'nodes': [2],
        'magnus': [1],
        'nprob': 4,
        'tasks': 1,
        'basis': '',
        'molecule': '',
        'exact_dir': '',
        'base_dir': 'output',
        'verbose': False
    }
    return params


class PFASST:

    def __init__(self, home, **kwargs):
        self.params = kwargs
        self.home = home
        self.exe = self.home + 'main.exe'
        self.levels = self.params['levels']
        self.iters = self.params['iters']
        self.verbose = self.params['verbose']
        self.input_file = self.params['pfasst_input']
        self.tasks = self.params['tasks']
        self.tfinal = self.params['tfinal']
        self.sweeps = self.params['sweeps']
        self.base_dir = self.params['base_dir']
        self.output = ''
        if type(self.sweeps) is int:
            self.sweeps = [self.sweeps]
        self.magnus = self.params['magnus']
        if type(self.magnus) is int:
            self.magnus = [self.magnus]
        self.nodes = self.params['nodes']
        if type(self.nodes) is int:
            self.nodes = [self.nodes]
        if self.params['nsteps']:
            self.nsteps = self.params['nsteps']
            self.dt = self.tfinal / self.nsteps
        elif self.params['dt']:
            self.dt = self.params['dt']
            self.nsteps = self.tfinal * self.dt

        self.nprob = self.params['nprob']
        self.basis = self.params['basis']
        self.molecule = self.params['molecule']
        self.exact_dir = self.params['exact_dir']

        self.base_string = "&PF_PARAMS\n\tnlevels = {}\n\tniters = {}\n\tqtype = 1\n\t\
abs_res_tol = 0.d-12\n\trel_res_tol = 0.d-12\n\tPipeline_G = .true.\n\tPFASST_pred = .true.\n/\n\
&PARAMS\n\tnnodes = {}\n\tnsweeps_pred = 1\n\tnsweeps = {}\n\t\
magnus_order = {}\n\tTfin = {}\n\tnsteps = {}\
\n\tnprob = {}\n\tbasis = {}\n\tmolecule = {}\n\texact_dir = {}\n/\n"

        if self.input_file:
            with open(self.input_file, 'r') as f:
                content = f.read()
            self._override_default_params(content)

        if self.nprob == 1:
            self.filename = 'rabi.nml'
        elif self.nprob == 2:
            self.filename = 'tdrabi.nml'
        elif self.nprob == 3:
            mol_name = re.findall(r'([A-Z]+)', self.molecule, re.IGNORECASE)
            self.filename = ''.join(mol_name) + '_' + self.basis + '.nml'
        else:
            self.filename = 'toda.nml'

        if not self.base_dir:
            self.base_dir = 'output'

        self.pkl = self.base_dir + '/nprob_{}-tfinal_{}-dt_{}-levels_{}-coarsenodes_{}-coarsemagnus_{}-tasks_{}.pkl'
        self.params = self._save_params()

    def _create_pf_string(self):
        if type(self.magnus) is int:
            self.magnus = [self.magnus]
        if type(self.nodes) is int:
            self.nodes = [self.nodes]
        if type(self.sweeps) is int:
            self.sweeps = [self.sweeps]
        nodes = ' '.join(map(str, self.nodes))
        magnus = ' '.join(map(str, self.magnus))
        sweeps = ' '.join(map(str, self.sweeps))

        pfinp = self.base_string.format(self.levels, self.iters,\
                                        nodes, sweeps, magnus, self.tfinal, \
                                        self.nsteps, self.nprob, \
                                        "\'"+self.basis+"\'", \
                                        "\'"+self.molecule+"\'", \
                                        "\'"+self.exact_dir+"\'")

        return pfinp

    def _override_default_params(self, content):
        try:
            match = re.search(r'magnus_order = (.+)', content)
        except AttributeError:
            pass
        else:
            self.magnus = map(int, match.group(1).split())

        try:
            match = re.search(r'nsweeps = (.+)', content)
        except AttributeError:
            pass
        else:
            self.sweeps = map(int, match.group(1).split())

        try:
            match = re.search(r'nnodes = (.+)', content)
        except AttributeError:
            pass
        else:
            self.nodes = map(int, match.group(1).split())

        try:
            self.levels = int(re.search(r'nlevels\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.iters = int(re.search(r'niters\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.tfinal = float(re.search(r'Tfin\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.nsteps = int(re.search(r'nsteps\ =\ (.+)', content).group(1))
        except AttributeError:
            pass
        try:
            self.dt = self.tfinal / self.nsteps
        except AttributeError:
            pass
        try:
            self.nprob = int(re.search(r'nprob\ =\ ([0-9])', content).group(1))
        except AttributeError:
            pass

        if self.nprob == 3:
            try:
                self.basis = re.search(r'basis\ =\ \'(.+)\'', content).group(1)
            except AttributeError:
                pass
            try:
                self.molecule = re.search(r'molecule\ =\ \'(.+)\'',
                                          content).group(1)
            except AttributeError:
                pass
            try:
                self.exact_dir = re.search(r'exact\_dir\ =\ \'(.+)\'',
                                           content).group(1)
            except AttributeError:
                pass
        else:
            self.basis = ''
            self.molecule = ''
            self.exact_dir = ''

    def _save_params(self):
        params = {
            'tasks': self.tasks,
            'pfasst_input': self.input_file,
            'tfinal': self.tfinal,
            'iters': self.iters,
            'sweeps': self.sweeps,
            'magnus': self.magnus,
            'nodes': self.nodes,
            'nprob': self.nprob,
            'basis': self.basis,
            'molecule': self.molecule,
            'exact_dir': self.exact_dir,
            'base_dir': self.base_dir,
            'nsteps': self.nsteps,
            'dt': self.dt,
            'verbose': self.verbose,
            'pkl': self.pkl
        }

        self.params = params
        return params

    def write_to_file(self):
        self.pfinp = self._create_pf_string()

        with open(self.base_dir + '/' + self.filename, 'w') as f:
            f.write(self.pfinp)

    def print_params(self):
        params = {
            'tfinal': self.tfinal,
            'nsteps': self.nsteps,
            'dt': self.dt,
            'nodes': self.nodes,
            'magnus': self.magnus
        }
        pprint(params, width=1)

    def run(self):
        if type(self.nodes) is int:
            self.nodes = [self.nodes]
        if type(self.magnus) is int:
            self.magnus = [self.magnus]
        if type(self.sweeps) is int:
            self.sweeps = [self.sweeps]
        self._create_pf_string()
        self.write_to_file()
        pkl = self.pkl.format(self.nprob, self.tfinal, self.dt, self.levels,
                              self.nodes[0], self.magnus[0], self.tasks)

        try:
            trajectory = self._get_traj_from_pkl(pkl)
        except:
            if self.verbose:
                nodes = ' '.join(map(str, nodes))
                magnus = ' '.join(map(str, magnus))

                print '---- running pfasst: tasks={}, nodes={}, magnus={}, dt={} ----'.format(
                    self.tasks, self.nodes, self.magnus, self.dt)

            output = check_output([ 'mpirun', '-np', str(self.tasks), self.exe, \
                self.base_dir + '/' + self.filename])

            trajectory = self._get_traj_from_output(output)

            self._save_pickle(trajectory, pkl)
            self._cleanup()

        return trajectory

    @staticmethod
    def _save_pickle(traj, pkl_file):
        with open(pkl_file, 'wb') as pkl:
            traj.to_pickle(pkl)

    def _cleanup(self):
        for file in glob.iglob(self.base_dir + '/*_soln'):
            remove(file)

    def _get_traj_from_pkl(self, pkl):
        try:
            with open(pkl, 'rb') as p:
                traj = pd.read_pickle(p)
        except:
            raise
        else:
            print 'recovering results from pickle!'
            return traj

    def _get_traj_from_output(self, output):
        trajectory = []
        prog = re.compile(
            "resid: time: (.*) rank: (.*) step: (.*) iter: (.*) level: (.*) resid: (.*)"
        )

        for line in output.split('\n'):
            match = prog.search(line)
            if match:
                time = float(match.group(1))
                rank, step, iteration, level = map(int, match.group(2, 3, 4, 5))
                residual_value = float(match.group(6))

                path_to_solution = self.base_dir+'/'+\
                                   "time_{:05.3f}-rank_{:03d}-step_{:04d}-iter_{:02d}-level_{:01d}_soln".format(
                                       time, rank, step, iteration, level)
                solution = self._get_solution(path_to_solution)
                trajectory.append(
                    State(time, rank, step, iteration, level, residual_value,
                          solution))

        df = pd.DataFrame(trajectory, columns=State._fields)
        return df

    def _get_solution(self, path_to_solution):
        f = FortranFile(path_to_solution)
        solution = f.read_record(np.complex_)
        f.close()

        dim = int(np.sqrt(solution.shape))
        solution = np.reshape(solution, (dim, dim))

        return solution

    def compute_reference(self):
        self.nsteps = 2**10
        self.magnus = 2

        self.dt = self.tfinal / self.nsteps

        traj = self.run()

        return traj

    def get_final_solution(self, trajectory):
        solution = trajectory[(trajectory['rank'] == self.tasks - 1) & \
                            (trajectory['iter'] == self.iters) & \
                            (trajectory['time'] == self.tfinal) & \
                            (trajectory['level'] == self.levels) & \
                            (trajectory['step'] == self.nsteps)]['solution'].values

        return solution


class Experiment(pd.DataFrame):
    """A variety of different pfasst-related experiments to be performed
    on a PFASST object
    """

    def __init__(self, columns=['dt', 'nsteps', 'final_solution', 'error', 'trajectory']):
        pd.DataFrame.__init__(self, columns=columns)

    def nodes_exp(self, pf):
        nodes = [2, 3]
        magnus = [1, 2]
        solns = {}
        dts = []

        for i, node in enumerate(nodes):
            pf.nodes = node
            pf.magnus = magnus[i]
            slope = self.short_convergence_exp(pf)
            print slope

        return True

    def short_convergence_exp(self, pf, steps=[4, 7]):
        experiment = Experiment()

        ref_traj = pf.compute_reference()
        ref_solution = pf.get_final_solution(ref_traj)

        nsteps = [2**i for i in steps]
        for i, step in enumerate(nsteps):
            pf.nsteps = step
            pf.dt = pf.tfinal / pf.nsteps
            trajectory = pf.run()
            final_solution = pf.get_final_solution(trajectory)
            error = self.compute_error(final_solution, ref_solution)

            experiment.loc[i] = pf.dt, pf.nsteps, final_solution, error, trajectory

        experiment.astype({'dt': np.float_,
                           'nsteps': np.int_,
                           'error': np.float_})

        return experiment

    def convergence_exp(self, pf, steps=[4, 5, 6, 7]):
        experiment = Experiment()

        ref_traj = pf.compute_reference()
        ref_solution = pf.get_final_solution(ref_traj)

        nsteps = [2**i for i in steps]
        for i, step in enumerate(nsteps):
            pf.nsteps = step
            pf.dt = pf.tfinal / pf.nsteps
            trajectory = pf.run()
            final_solution = pf.get_final_solution(trajectory)
            error = self.compute_error(final_solution, ref_solution)

            experiment.loc[i] = pf.dt, pf.nsteps, final_solution, error, trajectory

        experiment.astype({'dt': np.float_,
                           'nsteps': np.int_,
                           'error': np.float_})
        return experiment

    @staticmethod
    def get_linear_fit_loglog(x, y):
        slope, intercept = np.polyfit(np.log10(x.values),
                                      np.log10(y.values), 1)

        return slope, intercept

    @staticmethod
    def compute_error(soln, ref_soln):
        error = np.linalg.norm(abs(soln - ref_soln))  # / np.linalg.norm(ref_soln)

        return error.max()

    def plot_convergence(self, x, y, **kwargs):
        self.plot(x, y, **kwargs)
