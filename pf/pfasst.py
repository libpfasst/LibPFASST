import argparse
import glob
import re
from collections import namedtuple
from pprint import pprint
from os import remove, getcwd, chdir, makedirs, devnull
from subprocess import call, check_output
from scipy.io import FortranFile
import numpy as np
import pandas as pd

State = namedtuple(
    'State', ['time', 'rank', 'step', 'iter', 'level', 'residual', 'solution'])
Timing = namedtuple('Timing', ['task', 'step', 'iter', 'level', 'time'])


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
        if type(self.sweeps) is int:
            self.sweeps = [self.sweeps]
        self.magnus = self.params['magnus']
        if type(self.magnus) is int:
            self.magnus = [self.magnus]
        self.nodes = self.params['nodes']
        if type(self.nodes) is int:
            self.nodes = [self.nodes]
        self.nprob = self.params['nprob']
        self.basis = self.params['basis']
        self.molecule = self.params['molecule']
        self.exact_dir = self.params['exact_dir']
        self.base_dir = self.params['base_dir']
        if self.params['nsteps']:
            self.nsteps = self.params['nsteps']
            self.dt = self.tfinal / self.nsteps
        elif self.params['dt']:
            self.dt = self.params['dt']
            self.nsteps = self.tfinal * self.dt
        self.residuals = []
        self.output = ''

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

        return params

    def write_to_file(self):
        self.pfinp = self._create_pf_string()

        with open(self.base_dir + '/' + self.filename, 'w') as f:
            f.write(self.pfinp)

    def print_params(self):
        params = {}
        params['tfinal'] = self.tfinal
        params['nsteps'] = self.nsteps
        params['dt'] = self.dt
        params['nodes'] = self.nodes
        params['magnus'] = self.magnus

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

            output = check_output([
                'srun', '-n', str(self.tasks), self.exe,
                self.base_dir + '/' + self.filename
            ])

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


class NWC:

    def __init__(self, home, params):
        print 'Setting up NWC'
        self.home = home
        self.params = params
        self.nwc = self.home + '/nwc.sh'
        self.exact_dir = params['exact_dir']
        self.cwd = getcwd() + '/'
        self.base_string = 'echo\nscratch_dir ./scratch\npermanent_dir ./perm\n\n'\
                'start calc\n\ngeometry "system" units angstroms nocenter noautoz noautosym\n {}\nend\n\n' \
                'set geometry "system" \n\nbasis\n * library {}\nend\n\n' \
                'dft\n xc HFexch\n end\n\ntask dft energy\n\n' \
                'rt_tddft\n  checklvl 2\n  tmax {}\n  dt {}\n  ' \
                'field "kick"\ntype delta\npolarization x\nmax 75.0\nend\n\n' \
                'excite "system" with "kick"\nend\ntask dft rt_tddft\n\n'

        try:
            just_made_input = self._create_nwc_inp()
        except:
            raise Exception('unable to create input file')

        if just_made_input:
            try:
                success = self._run_nwc()
            except:
                raise Exception('unsuccessful NWC run, check input!')
        else:
            print 'exact dirs already exist'

        self.xtrans, self.ytrans = self._parse_transforms()
        self.initial = self._get_nwc_dmat('initial_condition')
        self.final = self._get_nwc_dmat('final_solution')

        self._write_zndarray_to_fortranfile(
            self.exact_dir + '/x_transform.ints', self.xtrans)
        self._write_zndarray_to_fortranfile(
            self.exact_dir + '/y_transform.ints', self.ytrans)
        self._write_zndarray_to_fortranfile(
            self.exact_dir + '/initial_condition.b', self.initial)
        self._write_mo_and_ao_sizes(self.exact_dir, self.xtrans.shape)

    def _create_nwc_inp(self):
        geometry = '\n'.join(self.params['molecule'].split(';'))
        nwc_inp = self.base_string.format(geometry, self.params['basis'],
                                          self.params['tfinal'],
                                          self.params['tfinal'] / 2048.)

        try:
            makedirs(self.exact_dir)
            makedirs(self.exact_dir + '/perm')
            makedirs(self.exact_dir + '/scratch')
        except OSError:
            return False
        else:
            with open(self.exact_dir + '/nw.inp', 'w') as f:
                f.write(nwc_inp)

        return True

    def _run_nwc(self):
        FNULL = open(devnull, 'w')
        chdir(self.exact_dir)
        print '---- running nwc ----'
        call([self.nwc, 'nw'], stdout=FNULL)
        chdir(self.cwd)
        FNULL.close()

        return True

    def _parse_transforms(self):
        with open(self.exact_dir + 'nw.out', 'r') as f:
            content = f.readlines()

        dims = {}
        vals = {}
        count = 0
        started = False
        for i in content:
            if 'transform start' in i:
                started = True
                continue

            if started:
                if r'transform matrix' in i:
                    l = []
                    count += 1
                    key = re.search(r'\%\s.+\s\%', i).group(0)[2:-2]
                elif 'global' in i:
                    dim = re.findall(r'[0-9]+\:[0-9]+', i)
                    for j in dim:
                        indices = re.findall('[0-9]+', j)
                        low = int(indices[0])
                        high = int(indices[1])
                        dims[key] = high - low + 1
                elif r'%%%%' in i:
                    vals[key] = l
                    l = []
                    if count == 2:
                        break
                elif i == '\n':
                    pass
                else:
                    l.append(i)

        for k, v in vals.iteritems():
            vals[k] = self._get_complex_matrix_vals(dims[k], v)

        return vals.values()

    def _get_complex_matrix_vals(self, dim, raw_matrix_text):
        """get complex values from matrix"""
        matrix = np.zeros((dim, dim), dtype=np.complex128)
        real = np.zeros((dim, dim), dtype=np.float64)
        imag = np.zeros((dim, dim), dtype=np.float64)
        for i in raw_matrix_text:
            if r'.' not in i and r'---' not in i:  #column label
                cols = [int(j) for j in re.findall('[0-9]+', i)]
            elif r'---' not in i:
                row = int(re.search('[0-9]+', i).group(0))
                elems = [j.split(', ') for j in \
                        re.findall(r'\-*[0-9]*\.[0-9]+\,\s*\-*[0-9]*\.[0-9]*', i)]

                for j, e in enumerate(elems):
                    real[row - 1, cols[j] - 1] = float(e[0])
                    imag[row - 1, cols[j] - 1] = float(e[1])
            else:
                pass

        matrix.real = real
        matrix.imag = imag

        return matrix

    def _get_nwc_dmat(self, fname):
        with open(self.exact_dir + fname, 'r') as f:
            content = f.readlines()

        l = []
        dims = []
        for i in content:
            if 'global' in i:
                dim = re.findall(r'[0-9]+\:[0-9]+', i)
                for j in dim:
                    indices = re.findall('[0-9]+', j)
                    low = int(indices[0])
                    high = int(indices[1])
                    dims.append(high - low + 1)
            elif i == '\n' or '=========' in i or r'<rt_tddft>' in i:
                pass
            else:
                l.append(i)

        dmat = self._get_complex_matrix_vals(dims[0], l)

        return dmat

    def _write_mo_and_ao_sizes(self, dirname, mo_ao_shape):
        with open(dirname + '/size_mo', 'w') as f:
            f.write(str(mo_ao_shape[0]))

        with open(dirname + '/size_ao', 'w') as f:
            f.write(str(mo_ao_shape[1]))

    def _write_zndarray_to_fortranfile(self, fname, zndarray):
        f = FortranFile(fname, 'w')
        f.write_record(zndarray)
        f.close()


class Experiments:
    """A variety of different pfasst-related experiments to be performed
    on a PFASST object
    """

    def __init__(self, plot=False):
        self.plot = plot

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
        nsteps = [2**i for i in steps]
        solns = {}
        trajs = {}
        dts = []

        for step in nsteps:
            pf.nsteps = step
            pf.dt = pf.tfinal / pf.nsteps
            dts.append(pf.dt)

            trajectory = pf.run()
            trajs[pf.dt] = trajectory
            solns[pf.dt] = get_final_solution(pf, trajectory)

        ref_traj = pf.compute_reference()
        ref_solution = get_final_solution(pf, ref_traj)

        errors = [compute_error(solns[dt], ref_solution) for dt in dts]

        return solns, dts, errors, trajs

    def convergence_exp(self, pf, steps=[4, 5, 6, 7]):
        nsteps = [2**i for i in steps]
        solns = {}
        trajs = {}
        dts = []

        ref_traj = pf.compute_reference()
        ref_solution = get_final_solution(pf, ref_traj)

        for step in nsteps:
            pf.nsteps = step
            pf.dt = pf.tfinal / pf.nsteps
            dts.append(pf.dt)

            trajectory = pf.run()
            trajs[pf.dt] = trajectory
            solns[pf.dt] = get_final_solution(pf, trajectory)

        errors = [compute_error(solns[dt], ref_solution) for dt in dts]

        slope, _ = get_linear_fit_loglog(dts, errors)

        return solns, dts, errors, trajs


def get_final_solution(pf, trajectory):
    solution = trajectory[(trajectory['rank'] == pf.tasks - 1) & (trajectory[
        'iter'] == pf.iters) & (trajectory['time'] == pf.tfinal) & (trajectory[
            'level'] == pf.levels) & (trajectory['step'] == pf.nsteps)][
                'solution'].values

    return solution


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


def get_linear_fit_loglog(timesteps, errors):
    slope, intercept = np.polyfit(np.log10(timesteps), np.log10(errors), 1)

    return slope, intercept


def plot_errors(timesteps, errors):
    """ logE v logdt """

    plt.scatter(timesteps, errors)
    plt.show()


def compute_error(soln, ref_soln):
    error = np.linalg.norm(soln - ref_soln)  # / np.linalg.norm(ref_soln)

    return error.max()
