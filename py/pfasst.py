import matplotlib
matplotlib.use('pdf')

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
from os import getcwd, chdir, makedirs, devnull
from subprocess import call, check_output
from scipy.interpolate import interp1d
from scipy.io import FortranFile


def create_argparser():
    parser = argparse.ArgumentParser(
        description='From nml file for PFASST, generate exact/ref solutions')
    parser.add_argument(
        'pfasst_input', metavar='inp', type=str, help='PFASST input filename')
    parser.add_argument('--tfinal', metavar='tfinal', type=float)
    parser.add_argument('--dt', metavar='dt', type=float)
    parser.add_argument('--nsteps', metavar='nsteps', type=float)
    parser.add_argument('--nodes', metavar='nodes', type=int)

    args = parser.parse_args()
    return args


class PFASST():
    def __init__(self, home, pfinp, **kwargs):
        self.pfasst = home + 'dev/main.exe'

        with open(pfinp, 'r') as f:
            content = f.read()
        self.base_string = content

        self.input = pfinp
        self.params = kwargs
        self.nprob = int(re.search(r'nprob\ =\ ([0-9])', content).group(1))

        self._set_optional_params()
        self.params = self._save_params()
        self._create_pf_string()

    def _save_params(self):
        params = {}
        params['tfinal'] = self.tfinal
        params['nsteps'] = self.nsteps
        params['dt'] = self.dt
        params['nodes'] = self.nnodes

        return params

    def _create_pf_string(self):
        for i in ['nnodes', 'Tfin', 'nsteps']:
            self.base_string = re.sub(r''+i+'\ =\ (.+)', \
                                    r''+i+' = {}', self.base_string)

    def _set_optional_params(self):
        if self.params['magnus']:
            self.magnus = int(self.params['magnus'])
        else:
            self.magnus = int(
                re.search(r'magnus_order\ = (.+)', self.base_string).group(1).split()[self.level])

        if self.params['nodes']:
            self.nnodes = int(self.params['nodes'])
        else:
            self.nnodes = int(
                re.search(r'nnodes\ = (.+)', self.base_string).group(1).split()[self.level])

        if self.params['tfinal']:
            self.tfinal = float(self.params['tfinal'])
        else:
            self.tfinal = float(re.search(r'Tfin\ =\ (.+)', self.base_string).group(1))

        if self.params['nsteps']:
            self.nsteps = int(self.params['nsteps'])
        else:
            self.nsteps = int(re.search(r'nsteps\ =\ (.+)', self.base_string).group(1))

        if self.params['dt']:
            self.dt = self.params['dt']
        else:
            self.dt = self.tfinal / self.nsteps

    def write_to_file(self):
        pfinp = self.base_string.format(self.nnodes, self.tfinal, int(self.nsteps))
        with open(self.input, 'w') as f:
            f.write(pfinp)

    def print_params(self):
        params = {}
        params['tfinal'] = self.tfinal
        params['nsteps'] = self.nsteps
        params['dt'] = self.dt
        params['nodes'] = self.nnodes
        params['magnus'] = self.magnus

        pprint(params, width=1)

    def run(self):
        print '---- running pfasst ----'
        pf.print_params()
        output = check_output([self.pfasst, pf.input])

        step_infos = re.findall(r'resid\:\ .+', output)
        error = float(step_infos[-2].split(':')[-1])

        f = FortranFile('final_solution')
        final_solution = f.read_record(np.complex_)
        f.close()

        dim = np.sqrt(final_solution.shape)
        final_solution = np.reshape(final_solution, (dim,dim))
        return final_solution

    def compute_reference(self):
        self.nsteps = 2**16

        self.dt = pf.tfinal/pf.nsteps

        self.write_to_file()
        ref_soln = pf.run()

        return ref_soln

def plot_errors(timesteps, errors):
    """ logE v logdt """

    dt = np.log(timesteps)
    err = np.log(errors)

    m, b = np.polyfit(dt, err, 1)

    print 'slope = {}, intercept = {}'.format(m,b)
    fig = plt.figure()
    ax = plt.gca()
    ax.scatter(dt, err)
    plt.savefig('conv')

def compute_error(soln, ref_soln):
    error = np.linalg.norm(soln-ref_soln)/np.linalg.norm(ref_soln)

    return error

if __name__ == '__main__':
    args = create_argparser()
    projhome = '/home/bkrull/devel/pfasst-nwchem'
    pf = PFASST(projhome, args.pfasst_input, **vars(args))

    FNULL = open(devnull, 'w')
    pf.ref = pf.compute_reference()

    stepsizes = [2**i for i in [4, 5, 6, 7]]
    solns = {}
    dts = []

    for stepsize in stepsizes:
        pf.nsteps = stepsize
        pf.dt = pf.tfinal/pf.nsteps
        dts.append(pf.dt)

        pf.write_to_file()
        soln = pf.run()
        solns[pf.dt] = soln

    FNULL.close()

    errors = []
    for dt in dts:
        error = compute_error(solns[dt], pf.ref)
        print 'nodes={}, dt={}, error={}'.format(pf.nnodes, dt, error)
        errors.append(error)

    plot_errors(dts, errors)
