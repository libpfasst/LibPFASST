import re
from os import getcwd, chdir, makedirs, devnull
from subprocess import call
from scipy.io import FortranFile
from pfasst import PFASST
import numpy as np

class NWC(PFASST):

    def __init__(self, home, **kwargs):

        PFASST.__init__(self, home, **kwargs)
        print 'Setting up NWC'
        self.params = kwargs
        self.nwc = self.home + '/nwc.sh'
        self.exact_dir = self.params['exact_dir']
        self.cwd = getcwd() + '/'
        self.nwc_base_string = 'echo\nscratch_dir ./scratch\npermanent_dir ./perm\n\n'\
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
        nwc_inp = self.nwc_base_string.format(geometry, self.params['basis'],
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


