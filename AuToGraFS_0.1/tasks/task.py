import sys
import optparse
import traceback
from time import time

import numpy as np

from ase.parallel import world
from ase.visualize import view
from ase.io import read, write
from ase.io import string2index
from ase.constraints import FixAtoms
from ase.optimize.lbfgs import LBFGS
from ase.utils import opencew, devnull, prnt
from ase.tasks.io import read_json, write_json
from ase.data import chemical_symbols, atomic_numbers
from ase.tasks.calcfactory import calculator_factory


class Task:
    taskname = 'generic-task'

    def __init__(self, calcfactory='emt',
                 tag=None, magmoms=None, gui=False,
                 write_summary=False, use_lock_files=False,
                 write_to_file=None, slice=slice(None),
                 logfile='-'):

        """Generic task object.

        This task will do a single of the energy and forces for the
        configurations that subcalsses define in their build_system()
        methods.

        calcfactory: CalculatorFactory object or str
            A calculator factory or the name of a calculator.

        For the meaning of the other arguments, see the add_options()
        and parse_args() methods."""

        self.set_calculator_factory(calcfactory)

        self.tag = tag
        self.magmoms = magmoms
        self.gui = gui
        self.write_summary = write_summary
        self.use_lock_files = use_lock_files
        self.write_to_file = write_to_file
        self.slice = slice

        if world.rank == 0:
            if logfile is None:
                logfile = devnull
            elif isinstance(logfile, str):
                if logfile == '-':
                    logfile = sys.stdout
                else:
                    logfile = open(logfile, 'w')
        else:
            logfile = devnull
        self.logfile = logfile
        
        self.write_funcs = [write_json]
        self.read_func = read_json
    
        self.data = {}  # data read from json files
        self.results = {}  # results from analysis of json files

        self.summary_header = [('name', ''), ('E', 'eV')]

        self.interactive_python_session = False
        self.contains = None
        self.modify = None

    def set_calculator_factory(self, calcfactory):
        if isinstance(calcfactory, str):
            calcfactory = calculator_factory(calcfactory)

        self.calcfactory = calcfactory
        
    def log(self, *args, **kwargs):
        prnt(file=self.logfile, *args, **kwargs)

    def get_filename(self, name=None, ext=''):
        filename = self.taskname + '-' + self.calcfactory.name.lower()
        if self.tag:
            filename += '-' + self.tag
        if name:
            filename = name + '-' + filename
        return filename + ext

    def expand(self, names):
        """Expand ranges like H-Li to H, He, Li."""
        if isinstance(names, str):
            names = [names]
            
        newnames = []
        for name in names:
            if name.count('-') == 1:
                s1, s2 = name.split('-')
                Z1 = atomic_numbers.get(s1)
                Z2 = atomic_numbers.get(s2)
                if Z1 is None or Z2 is None:
                    newnames.append(name)
                else:
                    newnames.extend(chemical_symbols[Z1:Z2 + 1])
            else:
                newnames.append(name)

        return newnames

    def exclude(self, names):
        newnames = []
        for name in names:
            atoms = self.create_system(name)
            if (self.contains is None or
                self.contains in atoms.get_chemical_symbols()):
                newnames.append(name)
        return newnames

    def run(self, names):
        """Run task far all names.

        The task will be one of these four:

        * Open ASE's GUI
        * Write configuration to file
        * Write summary
        * Do the actual calculation
        """

        names = self.expand(names)
        names = names[self.slice]
        names = self.exclude(names)

        if self.gui:
            for name in names:
                view(self.create_system(name))
            return
        
        if self.write_to_file:
            if self.write_to_file[0] == '.':
                for name in names:
                    filename = self.get_filename(name, self.write_to_file)
                    write(filename, self.create_system(name))
            else:
                assert len(names) == 1
                write(self.write_to_file, self.create_system(names[0]))
            return

        if self.write_summary:
            self.read(names)
            self.analyse()
            self.summarize(names)
            return

        atoms = None
        for name in names:
            if self.use_lock_files:
                lockfilename = self.get_filename(name, '.json')
                fd = opencew(lockfilename)
                if fd is None:
                    self.log('Skipping', name)
                    continue
                fd.close()
            atoms = self.run_single(name)
        
        return atoms

    def run_single(self, name):
        try:
            atoms = self.create_system(name)
        except Exception:
            self.log(name, 'FAILED')
            traceback.print_exc(file=self.logfile)
            return
            
        atoms.calc = self.calcfactory(self.get_filename(name), atoms)

        tstart = time()

        try:
            data = self.calculate(name, atoms)
        except KeyboardInterrupt:
            raise
        except Exception:
            self.log(name, 'FAILED')
            traceback.print_exc(file=self.logfile)
            return

        tstop = time()
        data['time'] = tstop - tstart

        for write in self.write_funcs:
            filenamebase = self.get_filename(name)
            write(filenamebase, atoms, data)
        
        return atoms

    def create_system(self, name):
        if '.' in name:
            system = read(name)
        else:
            system = self.build_system(name)

        if self.magmoms is not None:
            system.set_initial_magnetic_moments(
                np.tile(self.magmoms, len(system) // len(self.magmoms)))

        if self.modify:
            exec self.modify

        return system

    def calculate(self, name, atoms):
        e = atoms.get_potential_energy()
        f = atoms.get_forces()
        return {'energy': e, 'forces': f}

    def read(self, names):
        self.data = {}
        for name in names:
            filenamebase = self.get_filename(name)
            try:
                data = self.read_func(filenamebase)
            except (IOError, SyntaxError, ValueError):
                continue
            self.data[name] = data

    def analyse(self):
        for name, data in self.data.items():
            self.results[name] = [data['energy']]

    def summarize(self, names):
        self.log(' '.join('%10s' % x[0] for x in self.summary_header))
        self.log(' '.join('%10s' % x[1] for x in self.summary_header))
        for name in names:
            data = self.results.get(name, [])
            s = '%10s' % name
            for x in data:
                if x is None:
                    s += '           '
                else:
                    s += '%11.3f' % x
            self.log(s)

    def create_parser(self):
        calcname = self.calcfactory.name
        parser = optparse.OptionParser(
            usage='%prog [options] system(s)',
            description='Run %s calculation.' % calcname)
        self.add_options(parser)
        return parser

    def add_options(self, parser):
        general = optparse.OptionGroup(parser, 'General')
        general.add_option('-t', '--tag',
                            help='String tag added to filenames.')
        general.add_option('-M', '--magnetic-moment',
                           metavar='M1,M2,...',
                           help='Magnetic moment(s).  ' +
                           'Use "-M 1" or "-M 2.3,-2.3".')
        general.add_option('-G', '--gui', action='store_true',
                            help="Pop up ASE's GUI.")
        general.add_option('-s', '--write-summary', action='store_true',
                            help='Write summary.')
        general.add_option('--slice', metavar='start:stop:step',
                            help='Select subset of calculations using ' +
                            'Python slice syntax.  ' +
                            'Use "::2" to do every second calculation and ' +
                            '":-5" to do the last five.')
        general.add_option('-w', '--write-to-file', metavar='FILENAME',
                            help='Write configuration to file.')
        general.add_option('-i', '--interactive-python-session',
                            action='store_true',
                            help='Run calculation inside interactive Python ' +
                            'session.  A possible $PYTHONSTARTUP script ' +
                            'will be imported and the "atoms" variable ' +
                            'refers to the Atoms object.')
        general.add_option('-l', '--use-lock-files', action='store_true',
                            help='Skip calculations where the json ' +
                           'lock-file or result file already exists.')
        general.add_option('--contains', metavar='ELEMENT',
                            help='Run only systems containing specific ' +
                           'element.')
        general.add_option('--modify', metavar='...',
                            help='Modify system with Python statement.  ' +
                           'Example: "system.positions[-1,2]+=0.1". ' +
                           'Warning: no spaces allowed!')
        parser.add_option_group(general)
    
    def parse_args(self, args=None):
        if args is None:
            args = sys.argv[1:]

        parser = self.create_parser()
        self.calcfactory.add_options(parser)
        opts, args = parser.parse_args(args)

        if len(args) == 0:
            parser.error('incorrect number of arguments')

        self.parse(opts, args)
        self.calcfactory.parse(opts, args)

        return args

    def parse(self, opts, args):
        if opts.tag:
            self.tag = opts.tag
            
        if opts.magnetic_moment:
            self.magmoms = np.array(
                [float(m) for m in opts.magnetic_moment.split(',')])
        
        self.gui = opts.gui
        self.write_summary = opts.write_summary
        self.write_to_file = opts.write_to_file
        self.use_lock_files = opts.use_lock_files
        self.interactive_python_session = opts.interactive_python_session
        self.contains = opts.contains
        self.modify = opts.modify

        if opts.slice:
            self.slice = string2index(opts.slice)


class OptimizeTask(Task):
    taskname = 'opt'

    def __init__(self, fmax=None, constrain_tags=[], **kwargs):
        self.fmax = fmax
        self.constrain_tags = constrain_tags

        Task.__init__(self, **kwargs)
        
        self.summary_header.append(('E-E0', 'eV'))

    def optimize(self, name, atoms):
        mask = [t in self.constrain_tags for t in atoms.get_tags()]
        if mask:
            constrain = FixAtoms(mask=mask)
            atoms.constraints = [constrain]

        optimizer = LBFGS(atoms, trajectory=self.get_filename(name, '.traj'),
                          logfile=None)
        optimizer.run(self.fmax)
        
    def calculate(self, name, atoms):
        data = Task.calculate(self, name, atoms)

        if self.fmax is not None:
            self.optimize(name, atoms)

            data['minimum energy'] = atoms.get_potential_energy()
            data['minimum forces'] = atoms.get_forces()
        
        return data

    def analyse(self):
        Task.analyse(self)
        for name, data in self.data.items():
            if 'minimum energy' in data:
                self.results[name].append(data['energy'] -
                                          data['minimum energy'])
            else:
                self.results[name].append(None)

    def add_options(self, parser):
        Task.add_options(self, parser)

        optimize = optparse.OptionGroup(parser, 'Optimize')
        optimize.add_option('-R', '--relax', type='float', metavar='FMAX',
                            help='Relax internal coordinates using L-BFGS '
                            'algorithm.')
        optimize.add_option('--constrain-tags', type='str',
                            metavar='T1,T2,...',
                            help='Constrain atoms with tags T1, T2, ...')
        parser.add_option_group(optimize)

    def parse(self, opts, args):
        Task.parse(self, opts, args)

        self.fmax = opts.relax

        if opts.constrain_tags:
            self.constrain_tags = [int(t)
                                   for t in opts.constrain_tags.split(',')]
