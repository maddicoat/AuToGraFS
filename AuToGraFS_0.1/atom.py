"""This module defines the Atom object."""

import warnings

import numpy as np

from ase.data import atomic_numbers, chemical_symbols, atomic_masses


#         singular,    plural,     default value
names = {'position':        ('positions',           np.zeros(3)),
         'number':          ('numbers',             0),
         'tag':             ('tags',                0),
         'momentum':        ('momenta',             np.zeros(3)),
         'mass':            ('masses',              None),
         'magmom':          ('magmoms',             0.0),
         'charge':          ('charges',             0.0),
         'original_index':  ('original_indices',    0),
         'current_index':   ('current_indices',     0),
         'fragmentID':      ('fragmentIDs',         0),
         'mmtype':          ('mmtypes',             None),
         'bondlist':        ('bondlists',           None)
         }


def atomproperty(name, doc):
    """Helper function to easily create Atom attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def xyzproperty(index):
    """Helper function to easily create Atom XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')


class Atom(object):
    """Class for representing a single atom.

    Parameters:
    
    symbol: str or int
        Can be a chemical symbol (str) or an atomic number (int).
    position: sequence of 3 floats
        Atomi position.
    tag: int
        Special purpose tag.
    momentum: sequence of 3 floats
        Momentum for atom.
    mass: float
        Atomic mass in atomic units.
    magmom: float or 3 floats
        Magnetic moment.
    charge: float
        Atomic charge.
    """
    __slots__ = ['data', 'atoms', 'index']

    def __init__(self, symbol='X', position=(0, 0, 0),
                 tag=None, momentum=None, mass=None,
                 magmom=None, charge=None,
                 original_index=None, current_index=None,
                 fragmentID=None,mmtype=None, bondlist=None,
                 atoms=None, index=None): 

        self.data = d = {}

        if atoms is None:
            # This atom is not part of any Atoms object:
            if isinstance(symbol, str):
                d['number'] = atomic_numbers[symbol]
            else:
                d['number'] = symbol
            d['position'] = np.array(position, float)
            d['tag'] = tag
            if momentum is not None:
                momentum = np.array(momentum, float)
            d['momentum'] = momentum
            d['mass'] = mass
            if magmom is not None:
                magmom = np.array(magmom, float)
            d['magmom'] = magmom
            d['charge'] = charge
            d['original_index'] = original_index
            d['current_index'] = current_index
            d['fragmentID'] = 0
            d['mmtype'] = mmtype
            d['bondlist'] = {}

        self.index = index
        self.atoms = atoms

    def __repr__(self):
        s = "Atom('%s', %s" % (self.symbol, list(self.position))
        for name in ['tag', 'momentum', 'mass', 'magmom', 'charge', 'mmtype', 'original_index', 'current_index', 'fragmentID','bondlist']:
            value = self.get_raw(name)
            if value is not None:
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                s += ', %s=%s' % (name, value)
        if self.atoms is None:
            s += ')'
        else:
            s += ', index=%d)' % self.index
        return s

    def cut_reference_to_atoms(self):
        """Cut reference to atoms object."""
        for name in names:
            self.data[name] = self.get_raw(name)
        self.index = None
        self.atoms = None
        
    def get_raw(self, name):
        """Get attribute, return None if not explicitely set."""
        if name == 'symbol':
            return chemical_symbols[self.get_raw('number')]

        if self.atoms is None:
            return self.data[name]
        
        plural = names[name][0]
        if plural in self.atoms.arrays:
            return self.atoms.arrays[plural][self.index]
        else:
            return None

    def get(self, name):
        """Get attribute, return default if not explicitely set."""
        value = self.get_raw(name)
        if value is None:
            if name == 'mass':
                value = atomic_masses[self.number]
            else:
                value = names[name][1]
        return value

    def set(self, name, value):
        """Set attribute."""
        if name == 'symbol':
            name = 'number'
            value = atomic_numbers[value]

        if self.atoms is None:
            assert name in names
            self.data[name] = value
        else:
            plural, default = names[name]
            if plural in self.atoms.arrays:
                array = self.atoms.arrays[plural]
                if name == 'magmom' and array.ndim == 2:
                    assert len(value) == 3
                array[self.index] = value
            else:
                if name == 'magmom' and np.asarray(value).ndim == 1:
                    array = np.zeros((len(self.atoms), 3))
                elif name == 'mass':
                    array = self.atoms.get_masses()
                else:
                    default = np.asarray(default)
                    array = np.zeros((len(self.atoms),) + default.shape,
                                     default.dtype)
                array[self.index] = value
                self.atoms.new_array(plural, array)

    def delete(self, name):
        """Delete attribute."""
        assert self.atoms is None
        assert name not in ['number', 'symbol', 'position']
        self.data[name] = None

    symbol = atomproperty('symbol', 'Chemical symbol')
    number = atomproperty('number', 'Atomic number')
    position = atomproperty('position', 'XYZ-coordinates')
    tag = atomproperty('tag', 'Integer tag')
    momentum = atomproperty('momentum', 'XYZ-momentum')
    mass = atomproperty('mass', 'Atomic mass')
    magmom = atomproperty('magmom', 'Initial magnetic moment')
    charge = atomproperty('charge', 'Atomic charge')
    original_index = atomproperty('original_index', 'Index as read in from fragment')
    current_index = atomproperty('current_index', 'Current index of atom in supermolecule')
    fragmentID = atomproperty('fragmentID', 'I am in supermolecule. What fragment wad I originally from?')
    mmtype = atomproperty('mmtype', 'MM Type')
    bondlist = atomproperty('bondlist', 'Dict of bonded atoms and orders')
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)

    def _get(self, name):
        """Helper function for deprecated get methods."""
        warnings.warn('Use atom.%s' % name, stacklevel=3)
        return getattr(self, name)

    def _set(self, name, value):
        """Helper function for deprecated set methods."""
        warnings.warn('Use atom.%s = ...' % name, stacklevel=3)
        setattr(self, name, value)

    def get_symbol(self): return self._get('symbol')
    def get_atomic_number(self): return self._get('number')
    def get_position(self): return self._get('position')
    def get_tag(self): return self._get('tag')
    def get_momentum(self): return self._get('momentum')
    def get_mass(self): return self._get('mass')
    def get_initial_magnetic_moment(self): return self._get('magmom')
    def get_charge(self): return self._get('charge')
    def get_original_index(self): return self._get('original_index')
    def get_current_index(self): return self._get('current_index')
    def get_fragmentID(self): return self._get('fragmentID')
    def get_mmtype(self): return self._get('mmtype')
    def get_bondlist(self): return self._get('bondlist')

    def set_symbol(self, value): self._set('symbol', value)
    def set_atomic_number(self, value): self._set('number', value)
    def set_position(self, value): self._set('position', value)
    def set_tag(self, value): self._set('tag', value)
    def set_momentum(self, value): self._set('momentum', value)
    def set_mass(self, value): self._set('mass', value)
    def set_initial_magnetic_moment(self, value): self._set('magmom', value)
    def set_charge(self, value): self._set('charge', value)
    def set_original_index(self, value): self._set('original_index', value)
    def set_current_index(self, value): self._set('current_index', value)
    def set_fragmentID(self, value): self._set('fragmentID', value)
    def set_mmtype(self, value): self._set('mmtype', value)
    def set_bondlist(self, value): self._set('bondlist', value)
