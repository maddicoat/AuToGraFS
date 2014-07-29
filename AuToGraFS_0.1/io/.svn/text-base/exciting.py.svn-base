"""
This is the Implementation of the exciting I/O Functions
The functions are called with read write usunf the format "exi"

The module depends on lxml  http://codespeak.net/lxml/
"""
from math import pi, cos, sin, sqrt, acos

import numpy as np

from ase.atoms import Atoms
from ase.parallel import paropen
from ase.units import Bohr


def read_exciting(fileobj, index=-1):
    """Reads structure from exiting xml file.
    
    Parameters
    ----------
    fileobj: file object
        File handle from which data should be read.
        
    Other parameters
    ----------------
    index: integer -1
        Not used in this implementation.
    """
    from lxml import etree as ET
    #parse file into element tree
    doc = ET.parse(fileobj)
    root = doc.getroot()
    speciesnodes = root.find('structure').getiterator('species') 
 
    symbols = []
    positions = []
    basevects = []
    atoms = None
    #collect data from tree
    for speciesnode in speciesnodes:
        symbol = speciesnode.get('speciesfile').split('.')[0]
        natoms = speciesnode.getiterator('atom')
        for atom in natoms:
            x, y, z = atom.get('coord').split()
            positions.append([float(x), float(y), float(z)])
            symbols.append(symbol)
    # scale unit cell accorting to scaling attributes
    if doc.xpath('//crystal/@scale'):
        scale = float(str(doc.xpath('//crystal/@scale')[0]))
    else:
        scale = 1
        
    if doc.xpath('//crystal/@stretch'):
        a, b, c = doc.xpath('//crystal/@scale')[0].split()
        stretch = np.array([float(a),float(b),float(c)])
    else:    
        stretch = np.array([1.0, 1.0, 1.0])
    basevectsn = doc.xpath('//basevect/text()') 
    for basevect in basevectsn:
        x, y, z = basevect.split()
        basevects.append(np.array([float(x) * Bohr * stretch[0],
                                   float(y) * Bohr * stretch[1], 
                                   float(z) * Bohr * stretch[2]
                                   ]) * scale)
    atoms = Atoms(symbols=symbols, cell=basevects)
 
    atoms.set_scaled_positions(positions)
    if 'molecule' in root.find('structure').attrib.keys():
        if root.find('structure').attrib['molecule']:
            atoms.set_pbc(False)
    else:
        atoms.set_pbc(True)
        
    return atoms

def write_exciting(fileobj, images):
    """writes exciting input structure in XML
    
    Parameters
    ----------
    fileobj : File object
        Filehandle to which data should be written
    images : Atom Object or List of Atoms objects
        This function will write the first Atoms object to file 
    
    Returns
    -------
    """
    from lxml import etree as ET
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')
    root = atoms2etree(images)
    fileobj.write(ET.tostring(root, method='xml', 
                              pretty_print=True,
                              xml_declaration=True,
                              encoding='UTF-8'))

def atoms2etree(images):
    """This function creates the XML DOM corresponding
     to the structure for use in write and calculator
    
    Parameters
    ----------
    
    images : Atom Object or List of Atoms objects
        This function will create a 
    
    Returns
    -------
    root : etree object
        Element tree of exciting input file containing the structure
    """
    from lxml import etree as ET
    if not isinstance(images, (list, tuple)):
        images = [images]

    root = ET.Element('input')
    title = ET.SubElement(root, 'title')
    title.text = ''
    structure = ET.SubElement(root, 'structure')
    crystal= ET.SubElement(structure, 'crystal')
    atoms = images[0]
    for vec in atoms.cell:
        basevect = ET.SubElement(crystal, 'basevect')
        basevect.text = '%.14f %.14f %.14f' % tuple(vec / Bohr)
                             
    species = {}
    symbols = []
    for aindex, symbol in enumerate(atoms.get_chemical_symbols()):
        if symbol in species:
            species[symbol].append(aindex)
        else:
            species[symbol] = [aindex]
            symbols.append(symbol)
    scaled = atoms.get_scaled_positions()
    for symbol in symbols:
        speciesnode = ET.SubElement(structure, 'species',
                                    speciesfile='%s.xml' % symbol,
                                    chemicalSymbol=symbol)
        for a in species[symbol]:
           atom = ET.SubElement(speciesnode, 'atom',
                                coord='%.14f %.14f %.14f' % tuple(scaled[a]))
    return root
