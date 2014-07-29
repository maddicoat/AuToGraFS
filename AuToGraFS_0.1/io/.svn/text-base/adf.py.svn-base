from ase.atoms import Atoms
from ase.parallel import paropen


def write_adf(fileobj, images):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]


    fileobj.write('%s\n' % 'Units')
    fileobj.write('%s\n' % '    length Angstrom')
    fileobj.write('%s\n' % 'End')
        
    fileobj.write('%s\n' % 'System')
    fileobj.write('%s\n' % '    Atoms')

    symbols = images[0].get_chemical_symbols()

    for atoms in images:
        for s, (x, y, z), in zip(symbols, atoms.get_positions()):
            fileobj.write('%-8s %-4s %15.8f %15.8f %15.8f \n' % ('        ',s,x, y, z))
    
    fileobj.write('%s\n' % '    End')
    fileobj.write('%s\n' % '    Bonds')

    #write the bonding
    if any(images[0].get_bondlists()):
        bonding=images[0].get_bondlists()
        for index,b in enumerate(bonding):
            for k,v in b.items():
                fileobj.write('%s %-3d %-3d %3.2f\n' % ('        ', index+1,k,v))
    fileobj.write('%s\n' % '    End')
    if images[0].pbc.any():
        fileobj.write('%s\n' % '    Lattice')
        fileobj.write('%8.3f %8.3f %8.3f \n' %
        (images[0].get_cell()[0][0],
         images[0].get_cell()[0][1],
         images[0].get_cell()[0][2]))
        fileobj.write('%8.3f %8.3f %8.3f \n' %
        (images[0].get_cell()[1][0],
         images[0].get_cell()[1][1],
         images[0].get_cell()[1][2]))
        fileobj.write('%8.3f %8.3f %8.3f \n' %
        (images[0].get_cell()[2][0],
         images[0].get_cell()[2][1],
         images[0].get_cell()[2][2]))
        fileobj.write('%s\n' % '    End')
    fileobj.write('%s\n' % 'End')
    fileobj.write('%s\n' % '')
