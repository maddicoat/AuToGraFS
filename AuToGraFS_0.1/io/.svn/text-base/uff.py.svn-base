from ase.atoms import Atoms
from ase.parallel import paropen 


def read_uff(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    tag_signifier = 'Data'
    mStart = 'GEOMETRY CARTESIAN'
    mStop = 'END'
    ignoreLines = True
    lines = fileobj.readlines()
    natoms=0

    images = []
    symbols = []
    positions = []
    mmtypes = []
    original_indices = []
    current_indices = []
    bondlists = []
    tags = []
    info={}
    for line in lines:
        if tag_signifier in line:
            print line
            (pre,val) = line.split('=')
            (junk,key)=pre.split()
            key.strip()
            info[key.lower()] = val.strip()
        if mStart in line:
            ignoreLines= False
            continue
        if mStop in line:
            ignoreLines= True
        if not ignoreLines:
            #print line
            natoms+=1
            tmp = line.split()
            if len(tmp)==8:
                symbol, x, y, z, mm, qmmm, bonding, tag = line.split()
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
                mm = mm.rsplit('=',1)
                mmtype=mm[-1]
                mmtypes.append(mmtype)
                original_indices.append(natoms)
                current_indices.append(natoms)
                bonding=bonding.rsplit('=',1)
                bonding=bonding[-1]
                bonddict = dict((int(k.strip()), float(v.strip())) for k,v in
                        (item.split('/') for item in bonding.split(':')))
                bondlists.append(bonddict)
                #print bonddict
                tags.append(tag)
            elif len(tmp)==7:
                symbol, x, y, z, mm, junk, bonding = line.split()
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
                mm = mm.rsplit('=',1)
                mmtype=mm[-1]
                mmtypes.append(mmtype)
                original_indices.append(natoms)
                current_indices.append(natoms)
                bonding=bonding.rsplit('=',1)
                bonding=bonding[-1]
                #print tmp
                bonddict = dict((int(k.strip()), float(v.strip())) for k,v in
                        (item.split('/') for item in bonding.split(':')))
                bondlists.append(bonddict)
                #print bonddict
                tags.append(None)
            else:
                print "Line not parseable, wrong number of things"
                print line
                sys.exit()
    #print mmtypes
    print
    #images.append(Atoms(symbols=symbols, positions=positions))
    images.append(Atoms(symbols=symbols, positions=positions, original_indices=original_indices, current_indices=current_indices, mmtypes=mmtypes, bondlists=bondlists, tags=tags, info=info))
    return images[index]

            #print symbols
            #print len(tmp)
            #print "Number of atoms =",natoms
            #symbol, x, y, z, mm, junk, bonding, tag = line.split()
            #print positions
            #print mmtypes
            #print bondlists
            #print tags


def write_uff(fileobj, images):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]  

    

    #Add in something here to read standard headers from somewhere

    fileobj.write('%s\n' % 'GEOMETRY CARTESIAN')
    #for 1 ..N in current indices
    symbols = images[0].get_chemical_symbols()
    
    #stringify the bonding dicts. Remember they're stored in an numpy.ndarray (no iterkeys), so we do it the long way around
    bonding=images[0].get_bondlists()
    bondstrings=[]
    for bl in bonding:
        #tmp=dictjoin(bl,'/',':')
        bondstrings.append(dictjoin(bl,'/',':'))

    for atoms in images:
        for s, (x, y, z), m, b in zip(symbols, atoms.get_positions(), atoms.get_mmtypes(), bondstrings):
            fileobj.write('%-2s %15.8f %15.8f %15.8f %-11s%-5s %13s%-50s\n' % (s, x, y, z ,'    MMTYPE=',m, 'QMMM=MM BOND=',b))

    if images[0].pbc.any():
        fileobj.write('%s\n' % 'Periodic general')
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

    fileobj.write('%s\n' % 'END')

def genkvs(d, keys, joiner):
    for key in keys:
        yield '%s%s%s' % (key, joiner, d[key])


def dictjoin(_dict, joiner, sep):
    keys = sorted(_dict.iterkeys())
    return sep.join(genkvs(_dict, keys, joiner))

