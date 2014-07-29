#!/usr/bin/python
from ase.atom import Atom
from ase.atoms import Atoms
from ase.align import *
from ase.topology import make_model
from ase.io import read, write
from collections import deque
import argparse

parser = argparse.ArgumentParser(description="Generate a defined framework structure.")
parser.add_argument("-c", "--controlfile", nargs="?", type=str, default="control.txt", 
                   help="control file for the mof generator, default=control.txt")
parser.add_argument("-o", "--output-prefix", nargs="?", type=str, default="mof", 
                   help="output file for generated mof, default=mof")
parser.add_argument("-p", "--path", nargs="?", type=str, default="/home/maddicoat/src/ase/ase/database/", 
                   help="base path for centers, linkers and functional groups, default=/home/maddicoat/src/ase/ase/database/")
parser.add_argument("-e", "--output-extn", type=str, default="gin", 
                   help="extn (format) for the generated mof, default=gin")
parser.add_argument("-s", "--supercell", type=int, nargs="*", default=None, 
                   help="optional production of supercell. Default is None. Argument may be one integer or three.")
parser.add_argument('--leave_ghosts',dest='ghosts',action='store_true',
                   help="Leave extra ghost (Bq) atoms in structure for post-functionalisation. Default is False.")
parser.add_argument('--no-leave_ghosts',dest='ghosts',action='store_false',
                   help="Removes extra ghost (Bq) atoms in structure. This is the default. To leave them, use the --leave_ghosts option.")
parser.set_defaults(ghosts=False)


args = parser.parse_args()
#print args
output=args.output_prefix+'.'+args.output_extn
xyz_output=args.output_prefix+'.xyz'
cell_dim=args.supercell

if args.path:
    if not args.path.endswith("/"):
        args.path = args.path+'/'

#functionalise
def functionalise():
    for  frag_id,frag in fragments.iteritems():
        if "func" in frag.info:
            for anum,fgrp in frag.info["func"].iteritems():
                #print anum,fgrp
                filename = fgrp+'.inp'
                try:
                    this_fgrp=read(filename)
                except IOError:
                    path=args.path+fgrppath+filename #changed from basepath
                    this_fgrp=read(path)
                #print this_fgrp
                this_fgrp.fragmentIDs=frag_id
                mof_replace_index = mof.get_atom_index(frag_id,anum)
                #print "MRI = ", mof_replace_index
                fgrp_replace_index = this_fgrp.get_atom_indices('X') 
                if len(fgrp_replace_index) != 1:
                    raise RuntimeError('Functional group must contain only one dummy atom!')
                else:
                    fgrp_replace_index = fgrp_replace_index[0]
                #print replace_index
                mof.add_group(this_fgrp,mof_replace_index,fgrp_replace_index)

#Now set up some default locations for things. 
basepath="/home/maddicoat/src/ase/ase/database/"
modelpath="models/"
centerpath="centers/"
linkerpath="linkers/"
pillarpath="pillars/"
fgrppath="functional_groups/"

control = open(args.controlfile)
lines = control.readlines()
mol_spec=[]
frag_types=[]
info_spec=[]
#first line is either requesting a model or a topology. Handling is different in each case
if 'topology' in lines[0]:
    mol_ids = deque()
    mol_types = deque()
    for line in lines:
        (name,val) = line.split('=',1)
        name = name.strip().lower()
        val = val.strip() 
        if name == 'topology':
            topology = val
            #modelfile = val+'.model' #will be created
        elif name == 'center':
            mol_ids.appendleft(val)
            mol_types.appendleft(name)
        elif name == 'linker':
            mol_ids.append(val)
            mol_types.append(name)
        elif name =='supercell': # supercell read in here will carry but won't be written to new file
            cell_dim=[]
            if not "," in val:
                x = int(val.strip())
                cell_dim=[x,x,x]
            else:
                val.strip()
                (x,y,z) = val.split(',') 
                cell_dim = map(int,[x,y,z])
    make_model(topology,args.path,mol_types,mol_ids)        #args.path
    #Now we read in the newly created control file
    control = open('control-mofgen.txt')
    lines = control.readlines()

for line in lines:
    (name,val) = line.split('=',1)
    name = name.strip().lower()
    val = val.strip() #val may still have other stuff attached to it
    if name =='model':
        modelfile = val+'.model'
    elif name == 'center' or name == 'linker':
        frag_types.append(name)
        #do we have other stuff attached here?
        info_flag = len(val.split(' '))
        print "info flag:", info_flag
        if info_flag <= 1:
            mol_spec.append(val)
            info_spec.append(None)
        else:
            sbu_spec=[]
            sbu_spec=val.split(' ')
            mol_spec.append(sbu_spec.pop(0))
            info_spec.append(sbu_spec)
    elif name =='supercell':
        cell_dim=[]
        if not "," in val:
            x = int(val.strip())
            cell_dim=[x,x,x]
        else:
            val.strip()
            (x,y,z) = val.split(',') 
            cell_dim = map(int,[x,y,z])


#modelfile
try:
    model=read(modelfile)
except IOError:
    path=basepath+modelpath+modelfile
    model=read(path)

for slot in model:
    these_tags = model[slot].get_tags()



#controlfile
fragments={}
for obj in xrange(0, len(mol_spec)):
    name = obj
    filename = mol_spec[obj]+'.inp'
    try:
        fragments[name]=read(filename)
        fragments[name].fragmentIDs=obj
    except IOError:
        path=args.path+frag_types[obj]+'s/'+filename #changed from basepath
        fragments[name]=read(path)
        fragments[name].fragmentIDs=obj
    if info_spec[obj]:
        for attr in info_spec[obj]:
            #We may have either keywords with value (i.e. rotate=90), keywords with a dict (i.e. func=10,NH3:9,F): or bare keywords, which we can set to true
            option_type= len(attr.split('='))
            if option_type == 1: #binary
                fragments[name].info[attr.strip().lower()] = True
            if option_type == 2:
                (k,v) =attr.split('=')  
                #Test if val contains "," #even in the case of only one functionalisation, the list will contain a comma
                if "," not in v:
                    fragments[name].info[k.lower()] = float(v)
                else:
                    funcdict = dict((int(anum.strip()), fgrp.strip()) for anum,fgrp in
                                                (item.split(',') for item in v.split(':')))
                    fragments[name].info[k.lower()] = funcdict
    align(fragments[name],model[name])


test_mol = Atoms()
for obj in xrange(0, len(mol_spec)):
    test_mol += fragments[obj]

write('test_mol.xyz',test_mol)


#if not args.supercell:
if cell_dim is None: 
    mof=assemble(model,fragments)
    #write("assembled.xyz",mof)
    print mof.cell
    pbc=[]
    for x in range(3):
        if any(mof.get_cell()[x]):
            pbc.append(True)
        else:
            pbc.append(False)
    if model[0].info['build'].strip() == 'exact':
        cell=frame_cell(model,fragments)
        mof.set_cell(cell)
        mof.set_pbc(pbc)
        mof = bond_frame(mof,periodic=True)
    elif model[0].info['build'].strip() == 'systre':
        if 'extra_dummies' in model[0].info: #extra dummies only occurs with systre build - e.g. MIL-53
            mof.set_pbc(pbc)
            bond_extra_dummies(mof)
            mof = bond_frame(mof,periodic=True)
        else:
            mof.set_pbc(pbc)
            mof = bond_frame(mof,periodic=True)
    functionalise()
    if not args.ghosts:
        remove_banquo(mof)
    write(output,mof)    
    write(xyz_output,mof)    
    #for x in range(3):
    #    print mof.get_cell()[x]
    #print mof.get_pbc()
    #write('mof.cif',mof)
    frags_output='frags.out'
    f = open(frags_output, 'w')
    for n, o, s, (x, y, z), i in zip(mof.get_current_indices(),mof.get_original_indices(),mof.get_chemical_symbols(), mof.get_positions(), mof.get_fragmentIDs()):
        f.write('%-5d %-5d %-2s %15.8f %15.8f %15.8f %6d\n' % (n, o, s, x, y, z , i))
else:
    mof = assemble(model,fragments)
    pbc=[]
    for x in range(3):
        if any(mof.get_cell()[x]):
            pbc.append(True)
        else:
            pbc.append(False)
    if model[0].info['build'].strip() == 'exact':
        cell=frame_cell(model,fragments)
        mof.set_cell(cell)
        mof.set_pbc(pbc)
    elif model[0].info['build'].strip() == 'systre':
        mof.set_pbc(mof.pbc)
    functionalise()
    if not args.ghosts:
        remove_banquo(mof)
    bigmof = make_supercell(mof,cell_dim)
    write(output,bigmof)    
    write(xyz_output,bigmof)    
    frags_output='frags.out'
    f = open(frags_output, 'w')
    for n, o, s, (x, y, z), i in zip(bigmof.get_current_indices(),bigmof.get_original_indices(),bigmof.get_chemical_symbols(), bigmof.get_positions(), bigmof.get_fragmentIDs()):
        f.write('%-5d %-5d %-2s %15.8f %15.8f %15.8f %6d\n' % (n, o, s, x, y, z , i))
 
 



