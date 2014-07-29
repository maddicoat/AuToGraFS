#!/usr/bin/python
from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from ase.lattice.spacegroup import crystal
from ase.calculators.neighborlist import *
from ase.data import *
from math import pi
import numpy as np

def make_Td(names,sizes):

    #work out the expansion factor
    factor = sum(sizes)
    #factor = 3.0
        
    Td = Atoms()
    Td.append(Atom('F', position=(0,0,0)))
    Td.append(Atom('N', position=(0,0,factor)))
    #Atom 3
    Td.append(Atom('N', position=(1.0,0,0.0)))
    Td.set_angle([1,0,2],1.91062939)
    Td.set_distance(0,2, factor, fix=0)
    #Atom 4
    Td.append(Atom('N', position=(1.0,0,0.0)))
    Td.set_distance(0,3, factor, fix=0)
    Td.set_angle([1,0,3],1.91062939, mask=None)
    Td.set_dihedral([2,1,0,3],2.0943951, mask=None)
    #print Td.get_angle([1,0,2])
    #Atom 5
    Td.append(Atom('N', position=(1.0,0,0)))
    Td.set_angle([1,0,4],1.91062939)
    Td.set_dihedral([2,1,0,4],-2.0943951)
    Td.set_distance(0,4, factor, fix=0)
    
    #find side centroids
    tmp1 = Td[1:4]
    coside = tmp1.get_center_of_mass()
    Td.append(Atom('C', position = coside))
    #
    tmp1 = Td[2:5]
    coside = tmp1.get_center_of_mass()
    Td.append(Atom('C', position = coside))
    #
    tmp1 = Td[3:5]
    tmp1.append(Td[1])
    coside = tmp1.get_center_of_mass()
    Td.append(Atom('C', position = coside))
    #
    tmp1 = Td[1:3]
    tmp1.append(Td[4])
    coside = tmp1.get_center_of_mass()
    Td.append(Atom('C', position = coside))
    
    write('td.xyz',Td)
    
    # Now we make tags etc
    
    ## C -> triangles (faces), N => triangle caps
    
    eps = 0.05
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0
    
    #First detect global bonding
    cov_rad=[]
    for atom in Td:
        #cov_rad.append(covalent_radii[atom.number])
        cov_rad.append(factor / 2)
    
    nlist = NeighborList(cov_rad,skin=0.1,self_interaction=False,bothways=True)
    nlist.build(Td)
    
    nbond = 1
    bond_matrix = np.zeros( (len(Td),len(Td)) )
    pbc_bond_matrix = np.zeros( (len(Td),len(Td)) )
    
    bond_dict = {}
    
    for atom in Td:
        print "---------------------"
        indices, offsets = nlist.get_neighbors(atom.index)
        for index,offset in zip(indices,offsets):
            print atom.index, index, atom.symbol, Td[index].symbol, offset, Td.get_distance(atom.index,index,mic=True)
            if atom.index < index:
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                print bond
                bond_dict[bond] = nbond
                #Now we need to find the same bond the other direction to get the offsets right
                indices2,offsets2 = nlist.get_neighbors(index)
                for i2,o2 in zip(indices2,offsets2):
                    if i2 == atom.index:
                        #print "sum of offsets = ", offset, o2, sum(offset + o2)
                        if sum(offset + o2) == 0: #same bond
                            this_bond = [index, atom.index]
                            for o in o2:
                                this_bond.append(o)
                            bond = tuple(this_bond)
                            print bond
                            bond_dict[bond] = nbond
                nbond +=1
    
    print nbond
    for k,v in bond_dict.items():
        print k,v
    
    #Now we look for all the things with unit*factor bondlengths
    #Start with N (rectangles)
    for atom in Td:
        if atom.symbol == 'N':
            print'======================================='
            print 'N Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            model_molecule[n_obj][0].original_index = atom.index
            model_molecule[n_obj][0].symbol = 'F'
            indices, offsets = nlist.get_neighbors(atom.index)
            #Just checking the symbols here
            print indices
            print offsets
            symbols = ([Td[index].symbol for index in indices])
            symbol_string = ''.join(sorted([Td[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                print atom.index, index, offset
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                print bond_dict[bond]
                if Td[index].symbol == 'C':
                    #By definition, we can't be going over a periodic boundary (no pbc)
                    model_molecule[n_obj] += Td[index]
                    model_molecule[n_obj][-1].tag = bond_dict[bond]
                    model_molecule[n_obj].positions[-1] = Td.positions[index] 
    
    n_centers = n_obj
    
    ##Now we do the  C  (triangles)
    for atom in Td:
        if atom.symbol == 'C':
            #print'======================================='
            #print 'N Atom ',atom.index, " finding squares1"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                print index, offset, Td[index].symbol
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                #if not bond_dict.has_key(bond):
                    #Then
                if Td[index].symbol == 'N':
                    #By definition, we can't be going over a periodic boundary (no pbc)
                    model_molecule[n_obj] += Td[index]
                    model_molecule[n_obj][-1].tag = bond_dict[bond]
                    model_molecule[n_obj].positions[-1] = Td.positions[index] 
    
    
    f = open('Td.model','w')
    g = open('control-mofgen.txt','w')
    #Just for checking, now lets gather everything in the model_molecule dict into one big thing and print it
    test_mol = Atoms()
    for obj in model_molecule:
        test_mol += model_molecule[obj]
    
    write('test_model.xyz',test_mol)
    
    #print n_centers, n_model, n_obj
    #Headers and cell
    f.write('%-20s %-3d\n' %('Number of objects =',n_obj+1))
    f.write('%-20s\n' %('build = systre'))
    f.write('%5s\n' %('Cell:'))
    f.write('%8.3f %8.3f %8.3f \n' % (0.0, 0.0, 0.0))
    f.write('%8.3f %8.3f %8.3f \n' % (0.0, 0.0, 0.0))
    f.write('%8.3f %8.3f %8.3f \n' % (0.0, 0.0, 0.0))
    
    g.write('%-20s\n' %('model = Td'))
    
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Center: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = triangle'))
        #process positions to make it a bit more ideal
        for atom in model_molecule[obj]:
            #if atom.symbol == 'C':
            if atom.index != 0:
                model_molecule[obj].set_distance(atom.index,0,1.0,fix=1)
        for atom in model_molecule[obj]:
            (x,y,z) = atom.position
            #print atom.symbol, atom.position, atom.tag
            if atom.tag:
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
                #f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % (atom.symbol, x, y, z, atom.tag))
            else:
                f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Q', x, y, z))
                #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
        g.write('%-9s %-50s\n' %('center =', names[0]))
    
    for obj in xrange(n_centers+1,n_obj+1):
        f.write('\n%-8s %-3d\n' %('Linker: ', obj-n_centers))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = triangle'))
        #process positions to make it a bit more ideal
        for atom in model_molecule[obj]:
            #if atom.symbol == 'C':
            if atom.index != 0:
                model_molecule[obj].set_distance(atom.index,0,1.0,fix=1)
        for atom in model_molecule[obj]:
            (x,y,z) = atom.position
            #print atom.symbol, atom.position, atom.tag
            if atom.tag:
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
                #f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % (atom.symbol, x, y, z, atom.tag))
            else:
                f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Q', x, y, z))
                #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
        g.write('%-9s %-50s\n' %('linker =', names[1]))
    
    
    test_mol = Atoms()
    for obj in model_molecule:
        test_mol += model_molecule[obj]
    
    write('test_model2.xyz',test_mol)
    
    
