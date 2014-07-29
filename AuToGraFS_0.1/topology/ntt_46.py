#!/usr/bin/python
from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from ase.lattice.spacegroup import crystal
from ase.calculators.neighborlist import *
from ase.data import *
from math import pi, sqrt, cos
import numpy as np

def make_ntt_46(names,sizes):
    print "This currently produces incorrect connectivity"
    #factor = sum(sizes) #dist0 + dist1
    #subsuming the two triangle vertices into one hexagonal vertex
    #if len(sizes) == 3:
    #    hex_size = sqrt(sizes[1]**2 + sizes[2]**2 - 2*sizes[1]*sizes[2]*cos(2.0/3.0*pi)) #2/3pi = 120 deg
    #else:
    #    hex_size = sqrt(sizes[1]**2 + sizes[1]**2 - 2*sizes[1]*sizes[1]*cos(2.0/3.0*pi)) #both triangles the same

    #print "Hex_size = ",hex_size
    #factor = sizes[0] + hex_size
    print sizes
    #factor = sizes[0] + sizes[1]/2.0 #unit cell is 44A, meant to be 37
    factor = sizes[0] + sizes[1]/3.0

    a = 7.3485 * factor 
    
    # C -> rectangle (square), N => central triangles, O => outer triangles
    ntt_46 = crystal(['C', 'N', 'O'], [(0.0, 0.16667, 0.16667), (0.1668, 0.1668, 0.3332), (0.1111, 0.1111, 0.2222 )],
                   spacegroup=225,  #Fm-3m
                   cellpar=[a, a, a, 90, 90, 90])
    
    
    #write('test.xyz',ntt_46)
    
    eps = 0.05
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0
    
    #First detect global bonding
    cov_rad=[]
    for atom in ntt_46:
        #cov_rad.append(covalent_radii[atom.number])
        cov_rad.append(factor / 2)
    
    print cov_rad[0]
    
    nlist = NeighborList(cov_rad,skin=0.1,self_interaction=False,bothways=True)
    nlist.build(ntt_46)
    
    #To sort out tags, we need to label each bond
    nbond = 1
    bond_dict = {}
    
    for atom in ntt_46:
        print "---------------------"
        indices, offsets = nlist.get_neighbors(atom.index)
        for index,offset in zip(indices,offsets):
            print atom.index, index, atom.symbol, ntt_46[index].symbol, offset, ntt_46.get_distance(atom.index,index,mic=True)
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
    
    #Start with C (rectangles) 
    for atom in ntt_46:
        if atom.symbol == 'C':
            print'======================================='
            print 'C Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            model_molecule[n_obj][0].original_index = atom.index
            model_molecule[n_obj][0].symbol = 'F' 
            indices, offsets = nlist.get_neighbors(atom.index)
            #Just checking the symbols here
            print indices
            print offsets
            symbols = ([ntt_46[index].symbol for index in indices])
            symbol_string = ''.join(sorted([ntt_46[index].symbol for index in indices]))
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
                if ntt_46[index].symbol == 'O': #C is only bound to O, (outer triangles)
                    if any(o != 0 for o in offset):
                        #If we're going over a periodic boundary, we need to negate the tag
                        model_molecule[n_obj] += ntt_46[index]
                        model_molecule[n_obj].positions[-1] = ntt_46.positions[index] + np.dot(offset, ntt_46.get_cell())
                        model_molecule[n_obj][-1].tag = -bond_dict[bond]
                    else:
                        model_molecule[n_obj] += ntt_46[index]
                        model_molecule[n_obj][-1].tag = bond_dict[bond] 
                        model_molecule[n_obj].positions[-1] = ntt_46.positions[index] + np.dot(offset, ntt_46.get_cell())
    
    n_centers = n_obj
    #
    #
    ##Now we do the  N  (inner triangles)
    #So, what we need to do is keep the positions of the N, but use the bond info of the attached O
    for atom in ntt_46:
        if atom.symbol == 'N':
            #print'======================================='
            #print 'N Atom ',atom.index, " finding edges"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                print "N check", index, offset, ntt_46[index].symbol
                #this_bond = [atom.index, index]
                #for o in offset:
                #    this_bond.append(o)
                #bond = tuple(this_bond)
                if ntt_46[index].symbol == 'O':
                    #Now we find the Carbons attached to each O
                    O_indices, O_offsets = nlist.get_neighbors(index) #this is the index of the O
                    for O_index,O_offset in zip(O_indices, O_offsets):
                        if ntt_46[O_index].symbol == 'C':
                            #Now bond info
                            this_bond = [index, O_index]
                            for o in O_offset:
                                this_bond.append(o)
                            bond = tuple(this_bond)
                            if any(o != 0 for o in O_offset):
                                #If we're going over a periodic boundary, we need to negate the tag
                                model_molecule[n_obj] += ntt_46[O_index]
                                model_molecule[n_obj].positions[-1] = ntt_46.positions[O_index] + np.dot(O_offset, ntt_46.get_cell())
                                model_molecule[n_obj][-1].tag = -bond_dict[bond]
                            else:
                                model_molecule[n_obj] += ntt_46[O_index]
                                model_molecule[n_obj][-1].tag = bond_dict[bond] 
                                model_molecule[n_obj].positions[-1] = ntt_46.positions[O_index] + np.dot(O_offset, ntt_46.get_cell())
                    #Now move each pair of C towards each other
                    this_dist = model_molecule[n_obj].get_distance(-1, -2)
                    model_molecule[n_obj].set_distance(-1, -2, this_dist - 1.0) #fix = 0.5 implied

    
    n_hex = n_obj 
    #
    ###Now we do the  O  (outer triangles)
    #for atom in ntt_46:
    #    if atom.symbol == 'O':
    #        #print'======================================='
    #        #print 'N Atom ',atom.index, " finding edges"
    #        n_obj+=1
    #        model_molecule[n_obj] = Atoms()
    #        model_molecule[n_obj] += atom
    #        indices, offsets = nlist.get_neighbors(atom.index)
    #        for index,offset in zip(indices,offsets):
    #            print index, offset, ntt_46[index].symbol
    #            this_bond = [atom.index, index]
    #            for o in offset:
    #                this_bond.append(o)
    #            bond = tuple(this_bond)
    #            if ntt_46[index].symbol == 'N' or ntt_46[index].symbol =='C':
    #                if any(o != 0 for o in offset):
    #                    #If we're going over a periodic boundary, we need to negate the tag
    #                    model_molecule[n_obj] += ntt_46[index]
    #                    model_molecule[n_obj].positions[-1] = ntt_46.positions[index] + np.dot(offset, ntt_46.get_cell())
    #                    model_molecule[n_obj][-1].tag = -bond_dict[bond]
    #                else:
    #                    model_molecule[n_obj] += ntt_46[index]
    #                    model_molecule[n_obj][-1].tag = bond_dict[bond] 
    #                    model_molecule[n_obj].positions[-1] = ntt_46.positions[index] + np.dot(offset, ntt_46.get_cell())
    #
    #
    #
    f = open('ntt_46.model','w')
    g = open('control-mofgen.txt','w')
    #Just for checking, now lets gather everything in the model_molecule dict into one big thing and print it
    #test_mol = Atoms()
    #for obj in model_molecule:
    #    test_mol += model_molecule[obj]
    
    #write('test_model.xyz',test_mol)
    
    #print n_centers, n_model, n_obj
    #Headers and cell
    f.write('%-20s %-3d\n' %('Number of objects =',n_obj+1))
    f.write('%-20s\n' %('build = systre'))
    f.write('%5s\n' %('Cell:'))
    f.write('%8.3f %8.3f %8.3f \n' %
              (ntt_46.get_cell()[0][0],
               ntt_46.get_cell()[0][1],
               ntt_46.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (ntt_46.get_cell()[1][0],
               ntt_46.get_cell()[1][1],
               ntt_46.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (ntt_46.get_cell()[2][0],
               ntt_46.get_cell()[2][1],
               ntt_46.get_cell()[2][2])) 
    
    g.write('%-20s\n' %('model = ntt_46'))

    #Now write stuff          
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Center: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = rectangle')) #Yep it's a rectangle
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
    
    for obj in xrange(n_centers+1,n_hex+1):
        f.write('\n%-8s %-3d\n' %('Linker: ', obj-n_centers))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = hexagon'))
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
    
    #for obj in xrange(n_inner_triangle+1,n_obj+1):
    #    f.write('\n%-8s %-3d\n' %('Linker: ', obj-n_centers))
    #    f.write('%-3d\n' %(len(model_molecule[obj])))
    #    f.write('%-20s\n' %('type = triangle'))
    #    #process positions to make it a bit more ideal
    #    for atom in model_molecule[obj]:
    #        #if atom.symbol == 'C': 
    #        if atom.index != 0:
    #            model_molecule[obj].set_distance(atom.index,0,1.0,fix=1)
    #    for atom in model_molecule[obj]:
    #        (x,y,z) = atom.position
    #        #print atom.symbol, atom.position, atom.tag
    #        if atom.tag:
    #            f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
    #            #f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % (atom.symbol, x, y, z, atom.tag))
    #        else:
    #            f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Q', x, y, z))
    #            #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
    #    if(len(names)) == 3:
    #        g.write('%-9s %-50s\n' %('linker =', names[2]))
    #    else:
    #        g.write('%-9s %-50s\n' %('linker =', names[1]))
    
    #
    #
    test_mol = Atoms()
    for obj in model_molecule:
        test_mol += model_molecule[obj]
    
    write('test_model2.xyz',test_mol)
