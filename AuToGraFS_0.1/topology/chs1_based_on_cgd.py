#!/usr/bin/python
from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from ase.lattice.spacegroup import crystal
from ase.calculators.neighborlist import *
from ase.data import *
from math import pi
import numpy as np

#Note for when this becomes a subroutine
#factor = 6.0

def make_chs1(names,sizes):

#CRYSTAL
#  NAME "UHM-7 - simplified - 12 - classiified as chs-1"
#  GROUP I4/mmm
#  CELL 4.14300 4.14300 4.66569 90.0000 90.0000 90.0000
#  NODE 5 4  0.00000 0.26972 0.00000
#  NODE 1 3  0.06146 0.28112 0.20686
#  NODE 4 4  0.10690 0.10690 0.21827
#  NODE 2 3  0.00000 0.19156 0.39361
#  EDGE  0.06146 0.28112 0.20686   0.21888 0.43854 0.29314
#  EDGE  0.00000 0.26972 0.00000   0.06146 0.28112 0.20686
#  EDGE  0.00000 0.19156 0.39361   0.00000 0.19156 0.60639
#  EDGE  0.00000 0.19156 0.39361   0.10690 0.10690 0.21827
#  EDGE  0.10690 0.10690 0.21827   0.28112 -0.06146 0.20686
## EDGE_CENTER  0.14017 0.35983 0.25000
## EDGE_CENTER  0.03073 0.27542 0.10343
## EDGE_CENTER  0.00000 0.19156 0.50000
## EDGE_CENTER  0.05345 0.14923 0.30594
## EDGE_CENTER  0.19401 0.02272 0.21257
#END

# So the node numbering above sucks...
# The edges go between (using the above numbering)
# 1-1
# 5-1
# 2-2
# 2-4
# 4-1
# The same-same edges between triangles are where we need to dump the SiH2X2 units
# Even though SiH2X2 is "linear" (two dummy atoms), in reality it is tetrahedral, and we need to orient it, three points on a line isn't enough
# We'll use some H atoms to give the direction

    #work out the expansion factor
    #print "sizes =",sizes
    factor = sum(sizes) #* 0.95#
    factor = sizes[0] + sizes[2] + 2.0#/2.0
    #factor =  sizes[2]#/2.0
    #factor = 3.0  #testing
    print "factor = ",factor
    
    a = 4.14300 * factor 
    c = 4.66569 * factor
    
    # C => rectangles (PW), second one is more like a half-hexagon
    # N => triangles, first one corresponds to C1 linker, second to Cs linker
    # F => SiH2X2 locations, X => other edge centres
    # H => SiH2H2 directions (COM goes here)
    chs1 = crystal(['C', 'C', 'N', 'N', 'F', 'X', 'F', 'X', 'X', 'H', 'H'], 
                   [(0.0, 0.26972, 0.0), (0.10690, 0.10690, 0.21827),
                    (0.06146, 0.28112, 0.20686), (0.0, 0.19156, 0.39361), 
                    (0.14017, 0.35983, 0.25), (0.03073, 0.27542, 0.10343), (0.0, 0.19156, 0.5), (0.05345, 0.14923, 0.30594), (0.19401, 0.02272, 0.21257),
                    (0.0, 0.22, 0.5), (0.315898, 0.81554, 0.25)],
                   spacegroup=139,  #I4/mmm
                   cellpar=[a, a, c, 90, 90, 90])
    
    
    write('test.xyz',chs1)
    write('test.cif',chs1)
    
    eps = 0.10
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0
    
    #First detect global bonding
    cov_rad=[]
    for atom in chs1:
        if atom.symbol == 'H':
            #cov_rad.append(covalent_radii[atom.number])
            cov_rad.append(factor / 8)
        else: #Normal
            cov_rad.append(factor / 4)
    
    print cov_rad
    
    nlist = NeighborList(cov_rad,skin=0.1,self_interaction=False,bothways=True)
    nlist.build(chs1)
    
    #To sort out tags, we need to label each bond
    nbond = 1
    bond_matrix = np.zeros( (len(chs1),len(chs1)) )
    pbc_bond_matrix = np.zeros( (len(chs1),len(chs1)) )
    
    bond_dict = {}
    
    for atom in chs1:
        print "---------------------"
        if atom.symbol == 'X':
            indices, offsets = nlist.get_neighbors(atom.index)
            #Both bonds get the same number 
            #print "indices are",indices
            #print "offsets are",offsets
            
            #if offsets.any():  #This test fails because some other dummies are within radius
            #    print "An offset is not zero!"
            pbc_flag = False
            for index,offset in zip(indices,offsets):
                if chs1[index].symbol == 'C' or chs1[index].symbol == 'N':
                    if offset.any():
                        pbc_flag = True

            for index,offset in zip(indices,offsets):
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                print bond
                if pbc_flag: 
                    bond_dict[bond] = -nbond
                else:
                    bond_dict[bond] = nbond
                print atom.index, index, atom.symbol, chs1[index].symbol, offset, chs1.get_distance(atom.index,index,mic=True), bond_dict[bond]
                ##Now we need to find the same bond the other direction to get the offsets right
                indices2,offsets2 = nlist.get_neighbors(index)
                for i2,o2 in zip(indices2,offsets2):
                    if i2 == atom.index:
                        #print "sum of offsets = ", offset, o2, sum(offset + o2)
                        if sum(offset + o2) == 0: #same bond
                            this_bond = [index, atom.index]
                            for o in o2:
                                this_bond.append(o)
                            bond = tuple(this_bond)
                            if pbc_flag: 
                                bond_dict[bond] = -nbond
                            else:
                                bond_dict[bond] = nbond
                #            print bond
                #            bond_dict[bond] = nbond
        nbond +=1
    
    for atom in chs1:
        print "---------------------"
        if atom.symbol == 'F':
            #both bonds get different numbers
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                print atom.index, index, atom.symbol, chs1[index].symbol, offset, chs1.get_distance(atom.index,index,mic=True)
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
    
    #for atom in chs1:
    #    print "---------------------"
    #    indices, offsets = nlist.get_neighbors(atom.index)
    #    for index,offset in zip(indices,offsets):
    #        print atom.index, index, atom.symbol, chs1[index].symbol, offset, chs1.get_distance(atom.index,index,mic=True)
    #        if atom.index < index:
    #            this_bond = [atom.index, index]
    #            for o in offset:
    #                this_bond.append(o)
    #            bond = tuple(this_bond)
    #            print bond
    #            bond_dict[bond] = nbond
    #            #Now we need to find the same bond the other direction to get the offsets right
    #            indices2,offsets2 = nlist.get_neighbors(index)
    #            for i2,o2 in zip(indices2,offsets2):
    #                if i2 == atom.index:
    #                    #print "sum of offsets = ", offset, o2, sum(offset + o2)
    #                    if sum(offset + o2) == 0: #same bond
    #                        this_bond = [index, atom.index]
    #                        for o in o2:
    #                            this_bond.append(o)
    #                        bond = tuple(this_bond)
    #                        print bond
    #                        bond_dict[bond] = nbond
    #            nbond +=1
    
    
    print "Bond dict:"
    print nbond
    for k,v in bond_dict.items():
        print k,v
    print "End Bond dict:"
    #Now we want to delete bonds that are not of length factor
    for k,v in bond_dict.items():
        i1,i2,o1,o2,o3 = k
        position1 = chs1.positions[i1]
        position2 = chs1.positions[i2]
        position2_mic = chs1.positions[i2] + np.dot([o1,o2,o3], chs1.get_cell())
        this_dist = np.linalg.norm(position1 - position2)
        this_dist_mic = np.linalg.norm(position1 - position2_mic)
        print i1,i2, this_dist, this_dist_mic, (this_dist - factor), this_dist_mic - factor
        #calculate differences
        if abs(this_dist - factor/2) > eps and all(o == 0 for o in [o1,o2,o3]):
            print "deleting", k
            del bond_dict[k]
        elif abs(this_dist_mic - factor/2) > eps and any(o != 0 for o in [o1,o2,o3]) :
            print "deleting", k
            del bond_dict[k]
        ##Now factors of 2    
        #elif abs(this_dist - factor/2) > eps and all(o == 0 for o in [o1,o2,o3]):
        #    print "deleting", k
        #    del bond_dict[k]
        #elif abs(this_dist_mic - factor/2) > eps and any(o != 0 for o in [o1,o2,o3]) :
        #    print "deleting", k
        #    del bond_dict[k]

    print "New Bond dict:"
    print len(bond_dict)
    for k,v in sorted(bond_dict.items()):
        i1,i2,o1,o2,o3 = k
        position1 = chs1.positions[i1]
        position2 = chs1.positions[i2]
        position2_mic = chs1.positions[i2] + np.dot([o1,o2,o3], chs1.get_cell())
        this_dist = np.linalg.norm(position1 - position2)
        this_dist_mic = np.linalg.norm(position1 - position2_mic)
        print k,v, chs1[i1].symbol, chs1[i2].symbol, this_dist,this_dist_mic
    print "End New Bond dict:"

    
    #Now we look for all the things with unit*factor bondlengths to F or X
    #Start with C (rectangles) 
    for atom in chs1:
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
            #print indices
            #print offsets
            symbols = ([chs1[index].symbol for index in indices])
            symbol_string = ''.join(sorted([chs1[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                #print atom.index, index, offset
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                #print bond_dict[bond]
                if bond in bond_dict:#chs1[index].symbol == 'X': #Should only be X because the F atoms were put between triangles
                    model_molecule[n_obj] += chs1[index]
                    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                    #if any(o != 0 for o in offset):
                    #    #If we're going over a periodic boundary, we need to negate the tag
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    #    model_molecule[n_obj][-1].tag = -bond_dict[bond]
                    #else:
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
    
    n_centers = n_obj
    #
    #
    ##Now we do the  F  (SiX2H2)
    for atom in chs1:
        if atom.symbol == 'F':
            #print'======================================='
            #print 'F Atom ',atom.index, " finding SiH2X2 anchor points"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                print index, offset, chs1[index].symbol
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                #if not bond_dict.has_key(bond):
                    #Then 
                if bond in bond_dict and chs1[index].symbol == 'N':
                    model_molecule[n_obj] += chs1[index]
                    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                    #if any(o != 0 for o in offset):
                    #    #If we're going over a periodic boundary, we need to negate the tag
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    #    model_molecule[n_obj][-1].tag = -bond_dict[bond]
                    #else:
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                if chs1[index].symbol == 'H': #direction marker, always in cell, just add it
                    model_molecule[n_obj] += chs1[index]
   
    n_linkers = n_obj

 
    ##Now we do the  N  (triangles)
    for atom in chs1:
        if atom.symbol == 'N':
            #print'======================================='
            #print 'N Atom ',atom.index, " finding triangles"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                print index, offset, chs1[index].symbol
                this_bond = [atom.index, index]
                for o in offset:
                    this_bond.append(o)
                bond = tuple(this_bond)
                #if not bond_dict.has_key(bond):
                    #Then 
                if bond in bond_dict:
                    if chs1[index].symbol == 'F' or chs1[index].symbol == 'X':
                        model_molecule[n_obj] += chs1[index]
                        model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                        model_molecule[n_obj][-1].tag = bond_dict[bond] 
                        #if any(o != 0 for o in offset):
                        #    #If we're going over a periodic boundary, we need to negate the tag
                        #    model_molecule[n_obj] += chs1[index]
                        #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                        #    model_molecule[n_obj][-1].tag = -bond_dict[bond]
                        #else:
                        #    model_molecule[n_obj] += chs1[index]
                        #    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                        #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                if chs1[index].symbol == 'H': #direction marker, always in cell, just add it
                    model_molecule[n_obj] += chs1[index]
    
   
    
  
    #
    #
    #
    f = open('chs1.model','w')
    g = open('control-mofgen.txt','w')
    #Just for checking, now lets gather everything in the model_molecule dict into one big thing and print it
    #test_mol = Atoms()
    #for obj in model_molecule:
    #    test_mol += model_molecule[obj]
    #
    #write('test_model.xyz',test_mol)
    
    #print n_centers, n_model, n_obj
    #Headers and cell
    f.write('%-20s %-3d\n' %('Number of objects =',n_obj+1))
    f.write('%-20s\n' %('build = systre'))
    f.write('%5s\n' %('Cell:'))
    f.write('%8.3f %8.3f %8.3f \n' %
              (chs1.get_cell()[0][0],
               chs1.get_cell()[0][1],
               chs1.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (chs1.get_cell()[1][0],
               chs1.get_cell()[1][1],
               chs1.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (chs1.get_cell()[2][0],
               chs1.get_cell()[2][1],
               chs1.get_cell()[2][2])) 

    g.write('%-20s\n' %('model = chs1'))
    
    # 
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Center: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = rectangle'))
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
    
    for obj in xrange(n_centers+1,n_linkers+1):
        f.write('\n%-8s %-3d\n' %('Linker: ', obj-n_centers))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = linear'))
        #process positions to make it a bit more ideal
        #for atom in model_molecule[obj]:
        #    #if atom.symbol == 'C': 
        #    if atom.index != 0:
        #        model_molecule[obj].set_distance(atom.index,0,1.0,fix=1)
        for atom in model_molecule[obj]:
            (x,y,z) = atom.position
            #print atom.symbol, atom.position, atom.tag
            if atom.tag:
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
                #f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % (atom.symbol, x, y, z, atom.tag))
            elif atom.symbol == 'F':
                f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Q', x, y, z))
                #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
            elif atom.symbol == 'H':
                f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Bq', x, y, z))
                #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
        g.write('%-9s %-50s\n' %('linker =', names[1]))
    #
    for obj in xrange(n_linkers+1,n_obj+1):
        f.write('\n%-8s %-3d\n' %('Linker: ', obj-n_centers))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = triangle'))
        #process positions to make it a bit more ideal
        #for atom in model_molecule[obj]:
        #    #if atom.symbol == 'C': 
        #    if atom.index != 0:
        #        model_molecule[obj].set_distance(atom.index,0,1.0,fix=1)
        for atom in model_molecule[obj]:
            (x,y,z) = atom.position
            #print atom.symbol, atom.position, atom.tag
            if atom.tag:
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
                #f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % (atom.symbol, x, y, z, atom.tag))
            else:
                f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Q', x, y, z))
                #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
        g.write('%-9s %-50s\n' %('linker =', names[2]))
    #
    #
    test_mol = Atoms()
    for obj in model_molecule:
    #for obj in xrange(n_linkers+1,n_obj+1): #triangles only 
        test_mol += model_molecule[obj]
        test_mol.set_cell(chs1.get_cell())
        test_mol.set_pbc(chs1.get_pbc())
    
    write('test_model2.xyz',test_mol)
    write('test_model2.cif',test_mol)
