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

#So the cgd coordinates actually suck if you want ot build a framework...

    #work out the expansion factor
    #print "sizes =",sizes
    #factor =  sizes[2]#/2.0
    
    a = 28.7527 
    c = 35.4316  #uhm-7 crystal dimensions
    
    # Ne => centre of PW
    # C => corner of PW
    # N => triangles, first two corresponds to C1 linker, second one to Cs linker
    # Ca => corners of triangles (base end) i.e. adjacent to C_PW 
    # F => SiH2X2 locations, (COM)
    # H => SiH2H2 dummies 
    chs1 = crystal(['Ne', 'Ne', 'N', 'N', 'N', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'F', 'F', 'H', 'H', 'H'], 
                   [(0.451210, 0.820105, 0.5), (0.311875, 0.630620, 0.317373),
                    (0.5465, 0.0705, 0.8418), (0.3923, 0.1906, 0.8322), (0.22825, 0.05125, 0.344895),
                    (0.990100, 0.816000, 0.238990), (0.072300, 0.830500, 0.236450), (0.306260, 0.914300, 0.135060), (0.339440, 0.960790, 0.085750), (0.177400, 0.748400, 0.129220), (0.141300, 0.699100, 0.082930),
                    (0.429660, 0.821690, 0.288870), (0.648970, 0.890710, 0.5),
                    (0.093000, 0.684300, 0.457300), (0.115400, 0.430400, 0.304160), (0.185300, 0.478500, 0.254260)],
#   chs1 = crystal(['Ne', 'Ne', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'N', 'N', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'F', 'F', 'H', 'H', 'H'], 
#                  [(0.451210, 0.820105, 0.5), (0.311875, 0.630620, 0.317373),
#                   (0.893560, 0.736500, 0.153530), (0.991400, 0.666667, 0.052070), (0.845300, 0.885700, 0.215570), (0.795100, 0.815000, 0.146220), (0.692700, 0.891700, 0.050790), (0.944500, 0.813800, 0.219330),
#                   (0.5465, 0.0705, 0.8418), (0.3923, 0.1906, 0.8322), (0.22825, 0.05125, 0.344895),
#                   (0.990100, 0.816000, 0.238990), (0.072300, 0.830500, 0.236450), (0.306260, 0.914300, 0.135060), (0.339440, 0.960790, 0.085750), (0.177400, 0.748400, 0.129220), (0.141300, 0.699100, 0.082930),
#                   (0.429660, 0.821690, 0.288870), (0.648970, 0.890710, 0.5),
#                   (0.093000, 0.684300, 0.457300), (0.115400, 0.430400, 0.304160), (0.185300, 0.478500, 0.254260)],
                   spacegroup=87,  #I4/m
                   cellpar=[a, a, c, 90, 90, 90])
    
    
    write('test.xyz',chs1)
    write('test.cif',chs1)
    
    eps = 0.10
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0
    
    #First detect global bonding
    #Rigged such that each 
    cov_rad=[]
    for atom in chs1:
        #if atom.symbol == 'C':
        #    cov_rad.append(1.8)
        if atom.symbol == 'H':
            cov_rad.append(1.4)
        elif atom.symbol == 'N':
            cov_rad.append(3.7)
        elif atom.symbol == 'Ne':
            cov_rad.append(1.0)
        elif atom.symbol == 'Ca':
            cov_rad.append(3.0)
        elif atom.symbol == 'F':
            cov_rad.append(0.6)
        else: 
            cov_rad.append(covalent_radii[atom.number])
    
    print cov_rad
    
    nlist = NeighborList(cov_rad,skin=0.01,self_interaction=False,bothways=True)
    nlist.build(chs1)
    
    #To sort out tags, we need to label each bond, we get a bunch of extra bonds here that we need to filter out
    bond_matrix = np.zeros( (len(chs1),len(chs1)) , dtype=bool)

    for atom in chs1:
        indices, offsets = nlist.get_neighbors(atom.index)
        if atom.symbol == 'Ne': #or atom.symbol == 'F': #these two aren't affected by evil crossover business
            for index in indices:
                if chs1[index].symbol == 'Ca':
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True
        elif atom.symbol == 'F': #or atom.symbol == 'F': #these two aren't affected by evil crossover business
            for index in indices:
                if chs1[index].symbol == 'H':
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True
        elif atom.symbol == 'N': #Here we go, suck in both Ca
            for index in indices:
                if chs1[index].symbol == 'Ca' and chs1.get_distance(atom.index,index,mic=True) <5.0: #Fudgy, but rigging the covalent radii was a pain
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True
                elif chs1[index].symbol == 'H' and chs1.get_distance(atom.index,index,mic=True) >4.8: #This is a fudge that works. Better test would be take max H distance
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True
                else:
                    pass #No other bonds get counted


    ##Now we loop over the N, just to find which carbons go over PBC
    #for atom in chs1:
    #    if atom.symbol == 'N':
    #        indices, offsets = nlist.get_neighbors(atom.index)
    #        for index,offset in zip(indices,offsets): #Should only be C
    #            if chs1[index].symbol == 'C':
    #                if any(o != 0 for o in offset):
    #                    #If we're going over a periodic boundary, we need to negate the tag
    #                    chs1[index].tag = -chs1[index].tag
     
    factor = sum(sizes)/2.0 + 2#(sizes[0] + sizes[2])/2.0
    print "factor = ",factor
    a = 4.14300 * factor 
    c = 4.66569 * factor
    
    # Ne => centre of PW
    # N => triangles, first two corresponds to C1 linker, second one to Cs linker
    # Ca => corners of triangles (base end) i.e. adjacent to C_PW 
    # F => SiH2X2 locations, (COM)
    # H => SiH2H2 dummies 
    chs1 = crystal(['Ne', 'Ne', 'N', 'N', 'N', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'Ca', 'F', 'F', 'H', 'H', 'H'], 
                   [(0.451210, 0.820105, 0.5), (0.311875, 0.630620, 0.317373),
                    (0.5465, 0.0705, 0.8418), (0.3923, 0.1906, 0.8322), (0.22825, 0.05125, 0.344895),
                    (0.990100, 0.816000, 0.238990), (0.072300, 0.830500, 0.236450), (0.306260, 0.914300, 0.135060), (0.339440, 0.960790, 0.085750), (0.177400, 0.748400, 0.129220), (0.141300, 0.699100, 0.082930),
                    (0.429660, 0.821690, 0.288870), (0.648970, 0.890710, 0.5),
                    (0.093000, 0.684300, 0.457300), (0.115400, 0.430400, 0.304160), (0.185300, 0.478500, 0.254260)],
                   spacegroup=87,  #I4/m
                   cellpar=[a, a, c, 90, 90, 90])
   
    #Now that we've decided the bonds on the master structure, move them across to the scaled structure
    nbond = 1
    #Now our test will involve using indices,offsets and bond_matrix=True
    #In this case, the Ca and H atoms directly represent our dummy atoms, so we number them
    for atom in chs1:
        if atom.symbol == 'Ca' or atom.symbol == 'H':
            atom.tag = nbond
            nbond += 1

    #Now that we have our t/f matrix, we need to number the bonds, taking into account PBC
    bond_ids = np.zeros( (len(chs1),len(chs1)) )
    for atom in chs1:
        if atom.symbol == 'Ca' or atom.symbol == 'H':
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                if bond_matrix[index,atom.index]: #If we've already identified it a a bond 
                    if any(o != 0 for o in offset): #Then we're going over PBC
                        atom.tag = -atom.tag


    #for atom in chs1:
    #    print atom.index, atom.symbol, atom.position, atom.tag

    #Now we start assembling the model 
    #Start with Ne (rectangles) 
    for atom in chs1:
        if atom.symbol == 'Ne':
            print'======================================='
            print 'Ne Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            model_molecule[n_obj][0].original_index = atom.index
            model_molecule[n_obj][0].symbol = 'Ne' 
            indices, offsets = nlist.get_neighbors(atom.index)
            #print offsets
            symbols = ([chs1[index].symbol for index in indices])
            symbol_string = ''.join(sorted([chs1[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                if bond_matrix[index,atom.index]: 
                    model_molecule[n_obj] += chs1[index]
                    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    model_molecule[n_obj][-1].tag = chs1[index].tag 
                    print chs1[index].symbol, chs1[index].position, chs1[index].tag
                    print chs1[index].index, chs1[index].symbol, chs1[index].position, chs1[index].tag
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
                if bond_matrix[index,atom.index]: 
                    model_molecule[n_obj] += chs1[index]
                    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    model_molecule[n_obj][-1].tag = chs1[index].tag 
                    #if any(o != 0 for o in offset):
                    #    #If we're going over a periodic boundary, we need to negate the tag
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    #    model_molecule[n_obj][-1].tag = -bond_dict[bond]
                    #else:
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
   
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
                if bond_matrix[index,atom.index]: 
                    model_molecule[n_obj] += chs1[index]
                    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    model_molecule[n_obj][-1].tag = chs1[index].tag
                    #if any(o != 0 for o in offset):
                    #    #If we're going over a periodic boundary, we need to negate the tag
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
                    #    model_molecule[n_obj][-1].tag = -bond_dict[bond]
                    #else:
                    #    model_molecule[n_obj] += chs1[index]
                    #    model_molecule[n_obj][-1].tag = bond_dict[bond] 
                    #    model_molecule[n_obj].positions[-1] = chs1.positions[index] + np.dot(offset, chs1.get_cell())
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
        f.write('%-20s\n' %('type = linear_triangle'))
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
        for atom in model_molecule[obj]:
            #if atom.symbol == 'C': 
            if atom.index != 0:
                model_molecule[obj].set_distance(atom.index,0,2.0,fix=1)
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
