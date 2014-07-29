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

def make_sqp_B2(names,sizes):

    #ax110628-sr (ax110826.cif) personal communication from Prof. Dr. Jens Beckmann
    a = 10.9770
    b = 18.569
    c = 21.098
    alpha = 79.640
    beta =  76.060
    gamma = 85.670


    # P => central Si
    # Si => peripheral Si, actually using carbons from the original cif
    # N => centre of 4 Si/O
    # C => peripheral Si 

    sqp_B2 = crystal(['P', 'Si', 'Si', 'Si', 'Si', 'N', 'C', 'C', 'C', 'C'],
            [(0.65363, 0.82095, 0.14843), (0.4953, 0.8025, 0.2043), (0.7693, 0.7437, 0.1638), (0.7189, 0.9026, 0.1681), (0.6395, 0.8366, 0.0605),
             (0.240400, 0.263825, 0.376473), (-0.26884, 0.7000, 0.46936), (1.37428, 0.43072, 0.28200), (0.99881, 1.28479, 0.29492), (0.66767, 0.91355, -0.37629)],
            spacegroup=2,  #P-1
            cellpar=[a, b, c, alpha, beta, gamma])

    write('test.xyz',sqp_B2)
    write('test.cif',sqp_B2)

    eps = 0.10
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0

    # First detect global bonding
    # N-C distance is 3.5A
    cov_rad=[]
    for atom in sqp_B2:
        if atom.symbol == 'P':
            cov_rad.append(covalent_radii[atom.number])
        elif atom.symbol == 'Si':
            cov_rad.append(covalent_radii[atom.number])
        elif atom.symbol == 'N':
            cov_rad.append(2.0)
        elif atom.symbol == 'C':
            cov_rad.append(2.0)

    nlist = NeighborList(cov_rad,skin=0.01,self_interaction=False,bothways=True)
    nlist.build(sqp_B2)

    #To sort out tags, we need to label each bond, we get a bunch of extra bonds here that we need to filter out
    bond_matrix = np.zeros( (len(sqp_B2),len(sqp_B2)) , dtype=bool)

    for atom in sqp_B2:
        indices, offsets = nlist.get_neighbors(atom.index)
        if atom.symbol == 'P': #Here we go, suck in all Si
            for index,offset in zip(indices,offsets):
                if sqp_B2[index].symbol == 'Si':# and sqp_B2.get_distance(atom.index,index,mic=True) >10.0:
                    print "Atom ",atom.index, atom.symbol, " bonded to atom ", index, sqp_B2[index].symbol, " distance = ", sqp_B2.get_distance(atom.index,index,mic=True), "offset = ", offset
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True
        if atom.symbol == 'N': #Here we go, suck in all Si
            for index,offset in zip(indices,offsets):
                if sqp_B2[index].symbol == 'C': #and sqp_B2.get_distance(atom.index,index,mic=True) >10.0:
                    print "Atom ",atom.index, atom.symbol, " bonded to atom ", index, sqp_B2[index].symbol, " distance = ", sqp_B2.get_distance(atom.index,index,mic=True), "offset = ", offset
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True

    #Now we have a dirty trick, expand all the P-Si4 tetrahedra and find Si-C bonds
    for atom in sqp_B2:
        if atom.symbol == 'P':
            print'======================================='
            print 'P Atom ',atom.index
            indices, offsets = nlist.get_neighbors(atom.index)
            #print offsets
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                sqp_B2.set_distance(index,atom.index,10.0,fix=1)

    write('test2.xyz',sqp_B2)
    write('test2.cif',sqp_B2)
    nlist_CSi = NeighborList(cov_rad,skin=0.01,self_interaction=False,bothways=True)
    nlist_CSi.build(sqp_B2)

    print "Here"
    nbond = 1
    for atom in sqp_B2:
        indices, offsets = nlist_CSi.get_neighbors(atom.index)
        if atom.symbol == 'C': #Here we go, suck in all Si
            for index,offset in zip(indices,offsets):
                print "Atom ",atom.index, atom.symbol, " bonded to atom ", index, sqp_B2[index].symbol, " distance = ", sqp_B2.get_distance(atom.index,index,mic=True), "offset = ", offset
                if sqp_B2[index].symbol == 'Si':# and sqp_B2.get_distance(atom.index,index,mic=True) >10.0:
                    print "Adding bond"
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True
                    atom.tag = -nbond
                    sqp_B2[index].tag = -nbond
                    nbond += 1

    tags = sqp_B2.get_tags()
    print tags
    #Now build scaled model
    #base size 
    factor = sum(sizes) #(sizes[0] + sizes[2])/2.0
    print "factor = ",factor
    orig_a = a 
    scale = factor/orig_a
    b = b * scale
    c = c * scale

    # P => central Si
    # Si => peripheral Si, actually using carbons from the original cif
    # N => centre of 4 Si/O
    # C => peripheral Si 

    sqp_B2 = crystal(['P', 'Si', 'Si', 'Si', 'Si', 'N', 'C', 'C', 'C', 'C'],
            [(0.65363, 0.82095, 0.14843), (0.4953, 0.8025, 0.2043), (0.7693, 0.7437, 0.1638), (0.7189, 0.9026, 0.1681), (0.6395, 0.8366, 0.0605),
             (0.240400, 0.263825, 0.376473), (-0.26884, 0.7000, 0.46936), (1.37428, 0.43072, 0.28200), (0.99881, 1.28479, 0.29492), (0.66767, 0.91355, -0.37629)],
            spacegroup=2,  #P-1
            cellpar=[a, b, c, alpha, beta, gamma])

    sqp_B2.set_tags(tags)
        
    n_obj = -1 #re-initialise to -1 such that it can be indexed at start of add
    #Now we start assembling the model
    #We have only tetrahedra
    for atom in sqp_B2:
        if atom.symbol == 'P':
            print'======================================='
            print 'P Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            model_molecule[n_obj][0].original_index = atom.index
            model_molecule[n_obj][0].symbol = 'P'
            indices, offsets = nlist.get_neighbors(atom.index)
            #print offsets
            symbols = ([sqp_B2[index].symbol for index in indices])
            symbol_string = ''.join(sorted([sqp_B2[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                if bond_matrix[index,atom.index]:
                    model_molecule[n_obj] += sqp_B2[index]
                    model_molecule[n_obj].positions[-1] = sqp_B2.positions[index] + np.dot(offset, sqp_B2.get_cell())
                    model_molecule[n_obj][-1].tag = sqp_B2[index].tag
                    print sqp_B2[index].symbol, sqp_B2[index].position, sqp_B2[index].tag
                    print sqp_B2[index].index, sqp_B2[index].symbol, sqp_B2[index].position, sqp_B2[index].tag

    n_centers = n_obj

    # H-bond "square pyramids" are actually tetrahedra wrt centre of mass
    for atom in sqp_B2:
        if atom.symbol == 'N':
            print'======================================='
            print 'N Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            model_molecule[n_obj][0].original_index = atom.index
            model_molecule[n_obj][0].symbol = 'N'
            indices, offsets = nlist.get_neighbors(atom.index)
            #print offsets
            symbols = ([sqp_B2[index].symbol for index in indices])
            symbol_string = ''.join(sorted([sqp_B2[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                if bond_matrix[index,atom.index]:
                    model_molecule[n_obj] += sqp_B2[index]
                    model_molecule[n_obj].positions[-1] = sqp_B2.positions[index] + np.dot(offset, sqp_B2.get_cell())
                    model_molecule[n_obj][-1].tag = sqp_B2[index].tag
                    print sqp_B2[index].symbol, sqp_B2[index].position, sqp_B2[index].tag
                    print sqp_B2[index].index, sqp_B2[index].symbol, sqp_B2[index].position, sqp_B2[index].tag

    
    f = open('sqp_B2.model','w')
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
              (sqp_B2.get_cell()[0][0],
               sqp_B2.get_cell()[0][1],
               sqp_B2.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (sqp_B2.get_cell()[1][0],
               sqp_B2.get_cell()[1][1],
               sqp_B2.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (sqp_B2.get_cell()[2][0],
               sqp_B2.get_cell()[2][1],
               sqp_B2.get_cell()[2][2]))

    g.write('%-20s\n' %('model = sqp_B2'))

    #
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Center: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = tetrahedral'))
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
        f.write('%-20s\n' %('type = tetrahedral'))
        #process positions to make it a bit more ideal
        #this_dist = model_molecule[obj].get_distance(0,1)
        #for atom in model_molecule[obj]:
        #    #if atom.symbol == 'C':
        #    if atom.index != 0:
        #        model_molecule[obj].set_distance(0,atom.index,this_dist+1.0,fix=0)
        #for atom in model_molecule[obj]:
        #    #if atom.symbol == 'C':
        #    if atom.index != 0:
        #        model_molecule[obj].set_distance(0,atom.index,1.0,fix=1)
        for atom in model_molecule[obj]:
            (x,y,z) = atom.position
            print atom.symbol, atom.position, atom.tag
            if atom.tag:
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
                #f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % (atom.symbol, x, y, z, atom.tag))
            else:
                f.write('%-2s %15.8f %15.8f %15.8f\n' % ('Q', x, y, z))
                #f.write('%-2s %15.8f %15.8f %15.8f\n' % (atom.symbol, x, y, z))
        g.write('%-9s %-50s\n' %('linker =', names[1]))    
    #
    #
    test_mol = Atoms()
    for obj in model_molecule:
    #for obj in xrange(n_centers+1):#,n_obj+1): #triangles only
        test_mol += model_molecule[obj]
        test_mol.set_cell(sqp_B2.get_cell())
        test_mol.set_pbc(sqp_B2.get_pbc())

    write('test_model2.xyz',test_mol)
    write('test_model2.cif',test_mol)
