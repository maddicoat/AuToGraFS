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

def make_sqp_B1(names,sizes):

    #ax110628-sr (ax110826.cif) personal communication from Prof. Dr. Jens Beckmann
    a = 10.9770
    b = 18.569
    c = 21.098
    alpha = 79.640
    beta =  76.060
    gamma = 85.670


    # P => central Si
    # Si => peripheral Si, actually using carbons from the original cif

    sqp_b = crystal(['P', 'Si', 'Si', 'Si', 'Si'],
            [(0.65363, 0.82095, 0.14843), (0.4953, 0.8025, 0.2043), (0.7693, 0.7437, 0.1638), (0.7189, 0.9026, 0.1681), (0.6395, 0.8366, 0.0605)],
            spacegroup=2,  #P-1
            cellpar=[a, b, c, alpha, beta, gamma])

    write('test.xyz',sqp_b)
    write('test.cif',sqp_b)

    eps = 0.10
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0

    #First detect global bonding
    #The P-Si distance we need is 10.859 Ang
    #The P-Si distance from the other layer is 9.7503 Ang

    cov_rad=[]
    for atom in sqp_b:
        if atom.symbol == 'P':
            cov_rad.append(covalent_radii[atom.number])
        elif atom.symbol == 'Si':
            cov_rad.append(covalent_radii[atom.number])

    nlist = NeighborList(cov_rad,skin=0.01,self_interaction=False,bothways=True)
    nlist.build(sqp_b)

    #To sort out tags, we need to label each bond, we get a bunch of extra bonds here that we need to filter out
    bond_matrix = np.zeros( (len(sqp_b),len(sqp_b)) , dtype=bool)

    for atom in sqp_b:
        indices, offsets = nlist.get_neighbors(atom.index)
        if atom.symbol == 'P': #Here we go, suck in all Si
            for index,offset in zip(indices,offsets):
                if sqp_b[index].symbol == 'Si': #and sqp_b.get_distance(atom.index,index,mic=True) >10.0:
                    print "Atom ",atom.index, atom.symbol, " bonded to atom ", index, sqp_b[index].symbol, " distance = ", sqp_b.get_distance(atom.index,index,mic=True), "offset = ", offset
                    bond_matrix[index,atom.index] = True
                    bond_matrix[atom.index,index] = True


    #Now build scaled model
    #base size 
    factor = sizes[0] #(sizes[0] + sizes[2])/2.0
    print "factor = ",factor
    orig_a = a 
    scale = factor/orig_a
    b = b * scale
    c = c * scale

    # P => central Si
    # Si => peripheral Si

    sqp_b = crystal(['P', 'Si', 'Si', 'Si', 'Si'],
            [(0.65363, 0.82095, 0.14843), (0.4953, 0.8025, 0.2043), (0.7693, 0.7437, 0.1638), (0.7189, 0.9026, 0.1681), (0.6395, 0.8366, 0.0605)],
            spacegroup=2,  #P-1
            cellpar=[factor, b, c, alpha, beta, gamma])

    #All bonding info is incorrect for now, seeing as the structure has no bonding

    #Now that we've decided the bonds on the master structure, move them across to the scaled structure
    nbond = 1
    #Now our test will involve using indices,offsets and bond_matrix=True
    #In this case, the Si atoms directly represent our dummy atoms, so we number them
    for atom in sqp_b:
        if atom.symbol == 'Si':
            atom.tag = nbond
            nbond += 1

    #Now that we have our t/f matrix, we need to number the bonds, taking into account PBC
    bond_ids = np.zeros( (len(sqp_b),len(sqp_b)) )
    for atom in sqp_b:
        if atom.symbol == 'Si':
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                if bond_matrix[index,atom.index]: #If we've already identified it a a bond
                    if any(o != 0 for o in offset): #Then we're going over PBC
                        atom.tag = -atom.tag

    #Now we start assembling the model
    #We have only tetrahedra
    for atom in sqp_b:
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
            symbols = ([sqp_b[index].symbol for index in indices])
            symbol_string = ''.join(sorted([sqp_b[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                if bond_matrix[index,atom.index]:
                    model_molecule[n_obj] += sqp_b[index]
                    model_molecule[n_obj].positions[-1] = sqp_b.positions[index] + np.dot(offset, sqp_b.get_cell())
                    model_molecule[n_obj][-1].tag = sqp_b[index].tag
                    print sqp_b[index].symbol, sqp_b[index].position, sqp_b[index].tag
                    print sqp_b[index].index, sqp_b[index].symbol, sqp_b[index].position, sqp_b[index].tag

    n_centers = n_obj

    for atom in sqp_b:
        if atom.symbol == 'Si':
            print'======================================='
            print 'Si Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            #indices, offsets = nlist.get_neighbors(atom.index)
            #for index,offset in zip(indices,offsets):
            #    if bond_matrix[index,atom.index]:
            #        model_molecule[n_obj] += sqp_b[index]
            #        model_molecule[n_obj].positions[-1] = sqp_b.positions[index] + np.dot(offset, sqp_b.get_cell())
            #        model_molecule[n_obj][-1].tag = sqp_b[index].tag
            #        print sqp_b[index].symbol, sqp_b[index].position, sqp_b[index].tag
            #        print sqp_b[index].index, sqp_b[index].symbol, sqp_b[index].position, sqp_b[index].tag
            model_molecule[n_obj] += atom
            model_molecule[n_obj][-1].original_index = atom.index
            model_molecule[n_obj][-1].symbol = 'Si'
            
    
    f = open('sqp_B1.model','w')
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
              (sqp_b.get_cell()[0][0],
               sqp_b.get_cell()[0][1],
               sqp_b.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (sqp_b.get_cell()[1][0],
               sqp_b.get_cell()[1][1],
               sqp_b.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (sqp_b.get_cell()[2][0],
               sqp_b.get_cell()[2][1],
               sqp_b.get_cell()[2][2]))

    g.write('%-20s\n' %('model = sqp_B1'))

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
        f.write('%-20s\n' %('type = ZeroD_cap'))
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
        test_mol.set_cell(sqp_b.get_cell())
        test_mol.set_pbc(sqp_b.get_pbc())

    write('test_model2.xyz',test_mol)
    write('test_model2.cif',test_mol)
