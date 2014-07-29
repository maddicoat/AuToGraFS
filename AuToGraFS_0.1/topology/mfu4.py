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

def make_mfu4(names,sizes):

    #work out the expansion factor
    factor = sum(sizes)


    a = 20.0    #We'll form the bondlist and then enlarge the cell.
                #At a = 20Ang, Zn-N is 5Ang
    
    # Zn -> centre of mfu4, N -> centre of rectangle linker , C -> dummies 
    mfu4 = crystal(['Zn', 'N', 'C'], [(0.25, 0.25, 0.25), (0.0, 0.25, 0.25), (0.0539, 0.7248, 0.7752 )], spacegroup=225, #Fm-3m
                   cellpar=[a, a, a, 90, 90, 90])
    
    
    #write('test.xyz',mfu4)
    
    eps = 0.05
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    
    #First detect global bonding
    cov_rad=[]
    for atom in mfu4:
        if atom.symbol == 'Zn':
            cov_rad.append(covalent_radii[atom.number]*3) #Rigged to suck in C, but not N
        else:
            cov_rad.append(covalent_radii[atom.number])
    
    #print cov_rad
    
    nlist = NeighborList(cov_rad,skin=0.1,self_interaction=False,bothways=True)
    nlist.build(mfu4)
    
    #Now we exapnd the cell. The factor given is 1/4 the required unit cell
    
    a = factor * 4
    # Zn -> centre of mfu4, N -> centre of rectangle linker , C -> dummies 
    mfu4 = crystal(['Zn', 'N', 'C'], [(0.25, 0.25, 0.25), (0.0, 0.25, 0.25), (0.0539, 0.7248, 0.7752 )], spacegroup=225, #Fm-3m
                   cellpar=[a, a, a, 90, 90, 90])
    
    
    #To sort out tags, we need to label each bond
    nbond = 1
    #Now in this case, the C atoms directly represent our dummy atoms, so we number them
    for atom in mfu4:
        if atom.symbol == 'C':
            atom.tag = nbond
            nbond += 1
    
    
    #Now we loop over the N, just to find which carbons go over PBC 
    for atom in mfu4:
        if atom.symbol == 'N':
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets): #Should only be C
                if mfu4[index].symbol == 'C':
                    if any(o != 0 for o in offset):
                        #If we're going over a periodic boundary, we need to negate the tag
                        mfu4[index].tag = -mfu4[index].tag
    
    #print mfu4.get_tags()
    
    #Now we build the model  
    #Start with Zn
    for atom in mfu4:
        if atom.symbol == 'Zn':
            #print'======================================='
            #print 'Zn Atom ',atom.index
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            model_molecule[n_obj][0].original_index = atom.index
            indices, offsets = nlist.get_neighbors(atom.index)
            #Just checking the symbols here
            #print indices
            #print offsets
            #symbols = ([mfu4[index].symbol for index in indices])
            #symbol_string = ''.join(sorted([mfu4[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets): #Should only be C
                if mfu4[index].symbol == 'C':
                    model_molecule[n_obj] += mfu4[index]
                    model_molecule[n_obj][-1].tag = mfu4[index].tag 
                    model_molecule[n_obj].positions[-1] = mfu4.positions[index] + np.dot(offset, mfu4.get_cell())
    
    n_centers = n_obj
    #
    #
    ##Now we do the  N  (rectangles)
    for atom in mfu4:
        if atom.symbol == 'N':
            #print'======================================='
            #print 'N Atom ',atom.index, " finding squares1"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets): #Should only be C
                if mfu4[index].symbol == 'C':
                    model_molecule[n_obj] += mfu4[index]
                    model_molecule[n_obj][-1].tag = mfu4[index].tag 
                    model_molecule[n_obj].positions[-1] = mfu4.positions[index] + np.dot(offset, mfu4.get_cell())
    #
    #
    #
    f = open('mfu4.model','w')
    g = open('control-mofgen.txt','w')
    #Just for checking, now lets gather everything in the model_molecule dict into one big thing and print it
    test_mol = Atoms()
    for obj in model_molecule:
        test_mol += model_molecule[obj]
    
    #write('test_model.xyz',test_mol)
    
    #print n_centers, n_model, n_obj
    #Headers and cell
    f.write('%-20s %-3d\n' %('Number of objects =',n_obj+1))
    f.write('%-20s\n' %('build = systre'))
    f.write('%5s\n' %('Cell:'))
    f.write('%8.3f %8.3f %8.3f \n' %
              (mfu4.get_cell()[0][0],
               mfu4.get_cell()[0][1],
               mfu4.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (mfu4.get_cell()[1][0],
               mfu4.get_cell()[1][1],
               mfu4.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (mfu4.get_cell()[2][0],
               mfu4.get_cell()[2][1],
               mfu4.get_cell()[2][2])) 

    g.write('%-20s\n' %('model = mfu4')) 
    # Print the centers and linkers
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Center: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = mfu4'))
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
        g.write('%-9s %-50s\n' %('center =', names[0]))
    
    for obj in xrange(n_centers+1,n_obj+1):
        f.write('\n%-8s %-3d\n' %('Linker: ', obj-n_centers))
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
        g.write('%-9s %-50s\n' %('linker =', names[1]))
    
    #
    #
    #test_mol = Atoms()
    #for obj in model_molecule:
    #    test_mol += model_molecule[obj]
    #
    #write('test_model2.xyz',test_mol)
