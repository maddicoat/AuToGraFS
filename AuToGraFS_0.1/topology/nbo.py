from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from ase.lattice.spacegroup import crystal
from ase.calculators.neighborlist import *
from ase.data import *
from math import pi
import numpy as np

#Note for when this becomes a subroutine
#We should multiply the factor by 2 seeing as we're halving everything
#factor = 6.0
def make_nbo(names,sizes):

    #work out the expansion factor
    #dist0 = furthest_dummy(mol0)
    #dist1 = furthest_dummy(mol1)

    factor = sum(sizes) #dist0 + dist1

    factor = factor * 2

    a = 2.0 * factor
    #a = 4.0
    
    # C -> squares, N => midpoints
    nbo = crystal(['C', 'N'], [(0.0, 0.5, 0.5), (0.25,0.0,0.5)], spacegroup=229,
                   cellpar=[a, a, a, 90, 90, 90])
    
    
    #write('test.xyz',nbo)
    
    eps = 0.05
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    n_tags = 0
    
    #First detect global bonding
    cov_rad=[]
    for atom in nbo:
        #cov_rad.append(covalent_radii[atom.number])
        cov_rad.append(factor / 4)
    
    nlist = NeighborList(cov_rad,self_interaction=False,bothways=True)
    nlist.build(nbo)
    
    #To sort out tags, we need to label each bond
    nbond = 1
    bond_matrix = np.zeros( (len(nbo),len(nbo)) )
    
    for atom in nbo:
        indices, offsets = nlist.get_neighbors(atom.index)
        for index in indices:
            if atom.index < index:
                bond_matrix[index,atom.index] = nbond
                bond_matrix[atom.index,index] = nbond
                print atom.index, index, nbond
                nbond +=1
    
    print nbond
    print bond_matrix 
    
    #Now we look for all the things with unit*factor bondlengths
    #Start with C (squares) 
    for atom in nbo:
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
            symbols = ([nbo[index].symbol for index in indices])
            symbol_string = ''.join(sorted([nbo[index].symbol for index in indices]))
            #print symbol_string
            #for i,o in zip(indices, offsets):
            #    print i,o
            for index,offset in zip(indices,offsets):
                dist_mic = nbo.get_distance(atom.index,index,mic=True)
                dist_no_mic = nbo.get_distance(atom.index,index,mic=False)
                print index, dist_no_mic
                ##Cell expands by the factor, but we're halving it because the N represesnt midpoints of CC bonds. Default systre bondlength is 1.0
                if (abs(dist_mic - dist_no_mic) <= eps) and (abs(dist_mic - factor/2.0) < eps): 
                     model_molecule[n_obj] += nbo[index]
                     model_molecule[n_obj][-1].tag = bond_matrix[atom.index,index]
                elif abs(dist_mic - factor/2.0) < eps:
                    #If we're going over a periodic boundary, we need to negate the tag
                    print "Tag, ", nbo[index].tag, " goes over pbc"
                    nbo[index].tag = -(nbo[index].tag)
                    model_molecule[n_obj] += nbo[index]
                    model_molecule[n_obj].positions[-1] = nbo.positions[index] + np.dot(offset, nbo.get_cell())
                    model_molecule[n_obj][-1].tag = -bond_matrix[atom.index,index]
    
    n_centers = n_obj
    #
    #
    ##Now we do the edges / N/  (linear things)
    for atom in nbo:
        if atom.symbol == 'N':
            print'======================================='
            print 'N Atom ',atom.index, " finding edges"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets):
                dist_mic = nbo.get_distance(atom.index,index,mic=True)
                dist_no_mic = nbo.get_distance(atom.index,index,mic=False)
                ##Cell expands by the factor, but we're halving it because the N represesnt midpoints of CC bonds. Default systre bondlength is 1.0
                if (abs(dist_mic - dist_no_mic) <= eps) and (abs(dist_mic - factor/2.0) < eps): 
                     model_molecule[n_obj] += nbo[index]
                     model_molecule[n_obj][-1].tag = bond_matrix[atom.index,index]
                elif abs(dist_mic - factor/2.0) < eps:
                    #If we're going over a periodic boundary, we need to negate the tag
                    print "Tag, ", nbo[index].tag, " goes over pbc"
                    nbo[index].tag = -(nbo[index].tag)
                    model_molecule[n_obj] += nbo[index]
                    model_molecule[n_obj].positions[-1] = nbo.positions[index] + np.dot(offset, nbo.get_cell())
                    model_molecule[n_obj][-1].tag = -bond_matrix[atom.index,index]
    
    
    
    f = open('nbo.model','w')
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
              (nbo.get_cell()[0][0],
               nbo.get_cell()[0][1],
               nbo.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (nbo.get_cell()[1][0],
               nbo.get_cell()[1][1],
               nbo.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (nbo.get_cell()[2][0],
               nbo.get_cell()[2][1],
               nbo.get_cell()[2][2])) 
   
    g.write('%-20s\n' %('model = nbo')) 
    
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Centre: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = square'))
        #process positions to make it a bit more ideal
        for atom in model_molecule[obj]:
            if atom.symbol == 'N': 
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
        f.write('%-20s\n' %('type = linear'))
        #process positions to make it a bit more ideal
        for atom in model_molecule[obj]:
            if atom.symbol == 'C':
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
