from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from ase.lattice.spacegroup import crystal
from ase.calculators.neighborlist import *
from ase.data import *
from math import pi, sin, cos
import numpy as np

#Note for when this becomes a subroutine
#factor = 6.0

def make_mil53(names,sizes):

    #work out the expansion factor
    factor = sizes[1]   #linker is the only thing that matters #sum(sizes)

    #print sizes[0]
    a = 6.6085    #We'll form the bondlist and then enlarge the cell.
    b = 16.0      #16.675           
    c = 12.0      #12.813          
    #These cell dimensions are rigged to a linker length of 5 (circa size of benzene)
    
    # Al -> centre of mil53, N -> centre of rectangle linker , C -> dummies 
    mil53 = crystal(['Al', 'N', 'C'], [(0.25, 0.25, 0.25), (0.5, 0.0, 0.5), (0.5, 0.1292, 0.383 )], spacegroup=74, #Imma
                   cellpar=[a, b, c, 90, 90, 90])
    
    
    write('test.xyz',mil53)
    
    eps = 0.05
    model_molecule={}
    n_obj = -1 #initialise to -1 such that it can be indexed at start of add
    
    #First detect global bonding
    cov_rad=[]
    for atom in mil53:
        if atom.symbol == 'Al':
            cov_rad.append(covalent_radii[atom.number]*2) #Rigged to suck in C, but not N
        elif atom.symbol == 'N':
            cov_rad.append(covalent_radii[atom.number]*3) #Rigged to suck in C, but not Al
        else:
            cov_rad.append(covalent_radii[atom.number])
    
    #print cov_rad
    
    nlist = NeighborList(cov_rad,skin=0.1,self_interaction=False,bothways=True)
    nlist.build(mil53)
    
    #Now we exapnd the cell. The factor given is half the hypotenuse of a right triangle with angle 35 deg

    print "factor = " , factor 
    print b,c
    a = 2* sizes[0]
    b = b + (2*factor * sin(35*pi/180.0))
    c = c + (2*factor * cos(35*pi/180.0))
    print b,c
    # Al -> centre of mil53, N -> centre of rectangle linker , C -> dummies 
    mil53 = crystal(['Al', 'N', 'C'], [(0.25, 0.25, 0.25), (0.5, 0.0, 0.5), (0.5, 0.1292, 0.383 )], spacegroup=74, #Imma
                   cellpar=[a, b, c, 90, 90, 90]) #a (height) is fixed
    
    
    #To sort out tags, we need to label each bond
    nbond = 1
    #Now in this case, the C atoms directly represent our dummy atoms, so we number them
    for atom in mil53:
        if atom.symbol == 'C':
            atom.tag = nbond
            nbond += 1
    #Al also need to be tagged to each other
    for atom in mil53:
        if atom.symbol == 'Al':
            indices, offsets = nlist.get_neighbors(atom.index)
            #symbol_string = ''.join(sorted([mil53[index].symbol for index in indices]))
            #print symbol_string
            for index,offset in zip(indices,offsets):
                if mil53[index].symbol == 'Al' and all(o == 0 for o in offset) and index > atom.index: #interesting because it the same bond is both internal and pbc
                    #print "atom ", atom.index, " and atom ",index, " are bonded, non-pbc ", nbond
                    atom.tag = nbond
                    mil53[index].tag = nbond
                    nbond += 1

    
    
    #Now we loop over the N, just to find which carbons go over PBC 
    for atom in mil53:
        if atom.symbol == 'N':
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets): #Should only be C
                if mil53[index].symbol == 'C':
                    if any(o != 0 for o in offset):
                        #If we're going over a periodic boundary, we need to negate the tag
                        mil53[index].tag = -mil53[index].tag
    
    #print mil53.get_tags()
    
    used_carbon = []
    #Now we build the model  
    #Start with Al
    for atom in mil53:
        if atom.symbol == 'Al':
            #print'======================================='
            #print 'Al Atom ',atom.index
            indices, offsets = nlist.get_neighbors(atom.index)
            C_count = 0
            for index,offset in zip(indices,offsets): #Should only be C or Al
                if mil53[index].symbol == 'C'  and all(o == 0 for o in offset): # and mil53[index].tag not in  model_molecule[n_obj].get_tags():
                    C_count += 1
            #print "Al atom, ", atom.index, "has c-count in this cell", C_count
            if C_count == 2:
                n_obj+=1
                model_molecule[n_obj] = Atoms()
                model_molecule[n_obj] += atom
                for index,offset in zip(indices,offsets): #Should only be C or Al
                    if mil53[index].symbol == 'C'  and all(o == 0 for o in offset): # and mil53[index].tag not in  model_molecule[n_obj].get_tags():
                        model_molecule[n_obj] += mil53[index]
                        model_molecule[n_obj][-1].tag = mil53[index].tag
                        model_molecule[n_obj].positions[-1] = mil53.positions[index] + np.dot(offset, mil53.get_cell())
                        used_carbon.append(index)
                #We need to put the other Al as an anchor         
                for index,offset in zip(indices,offsets): 
                    if mil53[index].symbol =='Al' and all(o == 0 for o in offset):
                        model_molecule[n_obj] += mil53[index]
                        #model_molecule[n_obj][-1].tag = None #No tag
                        model_molecule[n_obj].positions[-1] = mil53.positions[index] + np.dot(offset, mil53.get_cell())

   #That gives us two Al centres, now the other two
   #These are the two that have 4 in-cell C
    for atom in mil53:
        if atom.symbol == 'Al':
            #print'======================================='
            #print 'Al Atom ',atom.index
            indices, offsets = nlist.get_neighbors(atom.index)
            C_count = 0
            for index,offset in zip(indices,offsets): #Should only be C or Al
                if mil53[index].symbol == 'C'  and all(o == 0 for o in offset): # and mil53[index].tag not in  model_molecule[n_obj].get_tags():
                    C_count += 1
            #print "Al atom, ", atom.index, "has c-count in this cell", C_count
            if C_count == 4:
                n_obj+=1
                model_molecule[n_obj] = Atoms()
                model_molecule[n_obj] += atom
                for index,offset in zip(indices,offsets): #Should only be C or Al
                    if mil53[index].symbol == 'C'  and all(o == 0 for o in offset) and index not in used_carbon:
                        model_molecule[n_obj] += mil53[index]
                        model_molecule[n_obj][-1].tag = mil53[index].tag
                        model_molecule[n_obj].positions[-1] = mil53.positions[index] + np.dot(offset, mil53.get_cell())
                #We need to put the other Al as an anchor         
                for index,offset in zip(indices,offsets): 
                    if mil53[index].symbol =='Al' and any(o != 0 for o in offset):
                        model_molecule[n_obj] += mil53[index]
                        #model_molecule[n_obj][-1].tag = None #No tag
                        model_molecule[n_obj].positions[-1] = mil53.positions[index] + np.dot(offset, mil53.get_cell())


#                    #We need to interrogate our C further as two of them belong to both Al in the same cell
#                    C_indices,C_offsets = nlist.get_neighbors(index)
#                    symbol_string = ''.join(sorted([mil53[index].symbol for index in C_indices]))
#                    print symbol_string
#
#                    model_molecule[n_obj] += mil53[index]
#                    model_molecule[n_obj][-1].tag = mil53[index].tag 
#                    model_molecule[n_obj].positions[-1] = mil53.positions[index] + np.dot(offset, mil53.get_cell())
    
    n_centers = n_obj
    #
    #
    ##Now we do the  N  (linear)
    for atom in mil53:
        if atom.symbol == 'N':
            #print'======================================='
            #print 'N Atom ',atom.index, " finding squares1"
            n_obj+=1
            model_molecule[n_obj] = Atoms()
            model_molecule[n_obj] += atom
            indices, offsets = nlist.get_neighbors(atom.index)
            for index,offset in zip(indices,offsets): #Should only be C
                if mil53[index].symbol == 'C':
                    model_molecule[n_obj] += mil53[index]
                    model_molecule[n_obj][-1].tag = mil53[index].tag 
                    model_molecule[n_obj].positions[-1] = mil53.positions[index] + np.dot(offset, mil53.get_cell())
    #
    #
    #
    f = open('mil53.model','w')
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
    f.write('%-13s\n' %('ExtraDummies'))
    f.write('%5s\n' %('Cell:'))
    f.write('%8.3f %8.3f %8.3f \n' %
              (mil53.get_cell()[0][0],
               mil53.get_cell()[0][1],
               mil53.get_cell()[0][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (mil53.get_cell()[1][0],
               mil53.get_cell()[1][1],
               mil53.get_cell()[1][2]))
    f.write('%8.3f %8.3f %8.3f \n' %
              (mil53.get_cell()[2][0],
               mil53.get_cell()[2][1],
               mil53.get_cell()[2][2])) 
    g.write('%-20s\n' %('model = mil53')) 
    # Print the centers and linkers
    for obj in xrange(n_centers+1):
        f.write('\n%-8s %-3d\n' %('Center: ', obj+1))
        f.write('%-3d\n' %(len(model_molecule[obj])))
        f.write('%-20s\n' %('type = mil53'))
        #process positions to make it a bit more ideal
        #for atom in model_molecule[obj]:
        #    #if atom.symbol == 'C':
        #    if atom.index != 0:
        #        model_molecule[obj].set_distance(atom.index,0,2.0,fix=1)
        for atom in model_molecule[obj]:
            (x,y,z) = atom.position
            #print atom.symbol, atom.position, atom.tag
            if atom.symbol == 'C': #can't use tags as a marker here, because we've given a tag to 'Q' as well
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('X', x, y, z, atom.tag))
            elif atom.symbol == 'Al':
                f.write('%-2s %15.8f %15.8f %15.8f %-4s\n' % ('Q', x, y, z, atom.tag))
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
    test_mol = Atoms()
    for obj in model_molecule:
        test_mol += model_molecule[obj]
    
    write('test_model2.xyz',test_mol)
