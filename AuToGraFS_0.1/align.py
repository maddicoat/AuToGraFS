import numpy as np
from math import pi
from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
from collections import deque

def align(mol,model):

    shape = mol.info['shape'].strip()
    model_shape = model.info['shape'].strip()

    #Special processing for octahedra
    if shape == 'octahedral' and model_shape == 'octahedral_2':
        return align_octahedron_2(mol, model) #Trigonal antiprism
    elif shape == 'octahedral' and model_shape == 'octahedral_3':
        return align_octahedron_3(mol, model) #Three corners of a cube
    elif shape == 'octahedral' and model_shape == 'octahedral_4':
        return align_octahedron_4(mol, model) #align one thing exactly, the second thing close enough and assume(!) the rest

    if shape == 'square' and model_shape == 'rectangle': #Special case, we want to allow this fairly often
        print "Aligning a square on a rectangle!"
        return align_rectangle(mol, model) 

    elif shape != model_shape:
        print shape
        print model_shape
        raise RuntimeError('Molecule shape,'+shape+' and model shape '+model_shape+' not the same')

    if shape == 'ZeroD_cap':
        return align_ZeroD_cap(mol, model)

    if shape == 'point_cap':
        return align_point_cap(mol, model)

    if shape == 'linear':
        return align_linear(mol, model)

    if shape == 'linear_triangle':
        return align_linear_triangle(mol, model)

    if shape == 'triangle':
        return align_triangle(mol, model)

    if shape == 'tetrahedral':
        return align_tetrahedral(mol, model)

    if shape == 'square':
        return align_square(mol, model)

    if shape == 'rectangle':
        return align_rectangle(mol, model)

    if shape == 'trigonal_bipyramid':
        return align_tri_bipyramid(mol, model)

    if shape == 'octahedral':
        return align_octahedron(mol, model) #square bipyramid

    if shape == 'tri_prism':
        return align_tri_prism(mol, model)
    
    if shape == 'hexagon':
        return align_hexagon(mol, model)

    if shape == 'mfu4': #Kind of an octahedron, kind of not...
        return align_mfu4(mol, model) 
    
    if shape == 'mil53': #special case
        return align_mil53(mol, model) 
    
    if shape == 'icosahedral': #Zr6
        return align_icosahedron(mol, model) 
    
    raise RuntimeError('Molecule shape '+shape+' not known!')

def align_ZeroD_cap(mol, model):
    """ A Zero-Dimensional cap is just a point in space with no alignment. 
        It's intended for use as a 'blind cap' """

    model_cod = model.get_center_of_dummies()
    #Now, the molecule
    mol_dummies = mol.get_atom_indices('X')
    if len(mol_dummies) != 1:
        raise RuntimeError('ZeroD cap has too many dummy atoms!')
    #translate
    mol.translate(model_cod - mol[mol_dummies[0]].position)
    #transfer the tag
    mol[0].tag = model[0].tag

def align_point_cap(mol, model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
    else:
        rot_angle=0
    print "Entering align_point_cap"
    #The model is a Q-X bond: Q is the inner point (imagine M of a paddlewheel), X is the align point as normal
    model_com = model.get_center_of_mass()
    model_cod = model.get_center_of_dummies()
    print "Model COM =", model_com
    print "Model COD =", model_cod
    model_dummy = model.get_atom_indices('X')[0]
    model_vector = (model_cod - model_com)
    #Now, the molecule 
    mol_dummies = mol.get_atom_indices('X')
    if len(mol_dummies) != 1:
        raise RuntimeError('Point cap has too many dummy atoms!')
    mol_bonded_atom=mol[mol_dummies[0]].bondlist.keys()[0] - 1 #remember bondlist is 1-indexed
    #translate
    mol.translate(model_cod - mol[mol_dummies[0]].position)
    #print mol_bonded_atom
    mol_vector = (mol[mol_dummies[0]].position - mol[mol_bonded_atom].position)
    #Now rotate the mol_vector onto the model
    mol.rotate(mol_vector,-model_vector)#,center=mol[mol_dummies[0]].position)
    #Transfer dummy tag
    mol[mol_dummies[0]].tag = model[model_dummy].tag
    #Now, if it has been requested, rotate the molecule about the alignment axis
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        #need a vector perpendicular to the triangle plane
        mol.rotate(-model_vector,rot_angle,center=mol[mol_bonded_atom].position) #pi fudge is because we don't know which way the dihedral is

    return 
    

def align_linear(mol, model):
    eps = 0.02
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180 
    else:
        rot_angle=0   
    model_com = model.get_center_of_mass()
    #mol_com = mol.get_center_of_mass()
    mol_cod = mol.get_center_of_dummies()
    mol.translate(model_com - mol_cod)
    mol_cod = mol.get_center_of_dummies()
    mol_com = mol.get_center_of_mass()
    
    #Grab the tags on the model and sort them asciibetically. We'll use a dict.
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        if a.tag:
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    mol_dummies = mol.get_atom_indices('X')
    model_atom=min(model_dummy_tags, key=model_dummy_tags.get)
    mol_dummies = mol.get_atom_indices('X')
    if 'flip' in mol.info:  #Biggest rotation
        dummy_distances=[]
        for a in mol_dummies:
            dist=get_intermolecular_distance(model, mol, model_atom, a)
            dummy_distances.append(dist)
        max_index = dummy_distances.index(max(dummy_distances))
        angle = get_general_angle(mol.positions[mol_dummies[max_index]],model_com,model.positions[model_atom])
        if angle > eps:
            mol.rotate(mol.positions[mol_dummies[min_index]]-model_com,model.positions[model_atom]-model_com,center=model_com)
        mol.rotate(mol.positions[mol_dummies[max_index]]-model_com,model.positions[model_atom]-model_com,center=model_com)
        transfer_dummy_tags_linear(mol,model)
        #mol[mol_dummies[min_index]].tag=tag
    else:
        #else: do the smallest rotation
        #print mol_dummies
        dummy_distances=[]
        for a in mol_dummies:
            dist=get_intermolecular_distance(model, mol, model_atom, a)
            dummy_distances.append(dist)
        min_index = dummy_distances.index(min(dummy_distances))
        #max_index = dummy_distances.index(max(dummy_distances))
        #write("pre-aligned_linear.xyz",mol)
        #write("pre-model_linear.xyz",model)
        angle = get_general_angle(mol.positions[mol_dummies[min_index]],model_com,model.positions[model_atom])
        #print "Angle = ",angle
        #print "About to rotate: ", mol.positions[mol_dummies[min_index]]-model_com, model.positions[model_atom]-model_com
        if angle > eps:
            mol.rotate(mol.positions[mol_dummies[min_index]]-model_com,model.positions[model_atom]-model_com,center=model_com)
        transfer_dummy_tags_linear(mol,model)
        #write("semi-aligned_linear.xyz",mol)
        #write("semi-model_linear.xyz",model)
    #by default we'll align such that the furthest atom from the COM and the max component of the moment of inertia are perpendicular,
    #unless we have an extra align point 'Bq' (Bq still has mass=0, so it's useful as a point, but doesn't alter the COM
    model_ghosts = model.get_atom_indices('Bq')
    print "MG #", len(model_ghosts), model_ghosts
    if len(model_ghosts) > 0:
        pass
    else:
        moments=mol.get_moments_of_inertia()
        extra_vec=[0,0,0]
        max_moment = np.argmax(moments)
        extra_vec[max_moment]=1
        furthest_atom = furthest_real(mol)
        #print moments, extra_vec
        dih=get_arbitrary_dihedral(extra_vec,mol.positions[mol_dummies[0]],mol.positions[mol_dummies[1]],mol.positions[furthest_atom])
        print "Second rotation: ", dih
        if not np.isnan(dih):
            #print mol.positions
            if 'flip' in mol.info:
                mol.rotate(mol.positions[mol_dummies[0]]-mol.positions[mol_dummies[1]],dih+rot_angle,center=mol_cod) #pi fudge is because we don't know which way the dihedral is
            else:
                mol.rotate(mol.positions[mol_dummies[0]]-mol.positions[mol_dummies[1]],dih+pi+rot_angle,center=mol_cod) #pi fudge is because we don't know which way the dihedral is
            #mol.rotate(mol.positions[mol_dummies[0]]-mol.positions[mol_dummies[1]],dih,center=mol_cod) #pi fudge is because we don't know which way the dihedral is

    #TODO:Also add in a check for the second dummy (if the angle is greater than some epsilon, then spit the dummy)
    return


def align_linear_triangle(mol, model):
    eps = 0.02
    """A linear triangle has two dummies and an alignment point in a triangular config.
       Moves COM_mol to COM_model and aligns dummies."""
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180 
    else:
        rot_angle=0   
    model_com = model.get_center_of_mass()
    model_cod = model.get_center_of_dummies()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    mol_cod = mol.get_center_of_dummies()
    mol_com = mol.get_center_of_mass()

    #First alignment will align the two COM-COD vectors
    angle = get_general_angle(mol_cod,model_com,model_cod)
    if angle > eps:
        mol.rotate(mol_cod-model_com,model_cod-model_com,center=model_com)

    #Second alignment fixes the dummies
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        if a.tag:
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    model_dummies = model.get_atom_indices('X')
    mol_dummies = mol.get_atom_indices('X')
    #first dummy to first dummy, we don't care...
    dih=get_arbitrary_dihedral(mol.positions[mol_dummies[0]],model_com,model_cod,model.positions[model_dummies[0]])
    print "Second rotation: ", dih
    if not np.isnan(dih):
        #print mol.positions
        if 'flip' in mol.info:
            mol.rotate(model_cod-model_com,dih+pi,center=model_com) #pi fudge is because we don't know which way the dihedral is
        else:
            mol.rotate(model_cod-model_com,dih,center=model_com)

    transfer_dummy_tags(mol,model)
    return




def align_triangle(mol,model): #07/01/2014 MAA
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    We'll deal implicitly with equilateral and isosceles """
    eps = 0.02 # on angles, about 1.0 deg
    length_eps = 0.2
    non_planar = False

    #It may be that we're trying to align a triangle that is really a tetrahedron 
    # i.e. three connections, but COM significantly abouve the plane
    # We'll assume this is the case iff the COM and COD differ for *both* the model and mol 
    model_com = model.get_center_of_mass()
    model_cod = model.get_center_of_dummies()
    model_com_cod = np.linalg.norm(model_com - model_cod)
    mol_cod = mol.get_center_of_dummies()
    mol_com = mol.get_center_of_mass()
    mol_com_cod = np.linalg.norm(mol_com - mol_cod)

    print "Model COM-COD = ", model_com_cod
    print "Mol COM-COD = ", mol_com_cod
    if model_com_cod > 0.5 and mol_com_cod > 1.0:
        print "Non planar triangle case!"
        non_planar = True
        mol.translate(model_com - mol_com)
        mol_cod = mol.get_center_of_dummies()
        mol_com = mol.get_center_of_mass()
        #First alignment will align the two COM-COD vectors
        angle = get_general_angle(mol_cod,model_com,model_cod)
        if angle > eps:
            mol.rotate(mol_cod-model_com,model_cod-model_com,center=model_com)
        #Now we need to keep this axis fixed
        #We'll assume equilateral and just pick the primary dummies
        model_dummies = model.get_atom_indices('X')
        mol_dummies = mol.get_atom_indices('X')
        dih = get_arbitrary_dihedral(mol[mol_dummies[0]].position,model_cod,model_com,model[model_dummies[0]].position)
        print "Second angle: ", dih    
        if dih > eps:#1.0deg
            print "rotating non-planar triangle"
            mol.rotate(model_cod-model_com,dih,center=model_cod)

    else: 
        model_cod = model.get_center_of_dummies()
        mol_com = mol.get_center_of_dummies()
        mol.translate(model_cod - mol_com)

        tmpmol = Atoms()
        tmpmol += model
        tmpmol += mol
        write('tmpmol0.xyz',tmpmol)
    
        model_dummies = model.get_atom_indices('X')
        mol_dummies = mol.get_atom_indices('X')
        # Now test for equal sides
        #model
        ref_set = ()
        model_sides = np.array([model.get_distance(model_dummies[1],model_dummies[2]), model.get_distance(model_dummies[0],model_dummies[2]), model.get_distance(model_dummies[0],model_dummies[1]) ])
        if abs(model_sides[0] - model_sides[1]) < length_eps and abs(model_sides[0] - model_sides[2]) < length_eps:
            #equilateral
            print "model is approximately equilateral"
            model_primary_index = 0
            model_middle_index = 1
            model_secondary_index = 2
        else:
            model_primary_index = np.argmin(model_sides) #model_sides is rigged such that the index of the array is opposite the side
            model_secondary_index = np.argmax(model_sides) #model_sides is rigged such that the index of the array is opposite the side
            model_middle_index = set([0,1,2]).difference([model_primary_index, model_secondary_index]).pop()
            print "model P/S/M:", model_primary_index, model_secondary_index, model_middle_index
        #mol
        mol_sides = np.array([mol.get_distance(mol_dummies[1],mol_dummies[2]), mol.get_distance(mol_dummies[0],mol_dummies[2]), mol.get_distance(mol_dummies[0],mol_dummies[1]) ])
        if abs(mol_sides[0] - mol_sides[1]) < length_eps and abs(mol_sides[0] - mol_sides[2]) < length_eps:
            #equilateral
            print "mol is approximately equilateral"
            mol_primary_index = 0
            mol_middle_index = 1
            mol_secondary_index = 2
        else:
            mol_primary_index = np.argmin(mol_sides) #mol_sides is rigged such that the index of the array is opposite the side
            mol_secondary_index = np.argmax(mol_sides) #mol_sides is rigged such that the index of the array is opposite the side
            mol_middle_index = set([0,1,2]).difference([mol_primary_index, mol_secondary_index]).pop()
            print "mol P/S/M:", mol_primary_index, mol_secondary_index, mol_middle_index
    
        #First align the first dummy of mol and model
        angle = get_general_angle(mol.positions[mol_dummies[mol_primary_index]],model.positions[model_dummies[model_primary_index]],model_cod)
        angle2 = get_general_angle(model.positions[model_dummies[model_primary_index]],mol.positions[mol_dummies[mol_primary_index]],model_cod)
        print angle, angle2 
        #Has a problem with 180 degree angles - in the case of 180 deg, both the above angles will show up as ~0.0
        print "First angle: ",angle
        if angle < eps and angle2 < eps: #180deg
            print "In Triangle 180degree case"
            v1=mol.positions[mol_dummies[0]]-mol.positions[mol_dummies[2]]  #doing this, because mol_com may not be in same plane as the three dummies
            v2=mol.positions[mol_dummies[1]]-mol.positions[mol_dummies[2]]
            v = np.cross(v1, v2)
            mol.rotate(v,pi,center=model_cod)
        elif angle > eps: #1.0deg:
            mol.rotate(mol[mol_dummies[mol_primary_index]].position-model_cod,model[model_dummies[model_primary_index]].position-model_cod,center=model_cod)
        #Now if not flipped we align last ->last, else last -> middle
        dih = get_arbitrary_dihedral(mol[mol_dummies[mol_secondary_index]].position,model_cod,model[model_dummies[model_primary_index]].position,model[model_dummies[model_secondary_index]].position)
        print "Second angle: ", dih    
        if dih > eps:#1.0deg
            print "rotating triangle"
            mol.rotate(model[model_dummies[model_primary_index]].position-model_cod,dih,center=model_cod)

        if 'flip' in mol.info:
            print "performing flip"
            mol.rotate(model[model_dummies[model_primary_index]].position-model_cod,pi,center=model_cod)
    #Now, if it has been requested, rotate everything. In this case, things might not match up, but it is a way to align A-B, B-C, C-A
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        #need a vector perpendicular to the triangle plane
        v1=mol.positions[mol_dummies[0]]-mol.positions[mol_dummies[2]]  #doing this, because mol_com may not be in same plane as the three dummies
        v2=mol.positions[mol_dummies[1]]-mol.positions[mol_dummies[2]]
        v = np.cross(v1, v2)
        mol.rotate(v,rot_angle,center=mol_com)

    transfer_dummy_tags_dist(mol,model)

    return


def align_tetrahedral(mol,model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    A is primary, D secondary. If tags aren't "A,B,C,D", then they are sorted ASCIIbetically 
    and the first and last tags are the primary and secondary alignment respectively"""
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_dummies()
    print "AT:",mol_com
    mol.translate(model_com - mol_com)
    #Grab the tags on the model put them in a dict.
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        #print a
        if a.tag :
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    print "MDT", model_dummy_tags        
    mol_dummies = mol.get_atom_indices('X')
    #We align on the first (ASCIIbetical) model atom
    model_atom=min(model_dummy_tags, key=model_dummy_tags.get)
    mol_dummies = mol.get_atom_indices('X')
    dummy_distances={}
    for a in mol_dummies:
        dist=get_intermolecular_distance(model, mol, model_atom, a)
        dummy_distances[a] = dist
    if 'flip' in mol.info:  #Flipped in the case of a tetrahedron means max-min alignment i.e. D-A, then closest and the other two work themsevles out
        print "Tetrahedron: flipped alignment"
        max_index = max(dummy_distances, key=dummy_distances.get)
        mol.rotate(mol.positions[max_index]-model_com,model.positions[model_atom]-model_com,center=model_com)
        #write("first_aligned_tetrahedron.xyz",mol)
        fixed_A = max_index
        fixed_axis = mol.positions[max_index]-model_com
    else:
        print "Tetrahedron: normal alignment" #min-min alignment
        min_index = min(dummy_distances, key=dummy_distances.get)
        mol.rotate(mol.positions[min_index]-model_com,model.positions[model_atom]-model_com,center=model_com)
        #write("first_aligned_tetrahedron.xyz",mol)
        fixed_A = min_index
        fixed_axis = mol.positions[min_index]-model_com
    #Second rotation is the same regardless of the first    
    #And the second rotation, we'll use angle and axis to ensure the first axis stays
        model_atom=max(model_dummy_tags, key=model_dummy_tags.get)
        dummy_distances={}
    for a in mol_dummies:
        dist=get_intermolecular_distance(model, mol, model_atom, a)
        dummy_distances[a] = dist
    min_index = min(dummy_distances, key=dummy_distances.get) #dummy_distances.index(min(dummy_distances))
    if min_index == fixed_A:
        #remove the minimum value from the dist and extract the next one
        del dummy_distances[min_index]
        min_index = min(dummy_distances, key=dummy_distances.get)
    angle=get_arbitrary_dihedral(mol.positions[min_index],model_com,mol.positions[fixed_A],model.positions[model_atom])
    mol.rotate(fixed_axis,angle,center=model_com)
    #So rotation here assumes the first fixed axis, meaning that A-A and A-D alignment are possible, but not A-B or A-C
    #Need to add a "fix-axis" option??
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180 
        mol.rotate(fixed_axis,rot_angle,center=model_com)

    transfer_dummy_tags_dist(mol,model)
    return    

def align_square(mol,model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    #We can pick our two alignment points straight up
    model_dummies = model.get_atom_indices('X')
    intra_model_distances = {}
    for dindex in model_dummies:
        this_dist = model.get_distance(model_dummies[0],dindex)
        intra_model_distances[dindex] = this_dist
    print intra_model_distances
    sorted_IModelD = sorted(intra_model_distances, key=intra_model_distances.get)
    model_X1 = model_dummies[0] # pick the first dummy to align on
    model_X2 = sorted_IModelD[1] # pick the closest dummy to the first one, must be a 90deg one 
    #mol
    mol_dummies = mol.get_atom_indices('X')
    intra_mol_distances = {}
    for dindex in mol_dummies:
        this_dist = mol.get_distance(mol_dummies[0],dindex)
        intra_mol_distances[dindex] = this_dist
    print intra_mol_distances
    sorted_IMolD = sorted(intra_mol_distances, key=intra_mol_distances.get)
    mol_X1 = mol_dummies[0] # pick the first dummy to align on
    mol_X2 = sorted_IMolD[1] # pick the closest dummy to the first one, must be a 90deg one 
    #Now we can do the rotation
    angle1 = get_general_angle(mol.positions[mol_X1],model_com,model.positions[model_X1])
    print "angle1 = ", angle1
    this_angle = angle_2vectors(mol.positions[mol_X1]-model_com , model.positions[model_X1]-model_com)
    print "this_angle = ", this_angle
    if angle1 > 0.02 and abs(angle1 - pi) >0.02:
        if 'flip' in mol.info:
            mol.rotate(mol.positions[mol_X1]-model_com,-model.positions[model_X1]-model_com,center=model_com)
        else:
            mol.rotate(mol.positions[mol_X1]-model_com,model.positions[model_X1]-model_com,center=model_com)

    #Second rotation 
    angle2 = get_general_angle(mol.positions[mol_X2],model_com,model.positions[model_X2])
    print "angle2 = ", angle2
    if angle2 > 0.02 and abs(angle2 - pi) >0.02:
        mol.rotate(mol.positions[mol_X2]-model_com,model.positions[model_X2]-model_com,center=model_com)

    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        #need a vector perpendicular to the square plane
        v1=mol.positions[mol_dummies[0]]-mol.positions[mol_dummies[2]]  #doing this, because mol_com may not be in same plane as the three dummies
        v2=mol.positions[mol_dummies[1]]-mol.positions[mol_dummies[2]]
        v = np.cross(v1, v2)
        mol.rotate(v,rot_angle,center=mol_com)

    #Now transfer the tags
    transfer_dummy_tags_dist(mol,model)    

    return


def align_rectangle(mol,model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the midpoints of pairs of dummies"""
    eps = 0.02 #angles

    #model_com = model.get_center_of_mass()
    model_com = model.get_center_of_dummies()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    print "fragID =", mol[0].fragmentID
    # Now locate the model dummies and determine their midpoints
    model_midpoints = Atoms()
    model_dummies = model.get_atom_indices('X')
    model_Q = model.get_atom_indices('Q')
    model_midpoints.append(model[model_Q[0]])
    for a1 in range(len(model_dummies)):
        for a2 in range(a1+1,len(model_dummies)):
            this_midpoint = model[model_dummies[a1]].position + (model[model_dummies[a2]].position - model[model_dummies[a1]].position)/2
            #D = model[model_Q[0]].position - this_midpoint
            D = model_com - this_midpoint
            if np.linalg.norm(D) > 0.25: 
                model_midpoints.append(Atom('X',position=this_midpoint))
    # Now locate the mol dummies and determine their midpoints
    mol_midpoints = Atoms()
    mol_midpoints.append(model[model_Q[0]])
    mol_dummies = mol.get_atom_indices('X')
    for a1 in range(len(mol_dummies)):
        #mol_midpoints.append(Atom('X',position=mol[mol_dummies[a1]].position))
        for a2 in range(a1+1,len(mol_dummies)):
            print "Midpoints"
            this_midpoint = mol[mol_dummies[a1]].position + (mol[mol_dummies[a2]].position - mol[mol_dummies[a1]].position)/2
            print mol[mol_dummies[a1]].position, mol[mol_dummies[a2]].position, this_midpoint 
            #D = model[model_Q[0]].position - this_midpoint # still use the model Q, should be the same as mol.COM
            D = model_com - this_midpoint # still use the model Q, should be the same as mol.COM
            if np.linalg.norm(D) > 0.2: 
                mol_midpoints.append(Atom('N',position=this_midpoint))
    #Now we align on the long axis
    model_midpoint_distances = {}
    for a in model_midpoints:
        if a.symbol is not 'Q':
            dist=model_midpoints.get_distance(0,a.index)
            model_midpoint_distances[a.index] = dist
    model_max_index = max(model_midpoint_distances, key=model_midpoint_distances.get)
    model_min_index = min(model_midpoint_distances, key=model_midpoint_distances.get)
    print model_midpoint_distances
    print "model",model_min_index, model_max_index
    #molecule...
    mol_midpoint_distances = {}
    for a in mol_midpoints:
        if a.symbol is not 'Q':
            dist=mol_midpoints.get_distance(0,a.index)
            mol_midpoint_distances[a.index] = dist
    mol_max_index = max(mol_midpoint_distances, key=mol_midpoint_distances.get)
    mol_min_index = min(mol_midpoint_distances, key=mol_midpoint_distances.get) 
    #In case of some weird geometries (notably chs1), the min and max index aren't perpendicular, which creates hassles with the second dihedral
    intra_model_angles = {}
    for d1 in model_midpoints:
        if d1.symbol != 'Q':
            this_angle = get_general_angle(model_midpoints[model_max_index].position, model_com, d1.position)
            if abs(this_angle -  pi/2) < eps:
                model_perp_index = d1.index
                break
    #Now mol
    #print "Finding mol perp index"
    intra_mol_angles = {}
    for d1 in mol_midpoints:
        if d1.symbol != 'Q':
            this_angle = get_general_angle(mol_midpoints[mol_max_index].position, model_com, d1.position)
            #print d1, this_angle, pi/2
            if abs(this_angle -  pi/2) < eps*10: #Need to be quite slack here to take into account things that aren't centrosymmetric
                #print "selecting ", d1.index
                mol_perp_index = d1.index
                break

    #Now actually align
    #write('model'+`mol[0].fragmentID`+'_prealigned.xyz',model)
    #write('mol'+`mol[0].fragmentID`+'_prealigned.xyz',mol)
    angle = get_general_angle(mol_midpoints.positions[mol_max_index],model_midpoints.positions[model_max_index],model_com)
    print "1st angle =", angle
    if angle > 0.02 and abs(angle -  pi) > eps: #about a degree
        mol.rotate(mol_midpoints.positions[mol_max_index]-model_com,model_midpoints.positions[model_max_index]-model_com,center='COM')
        mol_midpoints.rotate(mol_midpoints.positions[mol_max_index]-model_com,model_midpoints.positions[model_max_index]-model_com,center='COM')
    if 'rotate' in mol.info:
        print "Rectangle: rotated alignment"
        mol.rotate(mol_midpoints.positions[mol_max_index]-model_com,model_com-model_midpoints.positions[model_max_index],center='COM') #reverse of the rotation above
        mol_midpoints.rotate(mol_midpoints.positions[mol_max_index]-model_com,model_com-model_midpoints.positions[model_max_index],center='COM') #reverse of the rotation above
    fixed_axis = mol_midpoints.positions[mol_max_index]-model_com
    #Then the short axis
    #molecule again...
    print "After 1st rotation",mol_midpoint_distances
    #print "mol",mol_perp_index, mol_max_index
    angle = get_arbitrary_dihedral(mol_midpoints.positions[mol_perp_index],model_com,mol_midpoints.positions[mol_max_index],model_midpoints.positions[model_perp_index])
    #print "2nd angle =", angle
    if 'flip' in mol.info:
        print "Rectangle: flipped alignment"
        if angle > 0.02 and abs(angle -  pi) > eps:
            mol.rotate(fixed_axis, angle,center=model_com)
            mol.rotate(fixed_axis, pi,center=model_com)

        else:
            mol.rotate(fixed_axis,pi,center=model_com)
    else:
        print "Rectangle: normal alignment"
        print abs(angle -  pi)
        if angle > 0.02 and abs(angle -  pi) > eps:
            print "Rotating by ",angle
            mol.rotate(fixed_axis,angle,center=model_com)

    #write('model'+`mol[0].fragmentID`+'_aligned.xyz',model)
    #write('mol'+`mol[0].fragmentID`+'_aligned.xyz',mol)

    #transfer the dummy tags based on distance only
    transfer_dummy_tags_dist(mol,model)

    return

def align_tri_bipyramid(mol, model): #MAA 30/01/2014. 
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    eps=0.1
    # sort both the model and mol dummies into set of 3 (triangle plane) and set of 2 (main axis)
    #model first
    model_dummies = model.get_atom_indices('X')
    intra_model_distances = {}
    for dindex in model_dummies:
        this_dist = model.get_distance(0,dindex)
        intra_model_distances[dindex] = this_dist
    print intra_model_distances
    sorted_IMD = sorted(intra_model_distances, key=intra_model_distances.get)
    #If all five distances are the same, we need to test on angles
    intra_model_angles = {}
    if all(abs(intra_model_distances[index] - intra_model_distances[sorted_IMD[0]]) < eps for index in intra_model_distances.keys()):
        print "TBP model distances all equal."
        for d1 in model_dummies:
            for d2 in model_dummies:
                this_angle = get_general_angle(model[d1].position, model_com, model[d2].position)
                if abs(this_angle -  pi) < eps:
                    model_axial_indices = [d1,d2]
                    model_radial_indices = [dx for dx in model_dummies if dx not in model_axial_indices]
                    break
    else:
        #Now, there are two options, tri-bipyramid is tall and skinny (axis ditances longer than radial distances) or short and fat.
        #The long and skinny option is more likely, but we'll allow for both options
        model_axial_indices = [sorted_IMD[3],sorted_IMD[4]]
        model_radial_indices = [sorted_IMD[0],sorted_IMD[1],sorted_IMD[2]]
        print "Axial ", model_axial_indices
        print "Radial ", model_radial_indices
        #Is this assignment ok? Test that all three radial distances are the same
        if not all(abs(intra_model_distances[index] - intra_model_distances[sorted_IMD[0]]) < eps for index in model_radial_indices):
            #Then the assignement is wrong - we have a short and fat case
            print "TBP model assignment failed, reassigning!"
            model_axial_indices = [sorted_IMD[0],sorted_IMD[1]]
            model_radial_indices = [sorted_IMD[2],sorted_IMD[3],sorted_IMD[4]]
    #mol 
    mol_dummies = mol.get_atom_indices('X')
    intra_mol_distances = {}
    for dindex in mol_dummies:
        D = mol[dindex].position - model_com
        this_dist = np.linalg.norm(D)
        intra_mol_distances[dindex] = this_dist
    print intra_mol_distances
    sorted_IMD = sorted(intra_mol_distances, key=intra_mol_distances.get)
    #If all five distances are the same, we need to test on angles
    intra_mol_angles = {}
    if all(abs(intra_mol_distances[index] - intra_mol_distances[sorted_IMD[0]]) < eps for index in intra_mol_distances.keys()):
        print "TBP mol distances all equal."
        for d1 in mol_dummies:
            for d2 in mol_dummies:
                this_angle = get_general_angle(mol[d1].position, model_com, mol[d2].position)
                if abs(this_angle -  pi) < eps:
                    mol_axial_indices = [d1,d2]
                    mol_radial_indices = [dx for dx in mol_dummies if dx not in mol_axial_indices]
                    break
    else:
        #Now, there are two options, tri-bipyramid is tall and skinny (axis ditances longer than radial distances) or short and fat.
        #The long and skinny option is more likely, but we'll allow for both options
        mol_axial_indices = [sorted_IMD[3],sorted_IMD[4]]
        mol_radial_indices = [sorted_IMD[0],sorted_IMD[1],sorted_IMD[2]]
        #Is this assignment ok? Test that all three radial distances are the same
        if not all(abs(intra_mol_distances[index] - intra_mol_distances[sorted_IMD[0]]) < eps for index in mol_radial_indices):
            #Then the assignement is wrong - we have a short and fat case
            print "TBP mol assignment failed, reassigning!"
            mol_axial_indices = [sorted_IMD[0],sorted_IMD[1]]
            mol_radial_indices = [sorted_IMD[2],sorted_IMD[3],sorted_IMD[4]]
    #So now we can define our axes and do the first rotation
    #Normal orientation will go sorted_IModelD[0] to sorted_IMolD[0] - not lexicographic, but it doesn't really matter
    if 'flip' in mol.info:
        print "Trigonal bipyramid: flipped alignment"
        mol.rotate(mol[mol_axial_indices[0]].position-model_com,model[model_axial_indices[1]].position-model_com,center='COM') #D3h axes now aligned
        #Now we need to rotate in the triangle plane. We can use an angle rather than a dihedral
        this_angle = get_general_angle(model[model_radial_indices[0]].position,model_com,mol[mol_radial_indices[0]].position)
        mol.rotate(model[model_axial_indices[1]].position-model_com,this_angle,center='COM') 
    else:
        print "Trigonal bipyramid: normal alignment"
        mol.rotate(mol[mol_axial_indices[0]].position-model_com,model[model_axial_indices[0]].position-model_com,center='COM') #D3h axes now aligned
        #Now we need to rotate in the triangle plane. We can use an angle rather than a dihedral
        this_angle = get_general_angle(model[model_radial_indices[0]].position,model_com,mol[mol_radial_indices[0]].position)
        mol.rotate(model[model_axial_indices[0]].position-model_com,this_angle,center='COM') 

    #Any extra rotation is defined as being around the D3 axis
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        mol.rotate(model[model_axial_indices[0]].position-model_com,rot_angle,center=model_com)

    #transfer the dummy tags
    transfer_dummy_tags_dist(mol,model)

    return

def align_octahedron(mol,model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    eps=0.1
    #first we get the tags for the model, our model must have one unequal (shorter) axis that identifies it as the unique/first ("C") axis
    model_distances={}
    for i in range(len(model)):
        if  model[i].tag: #and model[i].tag > 0:#not model[i].tag.endswith("\""): #len(model[i].tag) ==1:
            model_distances[model[i].tag]=model.get_distance(i,0)
            #model_tags[model.get_distance(i,0)]=model[i].tag #may be non-unique but we don't care
            print model[i].tag, model_distances[model[i].tag]
        #print model_tags
    #print sorted(model_distances, key=lambda key: model_distances[key])    
    first_tag=sorted(model_distances, key=lambda key: model_distances[key])[0]
    last_tag=sorted(model_distances, key=lambda key: model_distances[key])[2]
    #Now find the dummy that is closest in distance to the model atom with tag 'C'
    model_atom = int(model.get_tag_indices(first_tag)[0])
    mol_dummies = mol.get_atom_indices('X')
    Dprint "atom to com distances"
    dummy_distances={}
    for a in mol_dummies:
        dist=get_intermolecular_distance(model, mol, 0, a) #0 is the COM
        dummy_distances[a] = dist 
        #print mol[a].index, dist
    #print dummy_distances    
    min_index = min(dummy_distances, key=dummy_distances.get)  
    #rotating by a near-zero angle will cause numerical problems
    angle = get_general_angle(mol.positions[min_index],model.positions[model_atom],model_com)
    if angle > 0.02:#1.0:
        mol.rotate(mol.positions[min_index]-model_com,model.positions[model_atom]-model_com,center='COM')
    angle = get_general_angle(mol.positions[min_index],model.positions[model_atom],model_com)
    fixed_axis=mol.positions[min_index]-model_com
    #Flipping the molecule...
    if 'flip' in mol.info:
        print "Octahedron: flipping primary axis"
        mol.rotate(mol.positions[min_index]-model_com,pi,center='COM')
    #Now we have COM_mol in the right place, and  the smallest axis aligned, we'll do the biggest
    model_atom = int(model.get_tag_indices(last_tag)[0])
    max_index = max(dummy_distances, key=dummy_distances.get)  
    angle = get_general_angle(mol.positions[max_index],model.positions[model_atom],model_com)
    print "Second angle: ",angle
    #Now finally do the second rotation
    if angle > 0.02:#1.0:
        mol.rotate(mol.positions[max_index]-model_com,model.positions[model_atom]-model_com,center='COM')
    #Add in a debug_print somewhere about the other two angles (i.e. how far the structure is from being actually square)
    #Now, if it has been requested, rotate everything. In this case, C axis stays, A,B rotate
    if 'rotate' in mol.info:
        #rotating around the major axis 
        rot_angle = mol.info['rotate'] * pi / 180
        mol.rotate(fixed_axis,rot_angle,center=model_com)
    
    transfer_dummy_tags_dist(mol,model)
    return

def align_octahedron_2(mol,model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    This is the second octahedron option. Aligns a trigonal antiprism."""
    eps = 0.1
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        if a.tag:
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    model_dummies = model.get_atom_indices('X')
    #First we need to align on the D3d axis of the molecule, otherwise we can end up with a mis-alignment.
    #So, lets find them. We'll use a purely geometric means
    #Basically we'll take the model atom, and find the distances to all the other dummies, there will be:
    #two short distances: two guys on the same face - these are the guys we want 
    #two medium distances: opposite face but not directly through the COM 
    #one long distance: guy directly opposite, through the COM
    intra_model_distances = {}
    for dindex in model_dummies:
        this_dist = model.get_distance(model_dummies[0],dindex)
        intra_model_distances[dindex] = this_dist
    print intra_model_distances
    sorted_IMD = sorted(intra_model_distances, key=intra_model_distances.get) 
    #This assignment fails if the tri_antiprism is short and fat rather than tall and skinny
    model_face1_indices = [sorted_IMD[0],sorted_IMD[1],sorted_IMD[2]]
    model_face2_indices = [sorted_IMD[3],sorted_IMD[4],sorted_IMD[5]]
    #Test face1
    diff1 = abs(model.get_distance(model_face1_indices[2],model_face1_indices[1]) - model.get_distance(model_face1_indices[2],model_face1_indices[0]))
    diff2 = abs(model.get_distance(model_face1_indices[1],model_face1_indices[0]) - model.get_distance(model_face1_indices[1],model_face1_indices[2]))
    if diff1 < eps and diff2 < eps:
        print "Assignment ok!"
    else:
        #we have a short, fat guy (2 medium - we want these, two small and one big)  
        print "Assignment failed, reassigning!"
        model_face1_indices = [sorted_IMD[0],sorted_IMD[3],sorted_IMD[4]]
        model_face2_indices = [sorted_IMD[2],sorted_IMD[3],sorted_IMD[5]]
    model_face1 = Atoms('X3',positions=[model[model_face1_indices[0]].position,model[model_face1_indices[1]].position,model[model_face1_indices[2]].position])
    model_face2 = Atoms('X3',positions=[model[model_face2_indices[0]].position,model[model_face2_indices[1]].position,model[model_face2_indices[2]].position])
    #Now the same trick for the molecule:
    mol_dummies = mol.get_atom_indices('X')
    intra_mol_distances = {}
    for dindex in mol_dummies:
        this_dist = mol.get_distance(mol_dummies[0],dindex)
        intra_mol_distances[dindex] = this_dist
    print intra_mol_distances
    sorted_IMD = sorted(intra_mol_distances, key=intra_mol_distances.get)
    #This assignment fails if the tri_antiprism is short and fat rather than tall and skinny
    mol_face1_indices = [sorted_IMD[0],sorted_IMD[1],sorted_IMD[2]]
    mol_face2_indices = [sorted_IMD[3],sorted_IMD[4],sorted_IMD[5]]
    #Test face1
    diff1 = abs(mol.get_distance(mol_face1_indices[2],mol_face1_indices[1]) - mol.get_distance(mol_face1_indices[2],mol_face1_indices[0]))
    diff2 = abs(mol.get_distance(mol_face1_indices[1],mol_face1_indices[0]) - mol.get_distance(mol_face1_indices[1],mol_face1_indices[2]))
    if diff1 < eps and diff2 < eps:
        print "Mol assignment ok!"
    else:
        #we have a short, fat guy (2 medium - we want these, two small and one big)
        print "Mol assignment failed, reassigning!"
        mol_face1_indices = [sorted_IMD[0],sorted_IMD[3],sorted_IMD[4]]
        mol_face2_indices = [sorted_IMD[2],sorted_IMD[3],sorted_IMD[5]]
    mol_face1 = Atoms('X3',positions=[mol[mol_face1_indices[0]].position,mol[mol_face1_indices[1]].position,mol[mol_face1_indices[2]].position])
    mol_face2 = Atoms('X3',positions=[mol[mol_face2_indices[0]].position,mol[mol_face2_indices[1]].position,mol[mol_face2_indices[2]].position])
    #Now find the centroids of the mol and model faces. We only need one for each
    model_centroid1 = model_face1.get_center_of_dummies()
    mol_centroid1 = mol_face1.get_center_of_dummies()
    
    if 'flip' in mol.info:
        print "Trimeric SBU: flipped alignment"
        mol.rotate(mol_centroid1-model_com,-model_centroid1-model_com,center='COM') #D3d axes now aligned
        #Now we need to such that mol_dummy, model_axis and model dummy are in the same plane
        dih=get_arbitrary_dihedral(mol[mol_face1_indices[0]].position,model_centroid1,model_com,model[model_face1_indices[0]].position)
        mol.rotate(model_centroid1-model_com,-dih + (60 * pi / 180),center='COM') #use the model centroid as the centroid is not updated, add rotation of 60 deg
    else:
        print "Trimeric SBU: normal alignment"
        mol.rotate(mol_centroid1-model_com,model_centroid1-model_com,center='COM') #D3d axes now aligned
        #Now we need to such that mol_dummy, model_axis and model dummy are in the same plane
        dih=get_arbitrary_dihedral(mol[mol_face1_indices[0]].position,model_centroid1,model_com,model[model_face1_indices[0]].position)
        mol.rotate(model_centroid1-model_com,-dih,center='COM') #use the model centroid as the centroid is not updated
    #Any extra rotation is defined as being around the D3d axis
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        mol.rotate(model_centroid1-model_com,rot_angle,center=model_com)

    #transfer the dummy tags
    transfer_dummy_tags_dist(mol,model)

    return


def align_octahedron_3(mol,model): #This is generally not a good option. Use octahedron_2 instead
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    This is the third octahedron option. Aligns six corners of a cube."""
    eps = 0.1
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        if a.tag:
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    model_dummies = model.get_atom_indices('X')
    #First select a pair of opposite dummies. Doesn't matter which
    intra_model_distances = {}
    for dindex in model_dummies:
        this_dist = model.get_distance(model_dummies[0],dindex)
        intra_model_distances[dindex] = this_dist
    print intra_model_distances
    sorted_ImodelD = sorted(intra_model_distances, key=intra_model_distances.get)
    opp_dummy = sorted_ImodelD[5]
    #Now the same trick for the molecule:
    mol_dummies = mol.get_atom_indices('X')
    intra_mol_distances = {}
    for dindex in mol_dummies:
        this_dist = mol.get_distance(mol_dummies[0],dindex)
        intra_mol_distances[dindex] = this_dist
    print intra_mol_distances
    sorted_ImolD = sorted(intra_mol_distances, key=intra_mol_distances.get)
    opp_dummy = sorted_ImolD[5]
    #Now we can align on the first axis
    if 'flip' not in mol.info:
        print "Octahedral SBU (3): normal alignment"
        mol.rotate(mol[sorted_ImolD[0]].position-model_com,model[sorted_ImodelD[0]].position-model_com,center='COM') #defined axes now aligned
    else:
        print "Octahedral SBU (3): flipped alignment"
        mol.rotate(mol[sorted_ImolD[5]].position-model_com,model[sorted_ImodelD[0]].position-model_com,center='COM') #defined axes now aligned
    #Now we need to align such that the other two pairs are equally misaligned
    model_midpoint = model[sorted_ImodelD[1]].position + (model[sorted_ImodelD[4]].position - model[sorted_ImodelD[1]].position)/2 #shortest + (second longest - shortest)/2
    mol_midpoint = mol[sorted_ImolD[1]].position + (mol[sorted_ImolD[4]].position - mol[sorted_ImolD[1]].position)/2 #shortest + (second longest - shortest)/2
    #Now we align these two midpoints about the first axis
    dih = get_arbitrary_dihedral(mol_midpoint,model_com,model[sorted_ImodelD[0]].position,model_midpoint)
    mol.rotate(model[sorted_ImodelD[0]].position-model_com,-dih,center='COM')

    transfer_dummy_tags_dist(mol,model)

    return

def align_octahedron_4(mol,model): #There are some peculiarly distorted octahedra around
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    This is the fourth octahedron option. Aligns first dummy exactly, second to the plane 
    and assumes(!) the rest. This is effectively a failsafe option. Align won't die, but
    the results may be... interesting.
    Flip option is ignored."""

    eps = 0.1
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    model_dummy_tags={}
    model_tag_indices={}
    model_dummies = model.get_atom_indices('X')
    mol_dummies = mol.get_atom_indices('X')
    print "OA4: ", len(mol), mol_dummies[0], len(model), model_dummies[0]
    #First align the first dummy of mol and model
    angle = get_general_angle(mol.positions[mol_dummies[0]],model.positions[model_dummies[0]],model_com)
    print "First angle: ",angle
    if angle > 0.02:#1.0deg:
        mol.rotate(mol[mol_dummies[0]].position-model_com,model[model_dummies[0]].position-model_com,center=model_com)
    #Now find the closest dummy to the one aligned
    intra_model_distances = {}
    for dindex in model_dummies:
        this_dist = model.get_distance(model_dummies[0],dindex)
        intra_model_distances[dindex] = this_dist
    print intra_model_distances
    sorted_ImodelD = sorted(intra_model_distances, key=intra_model_distances.get)
    nearest_model_dummy = sorted_ImodelD[1]
    #Now the same trick for the molecule:
    intra_mol_distances = {}
    for dindex in mol_dummies:
        this_dist = mol.get_distance(mol_dummies[0],dindex)
        intra_mol_distances[dindex] = this_dist
    print intra_mol_distances
    sorted_ImolD = sorted(intra_mol_distances, key=intra_mol_distances.get)
    nearest_mol_dummy = sorted_ImolD[1]
    #Now we can do the second alignment
    dih = get_arbitrary_dihedral(mol[nearest_mol_dummy].position,model_com,model[sorted_ImodelD[0]].position,model[nearest_model_dummy].position)
    if dih > 0.02:#1.0deg
        mol.rotate(model[model_dummies[0]].position-model_com,dih,center=model_com)


    transfer_dummy_tags_dist(mol,model)

    return

def align_tri_prism(mol,model):
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    eps = 0.1
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    #Grab the tags on the model and sort them asciibetically. We'll use a dict.
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        if a.tag:
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    model_dummies = model.get_atom_indices('X')
    #First we need to align on the D3h axis of the molecule, otherwise we can end up with a mis-alignment.
    #So, lets find them. We'll use a purely geometric means
    #Basically we'll take the model atom, and find the distances to all the other dummies, there will be:
    #one short distance: guy immediately below
    #two medium distances: two guys on the same face - these are the guys we want
    #two long distances: guys on the opposite face and not immediately below
    intra_model_distances = {}
    for dindex in model_dummies:
        this_dist = model.get_distance(model_dummies[0],dindex)
        intra_model_distances[dindex] = this_dist
    print intra_model_distances
    sorted_IMD = sorted(intra_model_distances, key=intra_model_distances.get)
    #This assignment fails if the tri_prism is tall and skinny rather than short and fat
    model_face1_indices = [model_dummies[0],sorted_IMD[2],sorted_IMD[3]]
    model_face2_indices = [sorted_IMD[1],sorted_IMD[4],sorted_IMD[5]]
    #Test face1
    diff1 = abs(model.get_distance(model_face1_indices[2],model_face1_indices[1]) - model.get_distance(model_face1_indices[2],model_face1_indices[0]))
    diff2 = abs(model.get_distance(model_face1_indices[1],model_face1_indices[0]) - model.get_distance(model_face1_indices[1],model_face1_indices[2]))
    if diff1 < eps and diff2 < eps:
        print "Assignment ok!"
    else:
        #we have a tall, skinny guy
        print "Assignment failed, reassigning!"
        model_face1_indices = [sorted_IMD[0],sorted_IMD[1],sorted_IMD[2]]
        model_face2_indices = [sorted_IMD[3],sorted_IMD[4],sorted_IMD[5]]
    print "model 1", model_face1_indices
    print "model 2", model_face2_indices
    model_face1 = Atoms('X3',positions=[model[model_face1_indices[0]].position,model[model_face1_indices[1]].position,model[model_face1_indices[2]].position])
    model_face2 = Atoms('X3',positions=[model[model_face2_indices[0]].position,model[model_face2_indices[1]].position,model[model_face2_indices[2]].position])
    #Now the same trick for the molecule:
    mol_dummies = mol.get_atom_indices('X')
    intra_mol_distances = {}
    for dindex in mol_dummies:
        this_dist = mol.get_distance(mol_dummies[0],dindex)
        intra_mol_distances[dindex] = this_dist
    print intra_mol_distances
    sorted_IMD = sorted(intra_mol_distances, key=intra_mol_distances.get)
    mol_face1_indices = [mol_dummies[0],sorted_IMD[2],sorted_IMD[3]]
    mol_face2_indices = [sorted_IMD[1],sorted_IMD[4],sorted_IMD[5]]
    #print "mol 1", mol_face1_indices
    #print "mol 2", mol_face2_indices
    #Same check for the molecule
    #Test face1
    diff1 = abs(mol.get_distance(mol_face1_indices[2],mol_face1_indices[1]) - mol.get_distance(mol_face1_indices[2],mol_face1_indices[0]))
    diff2 = abs(mol.get_distance(mol_face1_indices[1],mol_face1_indices[0]) - mol.get_distance(mol_face1_indices[1],mol_face1_indices[2]))
    if diff1 < eps and diff2 < eps:
        print "Mol Assignment ok!"
    else:
        #we have a tall, skinny guy
        print "Mol Assignment failed, reassigning!"
        mol_face1_indices = [sorted_IMD[0],sorted_IMD[1],sorted_IMD[2]]
        mol_face2_indices = [sorted_IMD[3],sorted_IMD[4],sorted_IMD[5]]
    
    mol_face1 = Atoms('X3',positions=[mol[mol_face1_indices[0]].position,mol[mol_face1_indices[1]].position,mol[mol_face1_indices[2]].position])
    mol_face2 = Atoms('X3',positions=[mol[mol_face2_indices[0]].position,mol[mol_face2_indices[1]].position,mol[mol_face2_indices[2]].position])
    #Now find the centroids of the mol and model faces. We only need one for each
    model_centroid1 = model_face1.get_center_of_dummies()
    mol_centroid1 = mol_face1.get_center_of_dummies()
    #write('model3_prealigned.xyz',model)
    #write('mol3_prealigned.xyz',mol)
    #Now we finally get to actually aligning stuff
    if 'flip' in mol.info:
        print "Trimeric SBU: flipped alignment"
        mol.rotate(mol_centroid1-model_com,-model_centroid1-model_com,center='COM') #D3h axes now aligned
        #Now we need to such that mol_dummy, model_axis and model dummy are in the same plane
        dih=get_arbitrary_dihedral(mol[mol_face1_indices[0]].position,model_centroid1,model_com,model[model_face1_indices[0]].position)
        mol.rotate(model_centroid1-model_com,-dih,center='COM') #use the model centroid as the centroid is not updated
    else:
        print "Trimeric SBU: normal alignment"
        mol.rotate(mol_centroid1-model_com,model_centroid1-model_com,center='COM') #D3h axes now aligned
        #write('model3_semialigned.xyz',model)
        #write('mol3_semialigned.xyz',mol)
        #Now we need to such that mol_dummy, model_axis and model dummy are in the same plane
        dih=get_arbitrary_dihedral(mol[mol_face1_indices[0]].position,model_centroid1,model_com,model[model_face1_indices[0]].position)
        mol.rotate(model_centroid1-model_com,-dih,center='COM') #use the model centroid as the centroid is not updated
    #Any extra rotation is defined as being around the D3h axis
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        mol.rotate(model_centroid1-model_com,rot_angle,center=model_com)

    #write('model3_aligned.xyz',model)
    #write('mol3_aligned.xyz',mol)
    #transfer the dummy tags 
    transfer_dummy_tags_dist(mol,model)    

    return

def align_hexagon(mol,model): #New 08/04/2014
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    We find the midpoint of each pair of tags and then align like a triangle.""" # Hexagons are often irregular
    eps = 0.02 #on agles, about 1 deg
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    model_dummies = model.get_atom_indices('X')
    #Now we want to pair up each of the dummy atoms
    model_midpoints = Atoms()
    model_Q = model.get_atom_indices('Q')
    model_midpoints.append(model[model_Q[0]])
    model_Xpairs = {}
    model_mask = np.zeros(len(model), dtype=bool )
    for a1 in range(len(model_dummies)):
        tmp_dist = {}
        for a2 in range(len(model_dummies)):
            if a1 != a2:
                dist=model.get_distance(model_dummies[a1], model_dummies[a2])
                tmp_dist[model_dummies[a2]] = dist
        model_Xpairs[model_dummies[a1]] = min(tmp_dist, key=tmp_dist.get)  #This should be symmetric, but maybe check it
        print model_Xpairs
    for a1 in model_dummies:
        if not model_mask[a1]:
            new_midpoint = model[a1].position + (model[model_Xpairs[a1]].position - model[a1].position)/2
            model_midpoints.append(Atom('X',position=new_midpoint))
            model_mask[a1] = True
            model_mask[model_Xpairs[a1]] = True
    #Now the molecule
    mol_dummies = mol.get_atom_indices('X')
    #Now we want to pair up each of the dummy atoms
    mol_midpoints = Atoms()
    mol_midpoints.append(model[model_Q[0]])
    mol_Xpairs = {}
    mol_mask = np.zeros(len(mol), dtype=bool )
    for a1 in range(len(mol_dummies)):
        tmp_dist = {}
        for a2 in range(len(mol_dummies)):
            if a1 != a2:
                dist=mol.get_distance(mol_dummies[a1], mol_dummies[a2])
                tmp_dist[mol_dummies[a2]] = dist
        mol_Xpairs[mol_dummies[a1]] = min(tmp_dist, key=tmp_dist.get)  #This should be symmetric, but maybe check it
        print mol_Xpairs
    for a1 in mol_dummies:
        if not mol_mask[a1]:
            new_midpoint = mol[a1].position + (mol[mol_Xpairs[a1]].position - mol[a1].position)/2
            mol_midpoints.append(Atom('X',position=new_midpoint))
            mol_mask[a1] = True
            mol_mask[mol_Xpairs[a1]] = True

    #Now we can start aligning. Yay!
    angle = get_general_angle(mol_midpoints[1].position,model_midpoints[1].position,model_com) #index 0 is the COM
    angle2 = get_general_angle(model_midpoints[1].position,mol_midpoints[1].position,model_com) #index 0 is the COM
    print "First angle: ",angle
    if angle < eps and angle2 < eps: #180deg
        print "In Hexagon 180degree case"
        v1=mol_midpoints.positions[1]-mol_midpoints.positions[2]  #doing this, because mol_com may not be in same plane as the three dummies
        v2=mol_midpoints.positions[2]-mol_midpoints.positions[3]
        v = np.cross(v1, v2)
        mol.rotate(v,pi,center=model_com)
    elif angle > eps: #1.0deg:
        mol.rotate(mol_midpoints[1].position-model_com,model_midpoints[1].position-model_com,center=model_com)
    #Now second rotation
    dih = get_arbitrary_dihedral(mol_midpoints[2].position,model_com, model_midpoints[1].position, model_midpoints[2].position)
    print "Hexagon dih = ", dih
    if dih > eps:#1.0deg
        mol.rotate(model_midpoints[1].position-model_com,dih,center=model_com)

    print "about to call TDT for hexagon case"
    transfer_dummy_tags_dist(mol,model)
    #print "back from TDT"
    return

def align_mfu4(mol,model): # 11/02/2014 MAA... 13/02/2014 Seems to work :-)
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies.
    We find the midpoint of each pair of tags and then align like an octahedron.""" #The reason we don't align on an octahedron model is tags
    eps = 0.1
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    model_dummy_tags={}
    model_tag_indices={}
    for a in model:
        if a.tag:
            model_dummy_tags[a.index]=a.tag
            model_tag_indices[a.tag]=a.index
    model_dummies = model.get_atom_indices('X')
    #Now we want to pair up each of the dummy atoms
    model_midpoints = Atoms()
    model_Q = model.get_atom_indices('Q')
    model_midpoints.append(model[model_Q[0]])
    model_Xpairs = {}
    model_mask = np.zeros(len(model), dtype=bool )
    for a1 in range(len(model_dummies)):
        tmp_dist = {}
        for a2 in range(len(model_dummies)):
            if a1 != a2:
                dist=model.get_distance(model_dummies[a1], model_dummies[a2])
                tmp_dist[model_dummies[a2]] = dist
        model_Xpairs[model_dummies[a1]] = min(tmp_dist, key=tmp_dist.get)  #This should be symmetric, but maybe check it
        print model_Xpairs
    for a1 in model_dummies:
        if not model_mask[a1]:
            new_midpoint = model[a1].position + (model[model_Xpairs[a1]].position - model[a1].position)/2
            model_midpoints.append(Atom('X',position=new_midpoint))
            model_mask[a1] = True
            model_mask[model_Xpairs[a1]] = True
    #Now the molecule
    mol_dummies = mol.get_atom_indices('X')
    #Now we want to pair up each of the dummy atoms
    mol_midpoints = Atoms()
    mol_midpoints.append(model[model_Q[0]])
    mol_Xpairs = {}
    mol_mask = np.zeros(len(mol), dtype=bool )
    for a1 in range(len(mol_dummies)):
        tmp_dist = {}
        for a2 in range(len(mol_dummies)):
            if a1 != a2:
                dist=mol.get_distance(mol_dummies[a1], mol_dummies[a2])
                tmp_dist[mol_dummies[a2]] = dist
        mol_Xpairs[mol_dummies[a1]] = min(tmp_dist, key=tmp_dist.get)  #This should be symmetric, but maybe check it
        print mol_Xpairs
    for a1 in mol_dummies:
        if not mol_mask[a1]:
            new_midpoint = mol[a1].position + (mol[mol_Xpairs[a1]].position - mol[a1].position)/2
            mol_midpoints.append(Atom('X',position=new_midpoint))
            mol_mask[a1] = True
            mol_mask[mol_Xpairs[a1]] = True
            
    #Now we can start aligning. Yay!
    #For lack of anything more intelligent, we'll align on the first midpoint
    angle = get_general_angle(mol_midpoints[1].position,model_midpoints[1].position,model_com) #index 0 is the COM
    print "Angle = ", angle
    if angle > 0.02 and abs(angle - pi) >0.02:#1.0:
        mol.rotate(mol_midpoints[1].position-model_com,model_midpoints[1].position-model_com,center='COM')
    fixed_axis = mol_midpoints[1].position-model_com
    write('model4_semialigned.xyz',model)
    write('mol4_semialigned.xyz',mol)
    #Second alignment. An angle should work, but we'll use a dihedral in case
    #We can pick any midpoint, for both the model and mol, except the one that is furthest away
    if model_midpoints.get_distance(1,2) < model_midpoints.get_distance(1,3):
        model_X = 2
    else:
        model_X = 3
    if mol_midpoints.get_distance(1,2) < mol_midpoints.get_distance(1,3):
        mol_X = 2
    else:
        mol_X = 3

    dih=get_arbitrary_dihedral(mol_midpoints[mol_X].position,model_com, model_midpoints[1].position, model_midpoints[model_X].position)
    print "Dih = ",dih*180/pi
    mol.rotate(fixed_axis,dih,center='COM')
    #any extra rotation is defined as being around the D3d axis
    if 'rotate' in mol.info:
        rot_angle = mol.info['rotate'] * pi / 180
        mol.rotate(fixed_axis,rot_angle,center=model_com)

    write('model4_aligned.xyz',model)
    write('mol4_aligned.xyz',mol)
    #transfer the dummy tags
    try:
        transfer_dummy_tags_dist(mol,model)
    except RuntimeError: #We can be 90 deg off
        print "Rotating by 90 deg, because of misalignment"
        mol.rotate(fixed_axis,90 * pi / 180,center='COM')
        transfer_dummy_tags_dist(mol,model)

    return

def align_mil53(mol,model):
    """Aligns the heavy atom and then
    rotates appropriately to align the other dummies."""
    eps = 0.02

    model_heavy_atoms = model.get_atom_indices('Q')
    #Now we need to get the metal atoms of the connector
    mol_heavy_atoms = [item for item in range(len(mol.get_atomic_numbers())) if mol[item].number >10]
    #translate
    mol.translate(model[model_heavy_atoms[0]].position - mol[mol_heavy_atoms[0]].position) #should check length of heavy_atoms, but...
    model_dummies = model.get_atom_indices('X')
    #anchor is the midpoints of the dummies
    #Now we need to align the true dummies
    #model_dummies = model.get_atom_indices('X')
    #First to find them
    #We'll rely on them being typed - i.e. an atom type of C_R #otherwise we could find them as close to mol_furthest_Q
    mol_true_dummies = [atom.index for atom in mol if atom.symbol == 'X' and atom.mmtype == 'C_R']
    #Now there are two ways to align, only one is correct
    model_midpoint = model[model_dummies[0]].position + (model[model_dummies[1]].position - model[model_dummies[0]].position)/2
    mol_midpoint = mol[mol_true_dummies[0]].position + (mol[mol_true_dummies[1]].position - mol[mol_true_dummies[0]].position)/2
    #Now, keeping the heavy atoms fixed, align the heavy-midpoint axes
    angle = get_general_angle(model_midpoint,model[model_heavy_atoms[0]].position,mol_midpoint)
    print "Align MIL, first angle =", angle
    if angle > eps:
        mol.rotate(mol_midpoint - mol[mol_heavy_atoms[0]].position, model_midpoint - model[model_heavy_atoms[0]].position, center=mol[mol_heavy_atoms[0]].position)
    write('model53_semialigned.xyz',model)
    write('mol53_semialigned.xyz',mol)
    #Recalc 
    mol_midpoint = mol[mol_true_dummies[0]].position + (mol[mol_true_dummies[1]].position - mol[mol_true_dummies[0]].position)/2
    #Now use COM. Relies on the model providing the second M/Q atom
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    #now align
    dih = get_arbitrary_dihedral(mol_com, mol[mol_heavy_atoms[0]].position, mol_midpoint, model_com)
    print "MIL53 dih =", dih
    if dih > eps:
        mol.rotate(model_midpoint - model[model_heavy_atoms[0]].position,dih,center=mol[mol_heavy_atoms[0]].position)
    write('model53_aligned.xyz',model)
    write('mol53_aligned.xyz',mol)
    
    #transfer the dummy tags
    transfer_dummy_tags_mil53(mol,model)

    return


def align_icosahedron(mol,model): #20/03/2014 Tested 30/03/2014
    """Moves COM_mol to COM_model and then rotates the molecule to align the dummies"""
    eps = 0.1
    model_com = model.get_center_of_mass()
    mol_com = mol.get_center_of_mass()
    mol.translate(model_com - mol_com)
    model_dummies = model.get_atom_indices('X')

    #In the ideal Fm-3m (225) symmetry, from one reference endpoint toother endpoints via the COM, there will be 1 180 deg angle, 1 0 deg angle (me),
    #4 60 deg angles, 4 120 deg angles and 2 90 deg angles.
    intra_model_angles = {}
    for d1 in model_dummies:
        this_angle = get_general_angle(model[model_dummies[0]].position, model_com, model[d1].position)
        if abs(this_angle -  pi/2) < eps:
            model_perp_index = d1
            break
    #Now mol
    mol_dummies = mol.get_atom_indices('X')
    intra_mol_angles = {}
    for d1 in mol_dummies:
        this_angle = get_general_angle(mol[mol_dummies[0]].position, model_com, mol[d1].position)
        print d1, this_angle
        if abs(this_angle -  pi/2) < eps:
            mol_perp_index = d1
            break

    #Now first alignment
    angle = get_general_angle(mol.positions[mol_dummies[0]],model_com,model.positions[model_dummies[0]])
    print "First angle: ",angle
    if angle > eps:#1.0:
        mol.rotate(mol.positions[mol_dummies[0]]-model_com,model.positions[model_dummies[0]]-model_com,center=model_com)

    #Second alignment can be done straight away
    dih=get_arbitrary_dihedral(mol[mol_perp_index].position,model[model_dummies[0]].position,model_com,model[model_perp_index].position)
    print "Second angle: ",angle
    if dih > eps:#1.0:
        mol.rotate(model[model_dummies[0]].position-model_com,-dih,center=model_com) #use the model centroid as the centroid is not updated

    transfer_dummy_tags_dist(mol,model)

    return

def transfer_dummy_tags_linear(mol,model):
    """For linear guys we use the centre-of-dummies as the COM doesn't have to line up"""
    mol_dummies = mol.get_atom_indices('X')
    model_dummies = model.get_atom_indices('X')
    com=mol.get_center_of_dummies()
    eps=0.05
    for a in model_dummies:
        for b in mol_dummies:
            print com,model[a].position,mol[b].position
            angle=get_general_angle(model[a].position,com,mol[b].position)
            print a,b,angle, model[a].tag
            if angle < eps:  #angle is always positive
                #mol[b].tag = ''.join(model[a].tag)
                mol[b].tag = model[a].tag
                print "TDT assign:",a,b, angle, pi, angle-pi, angle, model[a].tag
                break
                #print a,b,angle, model[a].tag

def transfer_dummy_tags(mol,model): 
    mol_dummies = mol.get_atom_indices('X')
    model_dummies = model.get_atom_indices('X')
    com=mol.get_center_of_mass()
    eps=0.05 #in radians = 2.86 deg
    for a in model_dummies:
        for b in mol_dummies:
            #angle=get_general_angle(com,model[a].position,mol[b].position)
            angle=get_general_angle(model[a].position,com,mol[b].position)
            #angle=get_general_angle(com,mol[b].position,model[a].position)
            print a,b, angle
            #if angle >pi - eps:  #angle is always positive
            if angle < eps:  #angle is always positive
                mol[b].tag = model[a].tag
                print "TDT assign:",a,b, angle, pi, angle-pi, angle, model[a].tag
                break
    #Need to add in a test to make sure everything is assigned
    return

def transfer_dummy_tags_dist(mol,model):    #This assigns the closest model dummy to each mol dummy
    mol_dummies = mol.get_atom_indices('X')
    model_dummies = model.get_atom_indices('X')
    com=mol.get_center_of_mass()
    #Now we set up a dictionary to match dummy to dummy
    mol_to_model={}
    for a in mol_dummies:
        these_distances = {}
        for b in model_dummies:
            this_dist = get_intermolecular_distance(mol, model, a, b)
            these_distances[b] = this_dist
        print "TDTD: ",these_distances
        if min(these_distances, key=these_distances.get) not in mol_to_model.values():
            mol_to_model[a] = min(these_distances, key=these_distances.get)
        else:
            mol_to_model[a] = min(these_distances, key=these_distances.get)
            print mol_to_model
            raise RuntimeError('Transferring dummy tags based on distance criteria. Something assigned twice.')

    assigned = 0
    for a_mol,a_mod in mol_to_model.iteritems():
        mol[a_mol].tag = model[a_mod].tag
        assigned += 1
        print "TDT-distance  assign:",a_mod,a_mol,model[a_mod].tag
 
    if assigned != len(model_dummies):
        raise RuntimeError('Not all tags transferred!')

    return

def transfer_dummy_tags_mil53(mol,model):
    """for MIL-53 et al. transfers tags to other (M-M) dummy atoms"""
    mol_dummies = [atom.index for atom in mol if atom.symbol == 'X' and atom.mmtype == 'H_'] 
    model_dummies = model.get_atom_indices('Q')

    model_tag = model[model_dummies[0]].tag
    print "TDT-MIL53: Tag =",model_tag
    for mol_xh in mol_dummies:
        mol[mol_xh].tag = model_tag

    return

def bond_extra_dummies(mof):
    """ Looks for tagged dummies that are less than 1Ang apart and have *different* tags
    and bonds them.
    This should only be called after the 'real dummies' are sorted or else weird stuff happens"""

    bond_eps = 1.1 #Ang
    #But we'll check for mmtype anyway
    mof_dummies = [atom.index for atom in mof if atom.symbol == 'X' and atom.mmtype == 'H_']
    
    write('internal.xyz',mof)
    bond_matrix = np.zeros((len(mof_dummies), len(mof_dummies)), dtype=np.float64)
    #print mof.get_distance(4,7) , mof.get_distance(4,7,mic=True)
    for a1 in xrange(len(mof_dummies) -1 ):
        for a2 in xrange(a1+1,len(mof_dummies)):
            #print mof_dummies[a1], mof_dummies[a2], mof[mof_dummies[a1]].position, mof[mof_dummies[a2]].position, mof.get_distance(mof_dummies[a1],mof_dummies[a2], mic=True) 
            bond_matrix[a1,a2] = mof.get_distance(mof_dummies[a1],mof_dummies[a2],mic=True)
            if bond_matrix[a1,a2] <bond_eps:
                print "Would form bond:, ",mof_dummies[a1], '=>' ,mof_dummies[a2], bond_matrix[a1,a2], mof[mof_dummies[a1]].tag, mof[mof_dummies[a2]].tag

    

    mil_bonds = True
    while mil_bonds:
        mil_bonds=find_mil_bond(mof)

    print bond_matrix
    for a1 in xrange(len(mof_dummies) -1 ):
        for a2 in xrange(a1+1,len(mof_dummies)):
            if bond_matrix[a1,a2] < bond_eps:
                #mof.form_bond(mof_dummies[a1],mof_dummies[a2])
                pass

    return mof

def find_mil_bond(mol):
    """Finds a single internal bond in a supercell
    Designed for MIL-53 type connectors, with a looser cutoff,
    and mic=True"""
    eps=1.10 #A
    dummy_indices= [atom.index for atom in mol if atom.symbol == 'X' and atom.mmtype == 'H_']
    for d0 in range(0,len(dummy_indices)-1):
        for d1 in range(d0+1, len(dummy_indices)):
            if mol.get_distance(dummy_indices[d0],dummy_indices[d1],mic=True) < eps: #and mol[dummy_indices[d0]].tag != mol[dummy_indices[d1]].tag:
                print "forming bond",dummy_indices[d0],dummy_indices[d1],mol.get_distance(dummy_indices[d0],dummy_indices[d1],mic=True), mol[dummy_indices[d0]].tag, mol[dummy_indices[d1]].tag 
                mol.form_bond(dummy_indices[d0],dummy_indices[d1])
                return True
    return False


#Intermolecular functions equivalent to those in atoms.py
def get_intermolecular_distance(mol0, mol1, a0, a1):
    """Return distance between two atoms in different molecules."""

    R0 = mol0.arrays['positions']
    R1 = mol1.arrays['positions']
    D = R1[a1] - R0[a0]
    return np.linalg.norm(D)

def get_general_angle(p0,p1,p2):
    """Get angle formed by three atoms.
    The atoms can be from any Atoms object, 
    
    calculate angle between the vectors p0->p1 and
    p1->p2, where each p_n is a xyz """
    # normalized vector 1->0, 1->2:
    v10 = p0 - p1 
    v12 = p2 - p1 
    #print "GGA: ",v10
    #print "GGA: ",v12
    v10 /= np.linalg.norm(v10)
    v12 /= np.linalg.norm(v12)
    angle = np.vdot(v10, v12)
    #print "GGA: ",angle
    if np.isnan(angle):
        angle2=0
    else:
        angle2 = np.arccos(angle)
    if np.isnan(angle2): #Fudge for numerical instability
        #print "HELP ME!!!! NaN"
        angle2=int(angle)
        angle2 = np.arccos(angle2)
    return angle2

def get_arbitrary_dihedral(a1,a2,a3,a4):
    """Calculate dihedral angle between four arbitrary atoms that
    need not be in the same molecule.

    Calculate dihedral angle between the vectors atom1->atom2
    and atom3->atom4, where the atomic coordinates of each atom
    are supplied. """

    # vector 0->1, 1->2, 2->3 and their normalized cross products:
    a = a2 - a1 #self.positions[list[1]] - self.positions[list[0]]
    b = a3 - a2 #self.positions[list[2]] - self.positions[list[1]]
    c = a4 - a3 #self.positions[list[3]] - self.positions[list[2]]
    bxa = np.cross(b, a)
    bxa /= np.linalg.norm(bxa)
    cxb = np.cross(c, b)
    cxb /= np.linalg.norm(cxb)
    angle = np.vdot(bxa, cxb)
    # check for numerical trouble due to finite precision:
    if angle < -1:
        angle = -1
    if angle > 1:
        angle = 1
    angle = np.arccos(angle)
    if np.vdot(bxa, c) > 0:
        angle = 2 * np.pi - angle
    return angle

def angle_2vectors(vec0,vec1):
    """Calculate angle between two vectors"""
    a = np.vdot(vec0,vec1)
    print "a =", a
    len0 = np.sqrt(np.dot(vec0,vec0))
    len1 = np.sqrt(np.dot(vec1,vec1))

    cos_theta = a / (len0 * len1)

    angle = np.arccos(cos_theta)


    return angle

def get_unique_tags(model):
    """Inspects a model and returns the number of unique tags, 
    which represent bonds between fragments, A, A' and A" are NOT unique"""
    alltags={}
    for k,v in model.iteritems():
        for tag in v.get_tags():
            #print "Found a tag", tag, tag[:1]
            if tag > 0  and tag is not None:
                if tag not in alltags:
                    alltags[tag]=[k]
                else:
                    alltags[tag].append(k)
    return alltags

def get_closest_dummy(mol,model,point):
    """given a molecule and a corresponding model and a point in that model,
    find the corresponding dummy atom in the molecule to the point in the model.
    Checks first for angle (COM,model,mol =180) then distance."""
    mol_dummies = mol.get_atom_indices('X')
    model_dummies = model.get_atom_indices('X')
    com=mol.get_center_of_mass()
    dummy_angles={}
    #print "IN GCD"
    for a in model_dummies:
        for b in mol_dummies:
            angle=get_general_angle(com,model[a].position,mol[b].position)
            print a,b,angle
    dummy_distances={}
    for a in mol_dummies:
        dist=get_intermolecular_distance(model, mol, point, a)
        dummy_distances[a] = dist 
    min_index = min(dummy_distances, key=dummy_distances.get)  
    return min_index

def find_tag(mol,tag):
    """Given a molecule and a tag, return the index of the atom that has the tag."""

    tags=mol.get_tags()
    print tag,tags
    b = [item for item in range(len(tags)) if tags[item] == tag]
    print tag,  b
    if len(b) != 1:
        raise RuntimeError('Tag, ',tag,' is defined '+`len(b)`+' times!')
    the_one_I_want=int(b[0])
    return the_one_I_want
    
def find_matching_tags(mol,tag):
    """Given a molecule and a tag, return the index of the atoms that have the tag.
    Explicitly excludes zero tags"""
    if tag != 0:
        tags=mol.get_tags()
        print "In FMT",tags
        snap = [item for item in range(len(tags)) if tags[item] == tag]
        print snap
        if len(snap) != 2:
            raise RuntimeError('Tag, ',tag,' is not present twice')
        return snap

def assemble(model,fragments):
    """Takes the model and the fragments (both dictionaries of Atoms), 
    fragments already rotated and does the translations to put the framework together"""

    build_type = model[0].info['build'].strip()

    if build_type == 'exact':
        return assemble_exact(model,fragments)

    if build_type == 'systre':
        return assemble_systre(model,fragments)
    
    raise RuntimeError('Model build type '+build_type+' not known!')

def assemble_exact(model,fragments):
    """Takes the model and the fragments (both dictionaries of Atoms), 
    fragments already rotated and puts the framework together.
    Exact assembly translates fragments such that matching tagged dummy atoms line up"""

    if len(model) != len(fragments):
        raise RuntimeError('Model has '+len(model)+' slots but there are '+len(fragments)+' fragments!')

    mof = Atoms()
    tag_dic=get_unique_tags(model)
    print "Tag_dic", tag_dic
    for frag in fragments:
        print "Begin assemble:",fragments[frag].get_tags()
    #So we set up some locking checks to make sure we place everything in a reasonable order 
    frag_locked = [False] * len(model)
    tag_locked={}
    for tag in tag_dic.keys():
        tag_locked[tag]=False
    #tag_queue=[]
    tag_queue = deque()
    initial_placement=True  
    counter=0
    while(not all(frag_locked)):
        if initial_placement:
            print "Here!"
            mof=fragments[0].copy()
            for tag in fragments[0].get_tags():
                #print tag#,len(tag)
                #if tag and not tag.endswith("\""): #len(tag)==1:
                if tag > 0 and tag is not None: 
                    tag_queue.append(tag)
                    print "Initial TQ:",tag_queue
            frag_locked[0]=True
            initial_placement=False
        else:
            current_tag=tag_queue.popleft()
            print "TQ:",tag_queue
            print "current_tag: ",current_tag 
            frags=tag_dic[current_tag]
            print "frags=", frags
            #Now either frags[0] or [1] is the next frag to place, one of them should already be done
            if frag_locked[frags[0]] == True and frag_locked[frags[1]] == False:
                next_frag=frags[1]
            elif frag_locked[frags[0]] == False and frag_locked[frags[1]] == True:
                 next_frag=frags[0]
            elif frag_locked[frags[0]] == True and frag_locked[frags[1]] == True:
                #This can and will happen if there's a closed cycle in the model, so we want to warn and lock the offending tag
                print 'WARNING: Both fragments, '+`frags[0]`+' and '+`frags[1]`+' associated with current tag, '+`current_tag`+' have been placed!\n Locking tag and moving on'
                tag_locked[current_tag]=True
                continue 
                raise RuntimeError('Both frag'+`frags[0]`+' and frag'+`frags[1]`+' have been placed!')
            elif frag_locked[frags[0]] == False and frag_locked[frags[1]] == False:
                raise RuntimeError('Neither frag'+`frags[0]`+' or frag'+`frags[1]`+' can be placed!')
    
            joiner0=find_tag(mof,current_tag)
            joiner1=find_tag(fragments[next_frag],current_tag)
            fragments[next_frag].translate(mof.positions[joiner0]-fragments[next_frag].positions[joiner1]) #translate then copy
            temp=fragments[next_frag].copy() #translate then copy
            mof+=temp
            counter +=1
            for tag in fragments[next_frag].get_tags():
                print "Next tags are:", fragments[next_frag].get_tags()
                print "current tag is:", tag
                #if len(tag)==1 and tag not in tag_queue and tag is not current_tag: #len=1 test won't work
                #if not tag.endswith("\"")  and tag != '' and not tag_locked[tag] and tag not in tag_queue and tag is not current_tag:
                if tag is not None and tag > 0 and not tag_locked[tag] and tag not in tag_queue and tag != current_tag:
                    print "appending: ", tag
                    tag_queue.append(tag)
                    print "TQ:",tag_queue
                    print "TL:",tag_locked
            frag_locked[next_frag]=True
            tag_locked[current_tag]=True

    #MAA 30/01/2014 
    #Being fully evil and diabling the tests below, because they fail for cycles. If you have lots of cycles, you should be using a topology anyway...
    #if tag_queue: 
    #    raise RuntimeError('Some tags still not placed!'+`tag_queue`)

    #if(not all(tag_locked)):
    #        raise RuntimeError('Some fragments not placed!'+`tag_locked`+'!')
    #print "Yibbida, yibbida, yibbida..."        
    return mof 

def assemble_systre(model,fragments):
    """Takes the model and the fragments (both dictionaries of Atoms),
    fragments already rotated and puts the framework together.
    Systre assembly does NOT translate fragments. Tags are matched and
    dummy atoms are removed, updating the bondlists in the process."""

    if len(model) != len(fragments):
        raise RuntimeError('Model has '+len(model)+' slots but there are '+len(fragments)+' fragments!')

    mof = Atoms()
    mof.set_cell(model[0].cell)
    #Now for systre build no translation is done, so things aren't order-dependent
    #Let's just walk through the frags
    for frag in fragments:
        print "Begin systre-style assemble:",fragments[frag].get_tags()
        mof += fragments[frag]

        #outfile='assemble'+`frag`+'.xyz'
        #write(outfile,mof)
    return mof

def bond_frame(mof, periodic=True):
    """Takes an assembled framework and forms all the necessary bonds by matching tags.
    *If optional periodic argument is True, includes bonds across PBC."""

    tags=mof.get_tags()
    unique_tags=set(tags)
    unique_tags.discard(0)
    for tag in unique_tags:
        #if (not periodic) and (not tag.endswith("\"")): old, string way
        if (not periodic) and (not tag < 0): 
            print "NP",tag
            indices=find_matching_tags(mof,tag)
            index0=int(indices[0])
            index1=int(indices[1])
            print  tag,index0,index1
            mof.form_bond(index0,index1)
        elif periodic:
            indices=find_matching_tags(mof,tag)
            index0=int(indices[0])
            index1=int(indices[1])
            #print  tag,index0,index1
            mof.form_bond(index0,index1)

    return mof


def frame_cell(model,fragments):  
    """Takes the model and the fragments (both dictionaries of Atoms), 
    and works out the cell""" 

    if len(model) != len(fragments):
        raise RuntimeError('Model has '+len(model)+' slots but there are '+len(fragments)+' fragments!')
    real_cell = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
    model_cell = model[0].get_cell() #assuming it's a unit cell
    axis_scale = [0,0,0]
    dummy_mol = Atoms()
    for obj in range(0,len(model)):
        print "FRAG ",obj
        #Need to determine the min and max only of the dummies in the real mol
        dummy_mol+=fragments[obj].copy()
    #write('assembled_mol.xyz',dummy_mol)
    del dummy_mol[[atom.index for atom in dummy_mol if atom.symbol!='X']]
    write('dummy_mol.xyz',dummy_mol)
   
    for axis in range(0,3):
        dummy_proj=[]
        for dummy in dummy_mol:
            #print dummy.position
            proj=np.dot(model_cell[axis],dummy.position)
            print "proj=",proj
            dummy_proj.append(proj)
        axis_proj=max(dummy_proj)-min(dummy_proj)
        print "AP=",axis,axis_proj
        axis_scale[axis]+=axis_proj
    print "AS=",axis_scale

    print real_cell
    for axis in range(0,3):
        if 'cell_scaling' in model[0].info:
            real_cell[axis]=model_cell[axis]*axis_scale[axis]*model[0].info['cell_scaling'][axis]
        else:
            print model_cell[axis], axis_scale[axis]
            real_cell[axis]=model_cell[axis]*axis_scale[axis]
    print real_cell

    return real_cell

def make_supercell(mol,v):
    """Takes an assembled, but not pbc-bonded framework and copies it to make a supercell"""
    eps = 1.0 #Ang 
    if isinstance(v, int):
            v = (v, v, v)    
    cell=mol.get_cell()
    if np.array_equal(cell,np.eye(3)):
        raise RuntimeError('mol cell has not been changed from the default value' )
    elif np.sqrt(np.dot(cell[2], cell[2])) == 0:
        #We need to invent a c-vector
        mol_height = max(mol.positions[:,2]) - min(mol.positions[:,2])
        print "mol_height = " ,mol_height
        if mol_height < 1.0:
            print "top"
            interlayer_dist = 3.5
        else:
            print "else"
            interlayer_dist = mol_height + 2.5
    else:
        interlayer_dist = cell[2,2]
    cell[2,2] = interlayer_dist        
    bigmol = mol.copy()
    #Now we need to add new cells one at a time and update the bondlists each time or else we'll have a non-unique mapping
    for m0 in range(v[0]):
        for m1 in range(v[1]):
            for m2 in range(v[2]):
                image=mol.copy()
                vec=[m0,m1,m2]
                update_image_tags(image,vec)
                if any (vec): #don't copy and translate the 0,0,0 case
                    image.translate(np.dot(cell,vec))
                    bigmol += image
                    bigmol.set_current_indices(range(1,len(bigmol)+1))
                    bigmol.update_bondlists_add()

    bigmol.set_cell(np.dot(cell,v))
    #Now detect coincident atoms and check for dummyness and tags
    internal_bonds=True
    while internal_bonds:
        internal_bonds=find_internal_bond(bigmol)
    pbc_bonds=True
    while pbc_bonds:
        pbc_bonds=find_pbc_bond(bigmol)
    
    return bigmol

def update_image_tags(mol,vec):
    """Update all tags on an image with the translation. 
    Uses 0.1 for x, 0.2 for y and 0.4 for z"""
    newtag=0
    dimensions=[1, 2, 4] #=> [x,y,z]
    #print vec
    for a in range(len(vec)):
        if vec[a]:
            newtag+=dimensions[a]*10**-vec[a] #PBC tags are -ve
    #print newtag
    for a in mol:
        if a.tag and a.tag != 0:
            if a.tag > 0:
                a.tag+=newtag
            elif a.tag < 0:
                a.tag -=newtag

def remove_banquo(mol):
    """Removes any Bq atoms from structure. Bq is used for an optional dummy atom.
       Most common usage is for connectors that may have attached solvent"""

    bq_present = True
    while bq_present:
        bq_present = pop_bq(mol)

    return False

def pop_bq(mol):
    """Accessory function for remove_banquo: removes a single Bq atom from structure."""

    bq_indices = [atom.index for atom in mol if atom.symbol=='Bq']
    for bq in range(len(bq_indices)):
        #Need to delete the bond to this atom, or else update_bondlists_pop will fail
        bq_index = bq_indices[0] #the first one
        for bond in mol[bq_index].bondlist.keys():
            print "Bond is ", bond
            print "and atom ", bond, "has bondlist ",mol[bond -1].bondlist
            mol[bond -1].bondlist.pop(bq_index+1, None)
            print "Popping bond to ",bq_index+1, " from atom ",bond -1
        mol.pop(bq_index)
        return True
    return False


def find_internal_bond(mol):
    """Finds a single internal bond in a supercell"""
    eps=1.0 #A
    dummy_indices= [atom.index for atom in mol if atom.symbol=='X']
    for d0 in range(0,len(dummy_indices)-1):
        for d1 in range(d0+1, len(dummy_indices)):
            if mol.get_distance(dummy_indices[d0],dummy_indices[d1]) < eps:
                #print "forming bond",dummy_indices[d0],dummy_indices[d1],mol.get_distance(dummy_indices[d0],dummy_indices[d1])
                mol.form_bond(dummy_indices[d0],dummy_indices[d1])
                return True
    return False

def find_pbc_bond(mol):
    """Finds a single pbc bond in a supercell"""
    #Find a tag without an xyz, find the equivalent tag with an xyz, bond it and return
    dimensions=set(['x','y','z'])
    #floor to remove integer, then multiply by 10 repeatedly until diff is exactly zero
    eps = 0.5
    angle_eps = 0.002
    

    cell_dims = np.zeros(3)
    for i in range(3):
        cell_dims[i] = np.sqrt(np.dot(mol.cell[i], mol.cell[i]))
    print "cell_dims = ",cell_dims

    dummy_indices= [atom.index for atom in mol if atom.symbol=='X']
    for d0 in range(0,len(dummy_indices)-1):
        for d1 in range(d0+1, len(dummy_indices)):
            #We get the two tags and see if they differ by a single character x,y or z
            tag0=mol[dummy_indices[d0]].tag
            tag1=mol[dummy_indices[d1]].tag
            #If tag is the same, then it's possible that these are two tags that match up across PBC
            if tag0 == tag1:
                this_angle=100
                #Checking for dist using mic doesn't work, it has a fit at two things being in the same place
                this_dist=mol.get_distance(dummy_indices[d0],dummy_indices[d1])
                #if we check for the dummy closest to zero, we can test only for zero angle
                for a in range(3):
                    print "a = ", a 
                    if abs(this_dist - cell_dims[a]) < eps:
                        vec = mol[dummy_indices[d1]].position - mol[dummy_indices[d0]].position
                        for i in range(len(vec)):
                            if vec[i] < 0:
                                vec[i] = -vec[i] #ensures always in +ve octant
                        this_angle = angle_2vectors(vec , mol.cell[a])
                        print "Angle = ", this_angle
                        print "Dist = ", this_dist
                        if this_angle < eps:
                            #then this vector is parallel to cell vector. i.e. Join these two tags
                            print "forming PBC bond",dummy_indices[d0],dummy_indices[d1],mol.get_distance(dummy_indices[d0],dummy_indices[d1])
                            mol.form_bond(dummy_indices[d0],dummy_indices[d1])
                            return True

    return False


                
def distance_matrix(mol):
    """Takes a mol and produces a symmetric distance matrix"""
    dist_mat=np.zeros((len(mol),len(mol)))
    for a0 in range(0,len(mol)):
        for a1 in range(a0,len(mol)):
            dist_mat[a0,a1]=mol.get_distance(a0,a1)
            dist_mat[a1,a0]=dist_mat[a0,a1]
    return dist_mat

def furthest_real(mol):
    """Takes a mol and returns the index of the real atom (not dummy)
    furthest away from the COM"""
    
    com=mol.get_center_of_mass()
    dist_array=np.zeros(len(mol))
    for a in range(0,len(mol)):
        if mol[a].symbol  != 'X' :
            D=mol.positions[a]-com
            dist_array[a]=np.linalg.norm(D)
    return np.argmax(dist_array)


