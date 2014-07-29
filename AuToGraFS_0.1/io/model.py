from ase.atoms import Atoms
from ase.parallel import paropen
import numpy as np

#This is pathetic, need to add in proper checking and parsing later
def read_model(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    scaled_cell = False    
    build_set = False    
    extra_dummies = False    
    lines = fileobj.readlines()
    L1 = lines[0].split("=")
    nobjects=int(L1[-1])
    del lines[:1]   
    if "build" in lines[0].lower():
        L1 = lines[0].split("=")
        build_type = L1[-1].strip().lower()
        build_set = True
        del lines[:1]   
    if "extra" in lines[0].lower():
        #print "Found the dummy line"
        extra_dummies = True
        del lines[:1]   
    if "cell" in lines[0].lower():
        x1,y1,z1=lines[1].split()
        x2,y2,z2=lines[2].split()
        x3,y3,z3=lines[3].split()
        del lines[:5]
    if "scale" in lines[0].lower():
        scaled_cell = True
        a_scale = int(lines[1])
        b_scale = int(lines[2])
        c_scale = int(lines[3])
        del lines[:5]
    model_molecule={}
    for obj in xrange(0, nobjects):
    #for obj in xrange(1, nobjects+1):
        info={}
        #name='sbu'
        #name+=`obj`
        name = obj
        #print name
        #Get  obj_id, natoms and geometry tag
        L1 = lines[0].split(":")
        obj_id = str(L1[0])
        obj_id.rstrip()
        L2 = lines[1].split()
        if len(L2) == 1:
            natoms=int(L2[0])
        else:
            print "Line should only contain the number of atoms"
            print lines
            exit()
        L3 = lines[2].split("=")
        geometry = L3[-1].rstrip()
        #print obj_id,natoms,geometry
        info={'SBUtype':obj_id,'shape':geometry}
        del lines[:3]
        images = []
        positions = []
        symbols = []
        tags = []
        for coord_line in range(natoms):
            tmp = lines[coord_line].split()
            if len(tmp)==5:
                symbol, x, y, z, tag = lines[coord_line].split()[:5]
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
                tags.append(int(tag))
            if len(tmp)==4:
                symbol, x, y, z = lines[coord_line].split()[:4]
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
                tags.append(None)
        print tags
        fragment = Atoms(symbols=symbols, positions=positions,tags=tags,info=info,fragmentIDs=obj)
        if obj ==0:
            fragment.set_cell([(float(x1), float(y1), float(z1)), (float(x2), float(y2), float(z2)), (float(x3), float(y3), float(z3))])
            if scaled_cell:
                fragment.info['cell_scaling'] = [a_scale,b_scale,c_scale]
            if extra_dummies:
                fragment.info['extra_dummies'] = True 
            if build_set:
                fragment.info['build'] = build_type
            else:
                print "Model doesn't have build type set. Shouls be 'exact' or 'systre'"
                exit()
        model_molecule[name]=fragment
        #print lines[:natoms]
        del lines[:natoms+1]
        #if not lines[0].strip():
        #    print lines[0]
        #    print "here"
        #    break
    for slot in model_molecule:
        print model_molecule[slot]
        print model_molecule[slot].get_tags()
    return model_molecule

def get_unique_tags(model):
    """Inspects a model and returns the number of unique tags, 
    which represent bonds between fragments, A, A' and A" are NOT unique"""
    alltags=[]
    for k,v in model.iteritems():
        alltags.append(v.get_tags())
        print alltags
    #    print blah
    #    print "======"
