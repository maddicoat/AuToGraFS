from ase.atom import Atom
from ase.atoms import Atoms
from ase.io import read, write
import numpy as np

def make_model(topology,basepath,mol_types,mol_ids):
    """Makes an appropriately sized model file using the given topology and
    molecules (things) supplied to it. Doesn't check if you have the right number
    of things or that they're of the correct shape as this check is built into AuToGraFS
    
    Known topologies:

    ==========  ======================  =============
    topology    center                  linker
    ==========  ======================  =============
    srs         triangle                linear
    dia         tetrahedral             linear
    sra         tetrahedral             linear
    nbo         square                  linear
    cds         square                  linear
    bnn         trigonal_bipyramid      linear
    pcu         octahedral              linear
    pts         tetrahedral             rectangle 
    ptt         tetrahedral             rectangle [rectangle] 
    bor         tetrahedral             triangle
    ctn         tetrahedral             triangle
    pth         tetrahedral             rectangle
    pto         square                  triangle 
    pyr         octahedral              triangle 
    stp         tri_prism               rectangle    
    soc         tri_prism               rectangle    
    tbo         square                  triangle 
    spn         octahedral              triangle
    gar         octahedral              tetrahedral
    ibd         octahedral              tetrahedral
    iac         octahedral              tetrahedral
    ifi         octahedral              tetrahedral #issue with octahedral actually being tri_prism, also RCSR and Acta Cryst. (2006). A62, 350-355 disagree on coordinates
    mtn_e       tri_prism               linear    #MIL-101 , actual symbol is mtn-e
    mfu4        mfu4                    rectangle #special centre type for mfu4 
    mil53       mil53                   linear #special centre type for mil53 
    rtl         octahedral              triangle 
    sod         tetrahedral             linear
    sqc19       icosahedral             linear # epinet naming
    pyr_c       octahedral              triangle #actual symbol is pyr-c , interpenetrated
    qom         octahedral              triangle #
    rhr         square                  linear 
    ntt         rectangle (square)      triangle [triangle] #actually the two triangles form a hexagon
    ntt_46      rectangle (square)      hexagon #pre-filled version of above
    eta_c       triangle                linear #interpenetrated, actual symbol is eta-c
    eta_c3      triangle                linear #interpenetrated, actual symbol is eta-c3
    eta_c4      triangle                linear #interpenetrated, actual symbol is eta-c4
    lig_c       triangle                linear #interpenetrated, actual symbol is lig-c
    srs_a_c     triangle                linear #interpenetrated, actual symbol uses hyphens, you should have learnt this already ;-)
    srs_c       triangle                linear #interpenetrated, actual symbol uses hyphens
    srs_c3      triangle                linear #interpenetrated, actual symbol uses hyphens
    srs_c4      triangle                linear #interpenetrated, actual symbol uses hyphens
    srs_c8      triangle                linear #interpenetrated, actual symbol uses hyphens
    srs_c4s     triangle                linear #interpenetrated, actual symbol is srs-c4*
    srs_c2s     triangle                linear #interpenetrated, actual symbol is srs-c2*
    ths_c       triangle                linear #interpenetrated
    twt_c       triangle                linear #interpenetrated 
    twt_c3      triangle                linear #interpenetrated 
    chs1        rectangle               triangle+linear
    ==========  ======================  =============
    Td          triangle                triangle
    ==========  ======================  =============
   
    sqp_B1      tetrahedral             
    sqp_B2      tetrahedral             tetrahedral
    dia_B2      tetrahedral             tetrahedral 

"""
    print basepath

    #Now set up some default locations for things.
    #basepath="/home/maddicoat/src/ase/ase/database/"
    if not basepath.endswith("/"):
        basepath = basepath+'/'
    centerpath="centers/"
    linkerpath="linkers/"
    pillarpath="pillars/"
    fgrppath="functional_groups/"
    
    #read "things" into a list of molecules
    molecules = []
    mol_sizes = []
    for mol_type,mol_id in zip(mol_types,mol_ids):
        filename = mol_id+'.inp'
        try:
            molecules.append(read(filename))
            molecules[-1].info['name'] = mol_id
            mol_sizes.append(furthest_dummy(molecules[-1]))
        except IOError:
            path=basepath+mol_type+'s/'+filename
            molecules.append(read(path))
            molecules[-1].info['name'] = mol_id
            mol_sizes.append(furthest_dummy(molecules[-1]))
   
    #print mol_sizes
    #now feed the molecules off to the topology
    if topology == 'srs':
        from ase.topology.srs import make_srs
        return make_srs(mol_ids,mol_sizes)

    if topology == 'dia':
        from ase.topology.dia import make_dia
        return make_dia(mol_ids,mol_sizes)

    if topology == 'sra':
        from ase.topology.sra import make_sra
        return make_sra(mol_ids,mol_sizes)

    if topology == 'nbo':
        from ase.topology.nbo import make_nbo
        return make_nbo(mol_ids,mol_sizes)
    
    if topology == 'cds':
        from ase.topology.cds import make_cds
        return make_cds(mol_ids,mol_sizes)

    if topology == 'pts':
        from ase.topology.pts import make_pts
        return make_pts(mol_ids,mol_sizes)

    if topology == 'ptt':
        from ase.topology.ptt import make_ptt
        return make_ptt(mol_ids,mol_sizes)

    if topology == 'bnn':
        from ase.topology.bnn import make_bnn
        return make_bnn(mol_ids,mol_sizes)

    if topology == 'pcu':
        from ase.topology.pcu import make_pcu
        return make_pcu(mol_ids,mol_sizes)

    if topology == 'bor':
        from ase.topology.bor import make_bor
        return make_bor(mol_ids,mol_sizes)

    if topology == 'ctn':
        from ase.topology.ctn import make_ctn
        return make_ctn(mol_ids,mol_sizes)

    if topology == 'pth':
        from ase.topology.pth import make_pth
        return make_pth(mol_ids,mol_sizes)

    if topology == 'pto':
        from ase.topology.pto import make_pto
        return make_pto(mol_ids,mol_sizes)
    
    if topology == 'pyr':
        from ase.topology.pyr import make_pyr
        return make_pyr(mol_ids,mol_sizes)

    if topology == 'stp':
        from ase.topology.stp import make_stp
        return make_stp(mol_ids,mol_sizes)

    if topology == 'soc':
        from ase.topology.soc import make_soc
        return make_soc(mol_ids,mol_sizes)

    if topology == 'tbo':
        from ase.topology.tbo import make_tbo
        return make_tbo(mol_ids,mol_sizes)

    if topology == 'spn':
        from ase.topology.spn import make_spn
        return make_spn(mol_ids,mol_sizes)

    if topology == 'gar':
        from ase.topology.gar import make_gar
        return make_gar(mol_ids,mol_sizes)
    
    if topology == 'ibd':
        from ase.topology.ibd import make_ibd
        return make_ibd(mol_ids,mol_sizes)

    if topology == 'iac':
        from ase.topology.iac import make_iac
        return make_iac(mol_ids,mol_sizes)

    if topology == 'mtn_e':
        from ase.topology.mtn_e import make_mtn_e
        return make_mtn_e(mol_ids,mol_sizes)

    if topology == 'mfu4':
        from ase.topology.mfu4 import make_mfu4
        return make_mfu4(mol_ids,mol_sizes)

    if topology == 'mil53':
        from ase.topology.mil53 import make_mil53
        return make_mil53(mol_ids,mol_sizes)

    if topology == 'rtl':
        from ase.topology.rtl import make_rtl
        return make_rtl(mol_ids,mol_sizes)

    if topology == 'sod':
        from ase.topology.sod import make_sod
        return make_sod(mol_ids,mol_sizes)

    #ifi fails due to tri_prism/octahedron issue. Worth solving?
    #if topology == 'ifi':
    #    from ase.topology.ifi import make_ifi
    #    return make_ifi(mol_ids,mol_sizes)

    if topology == 'sqc19':
        from ase.topology.sqc19 import make_sqc19
        return make_sqc19(mol_ids,mol_sizes)

    if topology == 'pyr_c':
        from ase.topology.pyr_c import make_pyr_c
        return make_pyr_c(mol_ids,mol_sizes)

    if topology == 'qom':
        from ase.topology.qom import make_qom
        return make_qom(mol_ids,mol_sizes)

    if topology == 'rhr':
        from ase.topology.rhr import make_rhr
        return make_rhr(mol_ids,mol_sizes)
    
    if topology == 'ntt':
        from ase.topology.ntt import make_ntt
        return make_ntt(mol_ids,mol_sizes)

    if topology == 'ntt_46':
        from ase.topology.ntt_46 import make_ntt_46
        return make_ntt_46(mol_ids,mol_sizes)

    if topology == 'eta_c3':
        from ase.topology.eta_c3 import make_eta_c3
        return make_eta_c3(mol_ids,mol_sizes)

    if topology == 'eta_c4':
        from ase.topology.eta_c4 import make_eta_c4
        return make_eta_c4(mol_ids,mol_sizes)

    if topology == 'eta_c':
        from ase.topology.eta_c import make_eta_c
        return make_eta_c(mol_ids,mol_sizes)
    
    if topology == 'lig_c':
        from ase.topology.lig_c import make_lig_c
        return make_lig_c(mol_ids,mol_sizes)

    if topology == 'srs_a_c':
        from ase.topology.srs_a_c import make_srs_a_c
        return make_srs_a_c(mol_ids,mol_sizes)
    
    if topology == 'srs_c':
        from ase.topology.srs_c import make_srs_c
        return make_srs_c(mol_ids,mol_sizes)

    if topology == 'srs_c3':
        from ase.topology.srs_c3 import make_srs_c3
        return make_srs_c3(mol_ids,mol_sizes)

    if topology == 'srs_c4':
        from ase.topology.srs_c4 import make_srs_c4
        return make_srs_c4(mol_ids,mol_sizes)

    if topology == 'srs_c8':
        from ase.topology.srs_c8 import make_srs_c8
        return make_srs_c8(mol_ids,mol_sizes)

    if topology == 'srs_c4s':
        from ase.topology.srs_c4s import make_srs_c4s
        return make_srs_c4s(mol_ids,mol_sizes)

    if topology == 'srs_c2s':
        from ase.topology.srs_c2s import make_srs_c2s
        return make_srs_c2s(mol_ids,mol_sizes)
    
    if topology == 'ths_c':
        from ase.topology.ths_c import make_ths_c
        return make_ths_c(mol_ids,mol_sizes)

    if topology == 'twt_c':
        from ase.topology.twt_c import make_twt_c
        return make_twt_c(mol_ids,mol_sizes)

    if topology == 'twt_c3':
        from ase.topology.twt_c3 import make_twt_c3
        return make_twt_c3(mol_ids,mol_sizes)

    if topology == 'chs1':
        from ase.topology.chs1 import make_chs1
        return make_chs1(mol_ids,mol_sizes)

    #-------------------------------------------------------#

    if topology == 'Td':
        from ase.topology.Td import make_Td
        return make_Td(mol_ids,mol_sizes)

    if topology == 'sqp_B1':
        from ase.topology.sqp_B1 import make_sqp_B1
        return make_sqp_B1(mol_ids,mol_sizes)

    if topology == 'sqp_B2':
        from ase.topology.sqp_B2 import make_sqp_B2
        return make_sqp_B2(mol_ids,mol_sizes)
    
    if topology == 'dia_B2':
        from ase.topology.dia_B2 import make_dia_B2
        return make_dia_B2(mol_ids,mol_sizes)



    raise RuntimeError('Topology: '+topology+' not implemented!')



def furthest_dummy(mol):
    """Takes a mol and returns the distance of the dummy
    furthest away from the COM"""

    com=mol.get_center_of_mass()
    dist_array=np.zeros(len(mol))
    for a in range(0,len(mol)):
        if mol[a].symbol  == 'X' :
            D=mol.positions[a]-com
            dist_array[a]=np.linalg.norm(D)
    return np.amax(dist_array)
