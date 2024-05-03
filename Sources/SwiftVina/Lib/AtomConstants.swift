//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/11/23.
//

import Foundation

// based on SY_TYPE_* but includes H
let EL_TYPE_H: sz     =  0
let EL_TYPE_C: sz     =  1
let EL_TYPE_N: sz     =  2
let EL_TYPE_O: sz     =  3
let EL_TYPE_S: sz     =  4
let EL_TYPE_P: sz     =  5
let EL_TYPE_F: sz     =  6
let EL_TYPE_Cl: sz    =  7
let EL_TYPE_Br: sz    =  8
let EL_TYPE_I: sz     =  9
let EL_TYPE_Si: sz    = 10 // Silicon
let EL_TYPE_At: sz    = 11 // Astatine
let EL_TYPE_Met: sz   = 12
let EL_TYPE_Dummy: sz = 13
let EL_TYPE_SIZE: sz  = 14

// AutoDock4
let AD_TYPE_C: sz     =  0
let AD_TYPE_A: sz     =  1
let AD_TYPE_N: sz     =  2
let AD_TYPE_O: sz     =  3
let AD_TYPE_P: sz     =  4
let AD_TYPE_S: sz     =  5
let AD_TYPE_H: sz     =  6 // non-polar hydrogen
let AD_TYPE_F: sz     =  7
let AD_TYPE_I: sz     =  8
let AD_TYPE_NA: sz    =  9
let AD_TYPE_OA: sz    = 10
let AD_TYPE_SA: sz    = 11
let AD_TYPE_HD: sz    = 12
let AD_TYPE_Mg: sz    = 13
let AD_TYPE_Mn: sz    = 14
let AD_TYPE_Zn: sz    = 15
let AD_TYPE_Ca: sz    = 16
let AD_TYPE_Fe: sz    = 17
let AD_TYPE_Cl: sz    = 18
let AD_TYPE_Br: sz    = 19
let AD_TYPE_Si: sz    = 20 // Silicon
let AD_TYPE_At: sz    = 21 // Astatine
let AD_TYPE_G0: sz    = 22 // closure of cyclic molecules
let AD_TYPE_G1: sz    = 23
let AD_TYPE_G2: sz    = 24
let AD_TYPE_G3: sz    = 25
let AD_TYPE_CG0: sz   = 26
let AD_TYPE_CG1: sz   = 27
let AD_TYPE_CG2: sz   = 28
let AD_TYPE_CG3: sz   = 29
let AD_TYPE_W: sz     = 30 // hydrated ligand
public let AD_TYPE_SIZE: sz  = 31

// X-Score
public let XS_TYPE_C_H: sz   =  0
public let XS_TYPE_C_P: sz   =  1
public let XS_TYPE_N_P: sz   =  2
public let XS_TYPE_N_D: sz   =  3
public let XS_TYPE_N_A: sz   =  4
public let XS_TYPE_N_DA: sz  =  5
public let XS_TYPE_O_P: sz   =  6
public let XS_TYPE_O_D: sz   =  7
public let XS_TYPE_O_A: sz   =  8
public let XS_TYPE_O_DA: sz  =  9
public let XS_TYPE_S_P: sz   = 10
public let XS_TYPE_P_P: sz   = 11
public let XS_TYPE_F_H: sz   = 12
public let XS_TYPE_Cl_H: sz  = 13
public let XS_TYPE_Br_H: sz  = 14
public let XS_TYPE_I_H: sz   = 15
public let XS_TYPE_Si: sz    = 16 // Silicon
public let XS_TYPE_At: sz    = 17 // Astatine
public let XS_TYPE_Met_D: sz = 18
public let XS_TYPE_C_H_CG0: sz = 19 // closure of cyclic molecules
public let XS_TYPE_C_P_CG0: sz = 20
public let XS_TYPE_G0: sz      = 21
public let XS_TYPE_C_H_CG1: sz = 22
public let XS_TYPE_C_P_CG1: sz = 23
public let XS_TYPE_G1: sz      = 24
public let XS_TYPE_C_H_CG2: sz = 25
public let XS_TYPE_C_P_CG2: sz = 26
public let XS_TYPE_G2: sz      = 27
public let XS_TYPE_C_H_CG3: sz = 28
public let XS_TYPE_C_P_CG3: sz = 29
public let XS_TYPE_G3: sz      = 30
public let XS_TYPE_W: sz       = 31 // hydrated ligand
public let XS_TYPE_SIZE: sz    = 32

// DrugScore-CSD
let SY_TYPE_C_3: sz   =  0
let SY_TYPE_C_2: sz   =  1
let SY_TYPE_C_ar: sz  =  2
let SY_TYPE_C_cat: sz =  3
let SY_TYPE_N_3: sz   =  4
let SY_TYPE_N_ar: sz  =  5
let SY_TYPE_N_am: sz  =  6
let SY_TYPE_N_pl3: sz =  7
let SY_TYPE_O_3: sz   =  8
let SY_TYPE_O_2: sz   =  9
let SY_TYPE_O_co2: sz = 10
let SY_TYPE_S: sz     = 11
let SY_TYPE_P: sz     = 12
let SY_TYPE_F: sz     = 13
let SY_TYPE_Cl: sz    = 14
let SY_TYPE_Br: sz    = 15
let SY_TYPE_I: sz     = 16
let SY_TYPE_Met: sz   = 17
let SY_TYPE_SIZE: sz  = 18

public struct AtomKind {
    let name: String
    public let radius: fl
    public let depth: fl
    public let hb_depth: fl // pair (i,j) is HB if hb_depth[i]*hb_depth[j] < 0
    public let hb_radius: fl
    let solvation: fl
    let volume: fl
    let covalent_radius: fl // from http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
}


// generated from edited AD4_parameters.data using a script, 
// then covalent radius added from en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
public let AtomKindData: [AtomKind] = [
    //       name,      radius,          depth,          hb_depth,      hb_r,           solvation,           volume,           covalent radius
    AtomKind(name: "C", radius: 2.00000, depth: 0.15000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77), //  0
    AtomKind(name: "A", radius: 2.00000, depth: 0.15000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00052, volume: 33.51030, covalent_radius: 0.77), //  1
    AtomKind(name: "N", radius: 1.75000, depth: 0.16000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00162, volume: 22.44930, covalent_radius: 0.75), //  2
    AtomKind(name: "O", radius: 1.60000, depth: 0.20000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00251, volume: 17.15730, covalent_radius: 0.73), //  3
    AtomKind(name: "P", radius: 2.10000, depth: 0.20000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 38.79240, covalent_radius: 1.06), //  4
    AtomKind(name: "S", radius: 2.00000, depth: 0.20000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00214, volume: 33.51030, covalent_radius: 1.02), //  5
    AtomKind(name: "H", radius: 1.00000, depth: 0.02000, hb_depth: 0.0, hb_radius: 0.0, solvation: 0.00051, volume: 0.00000, covalent_radius: 0.37), //  6
    AtomKind(name: "F", radius: 1.54500, depth: 0.08000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 15.44800, covalent_radius: 0.71), //  7
    AtomKind(name: "I", radius: 2.36000, depth: 0.55000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 55.05850, covalent_radius: 1.33), //  8
    AtomKind(name: "NA", radius: 1.75000, depth: 0.16000, hb_depth: -5.0, hb_radius: 1.9, solvation: -0.00162, volume: 22.44930, covalent_radius: 0.75), //  9
    AtomKind(name: "OA", radius: 1.60000, depth: 0.20000, hb_depth: -5.0, hb_radius: 1.9, solvation: -0.00251, volume: 17.15730, covalent_radius: 0.73), // 10
    AtomKind(name: "SA", radius: 2.00000, depth: 0.20000, hb_depth: -1.0, hb_radius: 2.5, solvation: -0.00214, volume: 33.51030, covalent_radius: 1.02), // 11
    AtomKind(name: "HD", radius: 1.00000, depth: 0.02000, hb_depth: 1.0, hb_radius: 0.0, solvation: 0.00051, volume: 0.00000, covalent_radius: 0.37), // 12
    AtomKind(name: "Mg", radius: 0.65000, depth: 0.87500, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 1.56000, covalent_radius: 1.30), // 13
    AtomKind(name: "Mn", radius: 0.65000, depth: 0.87500, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 2.14000, covalent_radius: 1.39), // 14
    AtomKind(name: "Zn", radius: 0.74000, depth: 0.55000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 1.70000, covalent_radius: 1.31), // 15
    AtomKind(name: "Ca", radius: 0.99000, depth: 0.55000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 2.77000, covalent_radius: 1.74), // 16
    AtomKind(name: "Fe", radius: 0.65000, depth: 0.01000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 1.84000, covalent_radius: 1.25), // 17
    AtomKind(name: "Cl", radius: 2.04500, depth: 0.27600, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 35.82350, covalent_radius: 0.99), // 18
    AtomKind(name: "Br", radius: 2.16500, depth: 0.38900, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 42.56610, covalent_radius: 1.14), // 19
    AtomKind(name: "Si", radius: 2.30000, depth: 0.20000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00143, volume: 50.96500, covalent_radius: 1.11), // 20
    AtomKind(name: "At", radius: 2.40000, depth: 0.55000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00110, volume: 57.90580, covalent_radius: 1.44), // 21
    AtomKind(name: "G0", radius: 0.00000, depth: 0.00000, hb_depth: 0.0, hb_radius: 0.0, solvation: 0.00000, volume: 0.00000, covalent_radius: 0.77), // 22
    AtomKind(name: "G1", radius: 0.00000, depth: 0.00000, hb_depth: 0.0, hb_radius: 0.0, solvation: 0.00000, volume: 0.00000, covalent_radius: 0.77), // 23
    AtomKind(name: "G2", radius: 0.00000, depth: 0.00000, hb_depth: 0.0, hb_radius: 0.0, solvation: 0.00000, volume: 0.00000, covalent_radius: 0.77), // 24
    AtomKind(name: "G3", radius: 0.00000, depth: 0.00000, hb_depth: 0.0, hb_radius: 0.0, solvation: 0.00000, volume: 0.00000, covalent_radius: 0.77), // 25
    AtomKind(name: "CG0", radius: 2.00000, depth: 0.15000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77), // 26
    AtomKind(name: "CG1", radius: 2.00000, depth: 0.15000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77), // 27
    AtomKind(name: "CG2", radius: 2.00000, depth: 0.15000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77), // 28
    AtomKind(name: "CG3", radius: 2.00000, depth: 0.15000, hb_depth: 0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77), // 29
    AtomKind(name: "W", radius: 0.00000, depth: 0.00000, hb_depth: 0.0, hb_radius: 0.0, solvation: 0.00000, volume: 0.00000, covalent_radius: 0.00)  // 30
]

let metal_solvation_parameter: fl = -0.00110
let metal_covalent_radius:     fl = 1.75 // for metals not on the list // FIXME this info should be moved to non_ad_metals

let atom_kinds_size = AtomKindData.count

struct AtomEquivalence {
    let name: String
    let to: String
}

let AtomEquivalenceData: [AtomEquivalence] = [
    AtomEquivalence(name: "Se", to: "S")
]

let atom_equivalence_size = AtomEquivalenceData.count

struct AcceptorKind {
    let ad_type: sz
    let radius: fl
    let depth: fl
}

let AcceptorKindData: [AcceptorKind] = [
    //           ad_type,             radius,      depth
    AcceptorKind(ad_type: AD_TYPE_NA, radius: 1.9, depth: 5.0),
    AcceptorKind(ad_type: AD_TYPE_OA, radius: 1.9, depth: 5.0),
    AcceptorKind(ad_type: AD_TYPE_SA, radius: 2.5, depth: 1.0)
]

let acceptor_kinds_size = AcceptorKindData.count

@inline(__always)
func ad_is_hydrogen(ad: sz) -> Bool {
    return ad == AD_TYPE_H || ad == AD_TYPE_HD
}

@inline(__always)
func ad_is_heteroatom(ad: sz) -> Bool {
    return ad != AD_TYPE_A && ad != AD_TYPE_C  && 
           ad != AD_TYPE_H && ad != AD_TYPE_HD && 
           ad < AD_TYPE_SIZE
}

func ad_type_to_el_type(t: sz) throws -> sz { 
    switch t {
    case AD_TYPE_C    : return EL_TYPE_C
    case AD_TYPE_A    : return EL_TYPE_C
    case AD_TYPE_N    : return EL_TYPE_N
    case AD_TYPE_O    : return EL_TYPE_O
    case AD_TYPE_P    : return EL_TYPE_P
    case AD_TYPE_S    : return EL_TYPE_S
    case AD_TYPE_H    : return EL_TYPE_H
    case AD_TYPE_F    : return EL_TYPE_F
    case AD_TYPE_I    : return EL_TYPE_I
    case AD_TYPE_NA   : return EL_TYPE_N
    case AD_TYPE_OA   : return EL_TYPE_O
    case AD_TYPE_SA   : return EL_TYPE_S
    case AD_TYPE_HD   : return EL_TYPE_H
    case AD_TYPE_Mg   : return EL_TYPE_Met
    case AD_TYPE_Mn   : return EL_TYPE_Met
    case AD_TYPE_Zn   : return EL_TYPE_Met
    case AD_TYPE_Ca   : return EL_TYPE_Met
    case AD_TYPE_Fe   : return EL_TYPE_Met
    case AD_TYPE_Cl   : return EL_TYPE_Cl
    case AD_TYPE_Br   : return EL_TYPE_Br
    case AD_TYPE_Si   : return EL_TYPE_Si
    case AD_TYPE_At   : return EL_TYPE_At
    case AD_TYPE_CG0  : return EL_TYPE_C
    case AD_TYPE_CG1  : return EL_TYPE_C
    case AD_TYPE_CG2  : return EL_TYPE_C
    case AD_TYPE_CG3  : return EL_TYPE_C
    case AD_TYPE_G0   : return EL_TYPE_Dummy
    case AD_TYPE_G1   : return EL_TYPE_Dummy
    case AD_TYPE_G2   : return EL_TYPE_Dummy
    case AD_TYPE_G3   : return EL_TYPE_Dummy
    case AD_TYPE_W    : return EL_TYPE_Dummy
    case AD_TYPE_SIZE : return EL_TYPE_SIZE
    default: throw VinaError("ad_type_to_el_type: unknown type \(t)")
    }
}

let xs_vdw_radii: [fl] = [ 
    1.9, // C_H
    1.9, // C_P
    1.8, // N_P
    1.8, // N_D
    1.8, // N_A
    1.8, // N_DA
    1.7, // O_P
    1.7, // O_D
    1.7, // O_A
    1.7, // O_DA
    2.0, // S_P
    2.1, // P_P
    1.5, // F_H
    1.8, // Cl_H
    2.0, // Br_H
    2.2, // I_H
    2.2, // Si
    2.3, // At
    1.2, // Met_D
    1.9, // C_H_CG0
    1.9, // C_P_CG0
    1.9, // C_H_CG1
    1.9, // C_P_CG1
    1.9, // C_H_CG2
    1.9, // C_P_CG2
    1.9, // C_H_CG3
    1.9, // C_P_CG3
    0.0, // G0
    0.0, // G1
    0.0, // G2
    0.0, // G3
    0.0  // W
]

public let xs_vinardo_vdw_radii: [fl] = [
    2.0, // C_H
    2.0, // C_P
    1.7, // N_P
    1.7, // N_D
    1.7, // N_A
    1.7, // N_DA
    1.6, // O_P
    1.6, // O_D
    1.6, // O_A
    1.6, // O_DA
    2.0, // S_P
    2.1, // P_P
    1.5, // F_H
    1.8, // Cl_H
    2.0, // Br_H
    2.2, // I_H
    2.2, // Si
    2.3, // At
    1.2, // Met_D
    2.0, // C_H_CG0
    2.0, // C_P_CG0
    2.0, // C_H_CG1
    2.0, // C_P_CG1
    2.0, // C_H_CG2
    2.0, // C_P_CG2
    2.0, // C_H_CG3
    2.0, // C_P_CG3
    0.0, // G0
    0.0, // G1
    0.0, // G2
    0.0, // G3
    0.0  // W
]

let non_ad_metal_names: [String] = [ // expand as necessary
    "Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
]

@inline(__always)
func is_non_ad_metal_name(name: String) -> Bool {
    return non_ad_metal_names.contains(name)
}

@inline(__always)
func xs_is_hydrophobic(_ xs: sz) -> Bool {
    return xs == XS_TYPE_C_H ||
           xs == XS_TYPE_F_H ||
           xs == XS_TYPE_Cl_H ||
           xs == XS_TYPE_Br_H || 
           xs == XS_TYPE_I_H
}

@inline(__always)
func xs_is_acceptor(_ xs: sz) -> Bool {
    return xs == XS_TYPE_N_A ||
           xs == XS_TYPE_N_DA ||
           xs == XS_TYPE_O_A ||
           xs == XS_TYPE_O_DA
}

@inline(__always)
func xs_is_donor(_ xs: sz) -> Bool {
    return xs == XS_TYPE_N_D ||
           xs == XS_TYPE_N_DA ||
           xs == XS_TYPE_O_D ||
           xs == XS_TYPE_O_DA ||
           xs == XS_TYPE_Met_D
}

@inline(__always)
func xs_donor_acceptor(t1: sz, t2: sz) -> Bool {
    return xs_is_donor(t1) && xs_is_acceptor(t2)
}

@inline(__always)
func xs_h_bond_possible(t1: sz, t2: sz) -> Bool {
    return xs_donor_acceptor(t1: t1, t2: t2) || xs_donor_acceptor(t1: t2, t2: t1)
}

@inline(__always)
public func ad_type_property(i: sz) -> AtomKind {
    assert(AD_TYPE_SIZE == atom_kinds_size)
    assert(i < atom_kinds_size)
    return AtomKindData[i]
}

@inline(__always)
func string_to_ad_type(name: String) -> sz { // returns AD_TYPE_SIZE if not found (no exceptions thrown, because metals unknown to AD4 are not exceptional)
    for i in 0..<atom_kinds_size {
        if AtomKindData[i].name == name {
            return sz(i)
        }
    }
    for i in 0..<atom_equivalence_size {
        if AtomEquivalenceData[i].name == name {
            return string_to_ad_type(name: AtomEquivalenceData[i].to)
        }
    }
    return AD_TYPE_SIZE
}

@inline(__always)
func max_covalent_radius() -> fl {
    var tmp: fl = 0.0
    for i in 0..<atom_kinds_size {
        if AtomKindData[i].covalent_radius > tmp {
            tmp = AtomKindData[i].covalent_radius
        }
    }
    return tmp
}
