//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation

class AtomType: NSCopying {
    enum T {
        case EL
        case AD
        case XS
        case SY
    }
    
    @inline(__always)
    var el, ad, xs, sy: sz
        
    init() {
        self.el = EL_TYPE_SIZE
        self.ad = AD_TYPE_SIZE
        self.xs = XS_TYPE_SIZE
        self.sy = SY_TYPE_SIZE
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let newAtomType = AtomType()
        newAtomType.xs = self.xs
        newAtomType.el = self.el
        newAtomType.ad = self.ad
        newAtomType.sy = self.sy
        return newAtomType
    }
    
    func get(_ atom_typing_used: T) -> sz {
        switch atom_typing_used {
        case .EL:
            return self.el
        case .AD:
            return self.ad
        case .XS:
            return self.xs
        case .SY:
            return self.sy
        }
    }

    func is_hydrogen() -> Bool {
        return ad_is_hydrogen(ad: ad)
    }

    func is_heteroatom() -> Bool {
        return ad_is_heteroatom(ad: ad) || xs == XS_TYPE_Met_D
    }

    func acceptable_type() -> Bool {
        return ad < AD_TYPE_SIZE || xs == XS_TYPE_Met_D
    }

    func assign_el() {
        do {
            el = try ad_type_to_el_type(t: ad)
            if ad == AD_TYPE_SIZE && xs == XS_TYPE_Met_D {
                el = EL_TYPE_Met
            }
        } catch {
            fatalError(error.localizedDescription)
        }
    }

    func same_element(_ a: AtomType) -> Bool {
        return el == a.el
    }

    func covalent_radius() -> fl {
        if ad < AD_TYPE_SIZE {
            return ad_type_property(i: ad).covalent_radius
        } else if xs == XS_TYPE_Met_D {
            return metal_covalent_radius
        }
        return 0
    }

    func optimal_covalent_bond_length(_ x: AtomType) -> fl {
        return covalent_radius() + x.covalent_radius()
    }
    
}

@inline(__always)
func num_atom_types(_ atom_typing_used: AtomType.T) -> sz {
    switch atom_typing_used {
    case .EL:
        return EL_TYPE_SIZE
    case .AD:
        return AD_TYPE_SIZE
    case .XS:
        return XS_TYPE_SIZE
    case .SY:
        return SY_TYPE_SIZE
    }
}

@inline(__always)
func get_type_pair_index(_ atom_typing_used: AtomType.T, _ a: AtomType, _ b: AtomType) -> sz { // throws error if any arg is unassigned in the given typing scheme
    let n = num_atom_types(atom_typing_used)

    let i = a.get(atom_typing_used)
    assert(i < n, "i < n")
    let j = b.get(atom_typing_used)
    assert(j < n, "j < n")

    if i <= j {
        return triangular_matrix_index(n, i, j)
    } else {
        return triangular_matrix_index(n, j, i)
    }
}

