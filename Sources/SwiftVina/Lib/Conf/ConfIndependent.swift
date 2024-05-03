//
//  File.swift
//
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

@_transparent
fileprivate func conf_smooth_div(_ x: fl, _ y: fl) -> fl {
    if (abs(x) < epsilon_fl) { return 0 }
    if (abs(y) < epsilon_fl) { return ((x*y > 0) ? max_fl : -max_fl) } // FIXME I hope -max_fl does not become NaN
    return x / y
}

final class conf_independent_inputs {
    var torsdof: fl // from TORSDOF keyword in pdbqt file
    var num_tors: fl
    var num_rotors: fl
    var num_heavy_atoms: fl
    var num_hydrophobic_atoms: fl
    var ligand_max_num_h_bonds: fl
    var num_ligands: fl
    var ligand_lengths_sum: fl
    
    init(m: Model) {
        self.torsdof = 0
        self.num_tors = 0
        self.num_rotors = 0
        self.num_heavy_atoms = 0
        self.num_hydrophobic_atoms = 0
        self.ligand_max_num_h_bonds = 0
        self.num_ligands = fl(m.num_ligands)
        self.ligand_lengths_sum = 0

        for i in 0..<m.num_ligands {
            let lig = m.get_ligand(i)
            self.ligand_lengths_sum += fl(m.ligand_length(i))
            self.torsdof += fl(lig.degrees_of_freedom)

            for j in lig.begin..<lig.end {
                let a = m.get_atom(j)
                if a.el != EL_TYPE_H {
                    let ar = atom_rotors(m, AtomIndex(i: j, in_grid: false))
                    self.num_tors += 0.5 * fl(ar)
                    if ar > 2 { self.num_rotors += 0.5 }
                    else { self.num_rotors += 0.5 * fl(ar) }

                    if xs_is_hydrophobic(a.xs) {
                        self.num_hydrophobic_atoms += 1
                    }

                    if xs_is_acceptor(a.xs) || xs_is_donor(a.xs) {
                        self.ligand_max_num_h_bonds += 1
                    }

                    self.num_heavy_atoms += 1
                }
            }
        }
    }
    
    func flv() -> flv {
        var tmp: flv = []
        tmp.append(num_tors)
        tmp.append(num_rotors)
        tmp.append(num_heavy_atoms)
        tmp.append(num_hydrophobic_atoms)
        tmp.append(ligand_max_num_h_bonds)
        tmp.append(num_ligands)
        tmp.append(ligand_lengths_sum)
        return tmp
    }

    private func num_bonded_heavy_atoms(_ m: Model, _ i: AtomIndex) -> UInt {
        var acc: UInt = 0
        let bonds = m.get_atom(i).bonds
        for j in 0..<bonds.count {
            let b = bonds[j]
            let a = m.get_atom(b.connected_atom_index)
            if !a.is_hydrogen() {
                acc += 1
            }
        }
        return acc
    }

    private func atom_rotors(_ m: Model, _ i: AtomIndex) -> UInt { // the number of rotatable bonds to heavy ligand atoms
        var acc: UInt = 0
        let bonds = m.get_atom(i).bonds
        for j in 0..<bonds.count {
            let b = bonds[j]
            let a = m.get_atom(b.connected_atom_index)
            if b.rotatable && !a.is_hydrogen() && num_bonded_heavy_atoms(m, b.connected_atom_index) > 1 { // not counting CH_3, etc
                acc += 1
            }
        }
        return acc
    }

    func get_names() -> [String] { // FIXME: maybe this should be static
        var tmp: [String] = []
        tmp.append("num_tors")
        tmp.append("num_rotors")
        tmp.append("num_heavy_atoms")
        tmp.append("num_hydrophobic_atoms")
        tmp.append("ligand_max_num_h_bonds")
        tmp.append("num_ligands")
        tmp.append("ligand_lengths_sum")
        return tmp
    }
}


protocol ConfIndependent {
    
    @inlinable
    func size() -> sz
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl
}

// Vina
final class num_tors_sqr : ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 0.1 * i.next() // [-1 .. 1]
        return x + weight * sqr(input.num_tors) / 5
    }
}


final class num_tors_sqrt : ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 0.1 * i.next() // [-1 .. 1]
        return x + weight * sqrt(input.num_tors) / sqrt(5.0)
    }
}

final class num_tors_div: ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 0.1 * (i.next() + 1) // weight is in [0..0.2]
        return conf_smooth_div(x, 1 + weight * input.num_tors / 5.0)
    }
}

final class ligand_length: ConfIndependent {
    
    @inlinable
    func size() -> sz { return 1 }
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = i.next()
        return x + weight * input.ligand_lengths_sum
    }
}

final class num_ligands: ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
        
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 1 * i.next() // weight is in [-1.. 1]
        return x + weight * input.num_ligands
    }
}

final class num_heavy_atoms_div: ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
        
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 0.05 * i.next()
        return conf_smooth_div(x, 1 + weight * input.num_heavy_atoms)
    }
}

final class num_heavy_atoms: ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 0.05 * i.next()
        return x + weight * input.num_heavy_atoms
    }
}

final class num_hydrophobic_atoms: ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
        
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = 0.05 * i.next()
        return x + weight * input.num_hydrophobic_atoms
    }
}

// AD42
final class ad4_tors_add: ConfIndependent {
    
    @inlinable func size() -> sz { return 1 }
    
    @inlinable
    func eval(_ input: conf_independent_inputs, _ x: fl, _ i: inout Iterator<fl>) -> fl {
        let weight: fl = i.next()
        return x + weight * input.torsdof
    }
}
