//
//  ScoringFunction.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

public enum SFChoice {
    case SF_VINA
    case SF_AD42
    case SF_VINARDO
}

final class ScoringFunction {
    
    private var potentials: ContiguousArray<Potential> = ContiguousArray<Potential>()
    private var confIndependents: ContiguousArray<ConfIndependent> = ContiguousArray<ConfIndependent>()
    private var weights: ContiguousArray<fl>
    private var cutoff: fl
    private var maxCutoff: fl
    private var atomTyping: AtomType.T

    private var numConfIndependents: Int {
        confIndependents.count
    }
    
    private var numPotentials: Int {
        potentials.count
    }
    
    init(scoring_function_choice: SFChoice, weights: ContiguousArray<fl>) {        
        switch scoring_function_choice {
        case .SF_VINA:
            potentials.append(vina_gaussian(offset_: 0, width_: 0.5, cutoff_: 8.0))
            potentials.append(vina_gaussian(offset_: 3, width_: 2.0, cutoff_: 8.0))
            potentials.append(vina_repulsion(offset_: 0.0, cutoff_: 8.0))
            potentials.append(vina_hydrophobic(good_: 0.5, bad_: 1.5, cutoff_: 8.0))
            potentials.append(vina_non_dir_h_bond(good_: -0.7, bad_: 0, cutoff_: 8.0))
            potentials.append(linearattraction(cutoff_: 20.0))
            confIndependents.append(num_tors_div())
            atomTyping = .XS
            cutoff = 8.0
            maxCutoff = 20.0
        case .SF_AD42:
            potentials.append(ad4_vdw(smoothing_: 0.5, cap_: 100000, cutoff_: 8.0))
            potentials.append(ad4_hb(smoothing_: 0.5, cap_: 100000, cutoff_: 8.0))
            potentials.append(ad4_electrostatic(cap_: 100, cutoff_: 20.48))
            potentials.append(ad4_solvation(desolvation_sigma_: 3.6, solvation_q_: 0.01097, charge_dependent_: true, cutoff_: 20.48))
            potentials.append(linearattraction(cutoff_: 20.0))
            confIndependents.append(ad4_tors_add())
            atomTyping = .AD
            cutoff = 20.48
            maxCutoff = 20.48
        case .SF_VINARDO:
            potentials.append(vinardo_gaussian(offset_: 0, width_: 0.8, cutoff_: 8.0))
            potentials.append(vinardo_repulsion(offset_: 0, cutoff_: 8.0))
            potentials.append(vinardo_hydrophobic(good_: 0, bad_: 2.5, cutoff_: 8.0))
            potentials.append(vinardo_non_dir_h_bond(good_: -0.6, bad_: 0, cutoff_: 8.0))
            potentials.append(linearattraction(cutoff_: 20.0))
            confIndependents.append(num_tors_div())
            atomTyping = .XS
            cutoff = 8.0
            maxCutoff = 20.0
        }

        self.weights = weights
    }

    @inlinable
    func eval(_ a: Atom, _ b: Atom, _ r: fl) -> fl { // intentionally not checking for cutoff
        var sum: fl = 0
//        for (weight, potential) in zip(weights[..<numPotentials], potentials) {
//            sum += weight * potential.eval(a: a, b: b, r: r)
//        }
        for i in 0..<numPotentials {
            sum += weights[i] * potentials[i].eval(a: a, b: b, r: r)
        }
        return sum
    }

    @inlinable
    func eval(_ t1: sz, _ t2: sz, _ r: fl) -> fl {
        var sum: fl = 0
//        for (weight, potential) in zip(weights[..<numPotentials], potentials) {
//            sum += weight * potential.eval(t1: t1, t2: t2, r: r)
//        }
        for i in 0..<numPotentials {
            sum += weights[i] * potentials[i].eval(t1: t1, t2: t2, r: r)
        }
        return sum
    }
    
    func conf_independent(_ m: Model, _ e: fl) -> fl {
        var e = e
        var it = Iterator<fl>(ContiguousArray<fl>(weights[numPotentials...]))
        let in_ = conf_independent_inputs(m: m)
        for i in 0..<numConfIndependents {
            // We don't accumulate energy. Why? I don't know...maybe we should??
            e = confIndependents[i].eval(in_, e, &it)
        }
        return e
    }
    
    func get_cutoff() -> fl {
        return cutoff
    }
    
    func get_max_cutoff() -> fl {
        return maxCutoff
    }
    
    func get_atom_typing() -> AtomType.T {
        return atomTyping
    }
    
    func get_atom_types() -> szv {
        var tmp: szv = szv()
        for i in 0..<num_atom_types(atomTyping) {
            tmp.append(i)
        }
        return tmp
    }
    
    func get_num_atom_types() -> sz {
        return num_atom_types(atomTyping)
    }

    func get_weights() -> flv {
        return weights
    }

}
