//
//  Noncache.swift
//
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

class Noncache: igrid {
    
    var slope: fl
    private var sgrid: SZVGrid
    private var gd: grid_dims
    private let p: precalculate
    
    init(m: Model, gd_: grid_dims, p_: precalculate, slope_: fl) {
        self.sgrid = SZVGrid(m, szv_grid_dims(gd_), p_.cutoff_sqr())
        self.gd = gd_
        self.p = p_
        self.slope = slope_
    }
    
    func within(m: Model, margin: fl = 0.0001) -> Bool {
        for i in 0..<m.num_movable_atoms { 
            if m.atoms[i].is_hydrogen() { continue }
            let a_coords = m.coords[i]
            for j in 0..<gd.count {
                if gd[j].n_voxels > 0 {
                    if a_coords[j] < gd[j].begin - margin || a_coords[j] > gd[j].end + margin {
                        return false
                    }
                }
            }
        }
        return true
    }
    
    func eval(m: Model, v: fl) -> fl {
        var e: fl = 0
        let cutoff_sqr = p.cutoff_sqr()
        let n = num_atom_types(.XS)
        for i in 0..<m.num_movable_atoms { 
            var this_e: fl = 0
            var out_of_bounds_penalty: fl = 0
            let a = m.atoms[i]
            var t1 = a.get(.XS)
            if t1 >= n { continue }
            switch t1 { 
            case XS_TYPE_G0, XS_TYPE_G1, XS_TYPE_G2, XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0, XS_TYPE_C_H_CG1, XS_TYPE_C_H_CG2, XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H
            case XS_TYPE_C_P_CG0, XS_TYPE_C_P_CG1, XS_TYPE_C_P_CG2, XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P
            default:
                break
            }   
            let a_coords = m.coords[i]
            var adjusted_a_coords: vec = a_coords
            for j in 0..<gd.count {
                if gd[j].n_voxels > 0 {
                    if a_coords[j] < gd[j].begin {
                        adjusted_a_coords[j] = gd[j].begin
                        out_of_bounds_penalty += abs(a_coords[j] - gd[j].begin)
                    } else if a_coords[j] > gd[j].end {
                        adjusted_a_coords[j] = gd[j].end
                        out_of_bounds_penalty += abs(a_coords[j] - gd[j].end)
                    }
                }
            }
            out_of_bounds_penalty *= slope
            let possibilities = sgrid.possibilities(adjusted_a_coords)
            for possibilities_j in 0..<possibilities.count {
                let j = possibilities[possibilities_j]
                let b = m.grid_atoms[j]
                let t2 = b.xs
                if t2 >= n { continue }
                let r_ba: vec = adjusted_a_coords - b.coords
                let r2 = sqr(r_ba)
                if r2 < cutoff_sqr {
                    let type_pair_index = get_type_pair_index(.XS, a, b)
                    this_e += p.eval_fast(type_pair_index: type_pair_index, r2: r2)
                }
            }

            curl(e: &this_e, v: v)
            e += this_e + out_of_bounds_penalty
        }
        return e
    }

    func eval_deriv(m: Model, v: fl) -> fl {
        var e: fl = 0
        let cutoff_sqr = p.cutoff_sqr()
        let n = num_atom_types(.XS)
        for i in 0..<m.num_movable_atoms { 
            var this_e: fl = 0
            var deriv: vec = zero_vec
            var out_of_bounds_deriv: vec = zero_vec
            var out_of_bounds_penalty: fl = 0
            let a = m.atoms[i]
            var t1 = a.get(.XS)
            if t1 >= n { continue }
            switch t1 {
            case XS_TYPE_G0, XS_TYPE_G1, XS_TYPE_G2, XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0, XS_TYPE_C_H_CG1, XS_TYPE_C_H_CG2, XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H
            case XS_TYPE_C_P_CG0, XS_TYPE_C_P_CG1, XS_TYPE_C_P_CG2, XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P
            default:
                break
            }
            let a_coords = m.coords[i]
            var adjusted_a_coords: vec = a_coords
            for j in 0..<gd.count {
                if gd[j].n_voxels > 0 {
                    if a_coords[j] < gd[j].begin {
                        adjusted_a_coords[j] = gd[j].begin
                        out_of_bounds_deriv[j] = -1
                        out_of_bounds_penalty += abs(a_coords[j] - gd[j].begin)
                    } else if a_coords[j] > gd[j].end {
                        adjusted_a_coords[j] = gd[j].end
                        out_of_bounds_deriv[j] = 1
                        out_of_bounds_penalty += abs(a_coords[j] - gd[j].end)
                    }
                }
            }
            out_of_bounds_penalty *= slope
            out_of_bounds_deriv *= slope
            let possibilities = sgrid.possibilities(adjusted_a_coords)
            for possibilities_j in 0..<possibilities.count {
                let j = possibilities[possibilities_j]
                let b = m.grid_atoms[j]
                let t2 = b.xs
                if t2 >= n { continue }
                let r_ba: vec = adjusted_a_coords - b.coords
                let r2 = sqr(r_ba)
                if r2 < cutoff_sqr {
                    let type_pair_index = get_type_pair_index(.XS, a, b)
                    let e_dor = p.eval_deriv(type_pair_index: type_pair_index, r2: r2)
                    this_e += e_dor.0
                    deriv += e_dor.1 * r_ba
                }
            }
            curl(e: &this_e, deriv: &deriv, v: v)
            m.minus_forces[i] = deriv + out_of_bounds_deriv
            e += this_e + out_of_bounds_penalty
        }
        return e
    }
    
    func eval_intra(m: Model, v: fl) -> fl {
        var e: fl = 0
        let cutoff_sqr = p.cutoff_sqr()
        let n = num_atom_types(.XS)
        for i in 0..<m.num_movable_atoms {
            if m.is_atom_in_ligand(i) { continue } // we only want flex - rigid interactions
            var this_e: fl = 0
            var out_of_bounds_penalty: fl = 0
            let a = m.atoms[i]
            var t1 = a.xs
            if t1 >= n { continue }
            switch t1 {
            case XS_TYPE_G0, XS_TYPE_G1, XS_TYPE_G2, XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0, XS_TYPE_C_H_CG1, XS_TYPE_C_H_CG2, XS_TYPE_C_H_CG3:
                t1 = XS_TYPE_C_H
            case XS_TYPE_C_P_CG0, XS_TYPE_C_P_CG1, XS_TYPE_C_P_CG2, XS_TYPE_C_P_CG3:
                t1 = XS_TYPE_C_P
            default:
                break
            }
            let a_coords = m.coords[i]
            var adjusted_a_coords: vec = a_coords
            for j in 0..<gd.count {
                if gd[j].n_voxels > 0 {
                    if a_coords[j] < gd[j].begin {
                        adjusted_a_coords[j] = gd[j].begin
                        out_of_bounds_penalty += abs(a_coords[j] - gd[j].begin)
                    } else if a_coords[j] > gd[j].end {
                        adjusted_a_coords[j] = gd[j].end
                        out_of_bounds_penalty += abs(a_coords[j] - gd[j].end)
                    }
                }
            }
            out_of_bounds_penalty *= slope
            let possibilities = sgrid.possibilities(adjusted_a_coords)
            for possibilities_j in 0..<possibilities.count {
                let j = possibilities[possibilities_j]
                let b = m.grid_atoms[j]
                let t2 = b.xs
                if t2 >= n { continue }
                let r_ba: vec = adjusted_a_coords - b.coords
                let r2 = sqr(r_ba)
                if r2 < cutoff_sqr {
                    let type_pair_index = get_type_pair_index(.XS, a, b)
                    this_e += p.eval_fast(type_pair_index: type_pair_index, r2: r2)
                }
            }
            curl(e: &this_e, v: v)
            e += this_e + out_of_bounds_penalty
        }
        return e
    }

    
}
