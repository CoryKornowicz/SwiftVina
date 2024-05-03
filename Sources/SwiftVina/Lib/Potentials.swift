//
//  Potentials.swift
//
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd
// Vina common functions

@inline(__always)
fileprivate func xs_radius(xs: sz) -> fl {
    assert(xs < XS_TYPE_SIZE)
    return xs_vdw_radii[xs]
}

@inline(__always)
fileprivate func xs_vinardo_radius(xs: sz) -> fl {
    assert(xs < XS_TYPE_SIZE)
    return xs_vinardo_vdw_radii[xs]
}


@inlinable
func slope_step(_ x_bad: fl, _ x_good: fl, _ x: fl) -> fl {
    if (x_bad < x_good) {
        if (x <= x_bad) { return 0 }
        if (x >= x_good) { return 1 }
    }
    else {
        if (x >= x_bad) { return 0 }
        if (x <= x_good) { return 1 }
    }
    return (x - x_bad) / (x_good - x_bad)
}

@inlinable
func is_glue_type(_ xs_t: sz) -> Bool {
    return (xs_t==XS_TYPE_G0) || (xs_t==XS_TYPE_G1) || (xs_t==XS_TYPE_G2) || (xs_t==XS_TYPE_G3)
}

@inline(__always)
fileprivate func optimal_distance(_ xs_t1: sz, _ xs_t2: sz) -> fl {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) { return 0.0 } // G0, G1, G2 or G3
    return xs_radius(xs: xs_t1) + xs_radius(xs: xs_t2)
}

@inlinable
func smooth_div(_ x: fl, _ y: fl) -> fl {
     if (abs(x) < epsilon_fl) { return 0 }
     if (abs(y) < epsilon_fl) { return ((x*y > 0) ? max_fl : -max_fl) } // FIXME I hope -max_fl does not become NaN
    // Lets get rid of those checks in swift
    return x / y
}

// Vinardo common functions
@inline(__always)
fileprivate func optimal_distance_vinardo(_ xs_t1: sz, _ xs_t2: sz) -> fl {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) { return 0.0 } // G0, G1, G2 or G3
    return xs_vinardo_radius(xs: xs_t1) + xs_vinardo_radius(xs: xs_t2)
}

// AD42 common functions
fileprivate func smoothen(_ r: fl, _ rij: fl, _ smoothing: fl) -> fl {
    let smoothing = smoothing * 0.5
    if (r > rij + smoothing) {
        return r - smoothing
    } else if (r < rij - smoothing) {
        return r + smoothing
    } else {
        return rij
    }
}

@inline(__always)
fileprivate func ad4_hb_eps(_ a: sz) -> fl {
    if (a < AD_TYPE_SIZE) { return ad_type_property(i: a).hb_depth }
    assert(false)
    return 0
}

@inline(__always)
fileprivate func ad4_hb_radius(_ a: sz) -> fl {
    if (a < AD_TYPE_SIZE) { return ad_type_property(i: a).hb_radius }
    assert(false)
    return 0
}

@inline(__always)
fileprivate func ad4_vdw_eps(_ a: sz) -> fl {
    if (a < AD_TYPE_SIZE) { return ad_type_property(i: a).depth }
    assert(false)
    return 0
}

@inline(__always)
fileprivate func ad4_vdw_radius(_ a: sz) -> fl {
    if (a < AD_TYPE_SIZE) { return ad_type_property(i: a).radius }
    assert(false)
    return 0
}

// Macrocycle - Vina and AD42
fileprivate func is_glued(_ xs_t1: sz, _ xs_t2: sz) -> Bool {
    return ((xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_H_CG0) ||
            (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_P_CG0) ||
            (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_H_CG0) ||
            (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_P_CG0) ||
            
            (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_H_CG1) ||
            (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_P_CG1) ||
            (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_H_CG1) ||
            (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_P_CG1) ||
            
            (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_H_CG2) ||
            (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_P_CG2) ||
            (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_H_CG2) ||
            (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_P_CG2) ||
            
            (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_H_CG3) ||
            (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_P_CG3) ||
            (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_H_CG3) ||
            (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_P_CG3))
}

protocol Potential {
    func eval(a: Atom, b: Atom, r: fl) -> fl
    func eval(t1: sz, t2: sz, r: fl) -> fl
}

// Vina
final class vina_gaussian: Potential {
    
    private var offset: fl
    private var width: fl
    private var cutoff: fl
    
    init(offset_: fl, width_: fl, cutoff_: fl) {
        self.offset = offset_
        self.width = width_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        return gauss(r - (optimal_distance(a.xs, b.xs) + offset)) // hard-coded to XS atom type
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        return gauss(r - (optimal_distance(t1, t2) + offset)) // hard-coded to XS atom type
    }
        
    @inline(__always)
    private func gauss(_ x: fl) -> fl {
        return exp(-sqr(x / width))
    }
    
}

final class vina_repulsion: Potential {
    private var offset: fl
    private var cutoff: fl
    
    init(offset_: fl, cutoff_: fl) {
        self.offset = offset_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        let d = r - (optimal_distance(a.xs, b.xs) + offset) // hard-coded to XS atom type
        if d > 0.0 {
            return 0.0
        }
        return d * d
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        let d = r - (optimal_distance(t1, t2) + offset) // hard-coded to XS atom type
        if d > 0.0 {
            return 0.0
        }
        return d * d
    }
 }

final class vina_hydrophobic: Potential {
    
    private var good: fl
    private var bad: fl
    private var cutoff: fl
    
    init(good_: fl, bad_: fl, cutoff_: fl) {
        self.good = good_
        self.bad = bad_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        if xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs) {
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs))
        } else {
            return 0.0
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2) {
            return slope_step(bad, good, r - optimal_distance(t1, t2))
        } else {
            return 0.0
        }
    }
        
}

final class vina_non_dir_h_bond: Potential {
    private var good: fl
    private var bad: fl
    private var cutoff: fl
    
    init(good_: fl, bad_: fl, cutoff_: fl) {
        self.good = good_
        self.bad = bad_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        if xs_h_bond_possible(t1: a.xs, t2: b.xs) {
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs))
        } else {
            return 0.0
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if xs_h_bond_possible(t1: t1, t2: t2) {
            return slope_step(bad, good, r - optimal_distance(t1, t2))
        } else {
            return 0.0
        }
    }
}

// Vinardo
final class vinardo_gaussian: Potential {
    
    private var offset: fl
    private var width: fl
    private var cutoff: fl
    
    init(offset_: fl, width_: fl, cutoff_: fl) {
        self.offset = offset_
        self.width = width_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        return gauss(r - (optimal_distance_vinardo(a.xs, b.xs) + offset)) // hard-coded to XS atom type
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        return gauss(r - (optimal_distance_vinardo(t1, t2) + offset)) // hard-coded to XS atom type
    }
    
    @inline(__always)
    private func gauss(_ x: fl) -> fl {
        return exp(-sqr(x / width))
    }
}

final class vinardo_repulsion: Potential {
    private var offset: fl
    private var cutoff: fl
    
    init(offset_: fl, cutoff_: fl) {
        self.offset = offset_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        let d = r - (optimal_distance_vinardo(a.xs, b.xs) + offset) // hard-coded to XS atom type
        if d > 0.0 {
            return 0.0
        }
        return d * d
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        let d = r - (optimal_distance_vinardo(t1, t2) + offset) // hard-coded to XS atom type
        if d > 0.0 {
            return 0.0
        }
        return d * d
    }
}

final class vinardo_hydrophobic: Potential {
    private var good: fl
    private var bad: fl
    private var cutoff: fl
    
    init(good_: fl, bad_: fl, cutoff_: fl) {
        self.good = good_
        self.bad = bad_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        if xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs) {
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs))
        } else {
            return 0.0
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2) {
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2))
        } else {
            return 0.0
        }
    }
}

final class vinardo_non_dir_h_bond: Potential {
    var good: fl
    var bad: fl
    var cutoff: fl
    
    init(good_: fl, bad_: fl, cutoff_: fl) {
        self.good = good_
        self.bad = bad_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if a.xs >= XS_TYPE_SIZE || b.xs >= XS_TYPE_SIZE {
            return 0.0
        }
        if xs_h_bond_possible(t1: a.xs, t2: b.xs) {
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs))
        } else {
            return 0.0
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if xs_h_bond_possible(t1: t1, t2: t2) {
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2))
        } else {
            return 0.0
        }
    }
}

// AD42
final class ad4_electrostatic: Potential {
    private var cap: fl
    private var cutoff: fl
    
    init(cap_: fl, cutoff_: fl) {
        self.cap = cap_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        let q1q2 = a.charge * b.charge * 332.0
        let B: fl = 78.4 + 8.5525
        let lB = -B * 0.003627
        let diel = -8.5525 + (B / (1 + 7.7839 * exp(lB * r)))
        if r < epsilon_fl {
            return q1q2 * cap / diel
        } else {
            return q1q2 * min(cap, 1.0 / (r * diel))
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        return 0.0
    }
}

final class ad4_solvation: Potential {
    private var desolvation_sigma: fl
    private var solvation_q: fl
    private var charge_dependent: Bool
    private var cutoff: fl
    
    init(desolvation_sigma_: fl, solvation_q_: fl, charge_dependent_: Bool, cutoff_: fl) {
        self.desolvation_sigma = desolvation_sigma_
        self.solvation_q = solvation_q_
        self.charge_dependent = charge_dependent_
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        let q1 = a.charge
        let q2 = b.charge
        assert(not_max(q1))
        assert(not_max(q2))
        let solv1 = solvation_parameter(a: a)
        let solv2 = solvation_parameter(a: b)
        let volume1 = volume(a: a)
        let volume2 = volume(a: b)
        let my_solv = charge_dependent ? solvation_q : 0
        let tmp = ((solv1 + my_solv * abs(q1)) * volume2 +
                   (solv2 + my_solv * abs(q2)) * volume1) * exp(-0.5 * sqr(r / desolvation_sigma))
        assert(not_max(tmp))
        return tmp
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        return 0.0
    }
    
    @inline(__always)
    private func volume(a: Atom) -> fl {
        if a.ad < AD_TYPE_SIZE {
            return ad_type_property(i: a.ad).volume
        } else if a.xs < XS_TYPE_SIZE {
            return ((4.0 / 3.0) * pi) * pow(xs_radius(xs: a.xs), 3)
        }
        assert(false)
        return 0.0
    }
    
    @inline(__always)
    private func solvation_parameter(a: Atom) -> fl {
        if a.ad < AD_TYPE_SIZE {
            return ad_type_property(i: a.ad).solvation
        } else if a.xs == XS_TYPE_Met_D {
            return metal_solvation_parameter
        }
        assert(false)
        return 0.0
    }
}

final class ad4_vdw: Potential {
    private var smoothing: fl
    private var cap: fl
    private var cutoff: fl
    
    init(smoothing_: fl, cap_: fl, cutoff_: fl) {
        self.smoothing = smoothing_
        self.cap = cap_
        self.cutoff = cutoff_
    }
    
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        let t1 = a.ad
        let t2 = b.ad
        let hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2)
        let vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2)
        let vdw_depth = sqrt(ad4_vdw_eps(t1) * ad4_vdw_eps(t2))
        if hb_depth < 0 {
            return 0.0 // interaction is hb, not vdw.
        }
        let r = smoothen(r, vdw_rij, smoothing)
        let c_12 = pow(vdw_rij, 12) * vdw_depth
        let c_6  = pow(vdw_rij, 6)  * vdw_depth * 2.0
        let r6   = pow(r, 6)
        let r12  = pow(r, 12)
        if r12 > epsilon_fl && r6 > epsilon_fl {
            return min(cap, (c_12 / r12) - (c_6 / r6))
        } else {
            return cap
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        return 0.0
    }
}

final class ad4_hb: Potential {
    private var smoothing: fl
    private var cap: fl
    private var cutoff: fl
    
    init(smoothing_: fl, cap_: fl, cutoff_: fl) {
        self.smoothing = smoothing_
        self.cap = cap_
        self.cutoff = cutoff_
    }
    
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        let t1 = a.ad
        let t2 = b.ad
        let hb_rij = ad4_hb_radius(t1) + ad4_hb_radius(t2)
        let hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2)
        if hb_depth >= 0 {
            return 0.0 // interaction is vdw, not hb.
        }
        let r = smoothen(r, hb_rij, smoothing)
        let c_12 = pow(hb_rij, 12) * -hb_depth * 10 / 2.0
        let c_10 = pow(hb_rij, 10) * -hb_depth * 12 / 2.0
        let r10  = pow(r, 10)
        let r12  = pow(r, 12)
        if r12 > epsilon_fl && r10 > epsilon_fl {
            return min(cap, c_12 / r12 - c_10 / r10)
        } else {
            return cap
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        return 0.0
    }
}

// Macrocycle - Vina and AD42
final class linearattraction: Potential {
    private var cutoff: fl
    
    init(cutoff_: fl) {
        self.cutoff = cutoff_
    }
    
    @inlinable
    func eval(a: Atom, b: Atom, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if is_glued(a.xs, b.xs) {
            return r
        } else {
            return 0.0
        }
    }
    
    @inlinable
    func eval(t1: sz, t2: sz, r: fl) -> fl {
        if r >= cutoff {
            return 0.0
        }
        if is_glued(t1, t2) {
            return r
        } else {
            return 0.0
        }
    }
}
