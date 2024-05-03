//
//  Functions.swift
//
//
//  Created by Cory Kornowicz on 11/14/23.
//

import Foundation
import simd
import Accelerate

@inlinable
func show_variable(name: String, value: Any) {
    print("\(name) = \(value)")
}

@inlinable
func pairwise_clash_penalty(r: fl, covalent_r: fl) -> fl {
    // r = 0          -> max_penalty
    // r = covalent_r -> 1
    // elsewhere      -> hyperbolic function
    assert(r >= 0)
    assert(covalent_r > epsilon_fl)
    let x: fl = r / covalent_r
    if x > 2 { return 0 }
    return 1 - x * x / 4
}


func eval_interacting_pairs(_ p: precalculate_byatom, _ v: fl, _ pairs: InteractingPairs, _ coords: vecv, with_max_cutoff: Bool = false) -> fl { // clean up
    var e: fl = 0
    var cutoff_sqr: fl = p.cutoff_sqr()
    
    if with_max_cutoff {
        cutoff_sqr = p.max_cutoff_sqr()
    }
    
    for i in 0..<pairs.count {
        let ip: InteractingPair = pairs[i]
        let r2: fl = vec_distance_sqr(coords[ip.a], coords[ip.b])
        if r2 < cutoff_sqr {
            var tmp: fl = p.eval_fast(i: ip.a, j: ip.b, r2: r2)
            curl(e: &tmp, v: v)
            e += tmp
        }
    }
    return e
}

func eval_interacting_pairs_deriv(_ p: precalculate_byatom, _ v: fl, _ pairs: InteractingPairs, _ coords: vecv, _ forces: inout vecv, with_max_cutoff: Bool = false) -> fl { // clean up
    var e: fl = 0
    var cutoff_sqr: fl = p.cutoff_sqr()
    
    if with_max_cutoff {
        cutoff_sqr = p.max_cutoff_sqr()
    }
    
    for i in 0..<pairs.count {
        let ip: InteractingPair = pairs[i]
        let r: vec = coords[ip.b] - coords[ip.a] // a -> b
        let r2: fl = sqr(r)
        if r2 < cutoff_sqr {
            var tmp: pr = p.eval_deriv(i: ip.a, j: ip.b, r2: r2)
            var force: vec = tmp.1 * r
            curl(e: &tmp.0, deriv: &force, v: v)
            e += tmp.0
            // FIXME inefficient, if using hard curl
            forces[ip.a] -= force // we could omit forces on inflex here
            forces[ip.b] += force
        }
    }
    return e
}
