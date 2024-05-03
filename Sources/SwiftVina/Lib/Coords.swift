//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

func rmsd_upper_bound(_ a: vecv, _ b: vecv) -> fl {
    assert(a.count == b.count)
    let acc: fl = zip(a, b).map { vec_distance_sqr($0, $1) }.reduce(0, +)
    return (a.count > 0) ? sqrt(acc / fl(a.count)) : 0
}

func find_closest(_ a: vecv, _ b: output_container) -> Pair<sz, fl> {
    var tmp: Pair<sz, fl> = (sz(b.count), fl.greatestFiniteMagnitude)
    for i in 0..<b.count {
        let res = rmsd_upper_bound(a, b[i].coords)
        if i == 0 || res < tmp.1 {
            tmp = (sz(i), res)
        }
    }
    return tmp
}

func add_to_output_container(_ out: inout output_container, _ t: OutputType, _ min_rmsd: fl, _ max_size: sz) {
    let closest_rmsd = find_closest(t.coords, out)
    if closest_rmsd.0 < sz(out.count) && closest_rmsd.1 < min_rmsd { // have a very similar one
        if t.e < out[closest_rmsd.0].e { // the new one is better, apparently
            out[closest_rmsd.0] = t // MARK: maybe? c++ -> FIXME? slow
        }
    }
    else { // nothing similar
        if out.count < max_size {
            out.append(t) // the last one had the worst energy - replacing
        }
        else {
            if !out.isEmpty && t.e < out.last!.e { // FIXME? - just changed
                out[out.count-1] = t // FIXME? slow
            }
        }
    }
    out.sort()
}
