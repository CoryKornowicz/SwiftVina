//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

struct precalculate_element {
    
    var smooth: prv 
    private var fast: flv
    private var factor: fl

    init(n: sz, factor: fl) {
        self.smooth = prv(repeating: pr(0, 0), count: Int(n))
        self.fast = flv(repeating: 0, count: Int(n))
        self.factor = factor
    }
    
    @inlinable
    func eval_fast(r2: fl) -> fl {
        assert(r2 * factor < fl(fast.count))
        let i = sz(factor * r2)
        assert(i < fast.count)
        return fast[Int(i)]
    }

    func eval_deriv(r2: fl) -> pr {
        let r2_factored = factor * r2
        assert(r2_factored + 1 < fl(smooth.count))
        let i1 = sz(r2_factored)
        let i2 = i1 + 1
        assert(i1 < smooth.count)
        assert(i2 < smooth.count)
        let rem = r2_factored - fl(i1)
        assert(rem >= -epsilon_fl)
        assert(rem < 1 + epsilon_fl)
        let p1 = smooth[Int(i1)]
        let p2 = smooth[Int(i2)]
        let e = p1.0 + rem * (p2.0 - p1.0)
        let dor = p1.1 + rem * (p2.1 - p1.1)
        return pr(e, dor)
    }

    mutating func init_from_smooth_fst(rs: flv) {
        let n = smooth.count
        assert(rs.count == n)
        assert(fast.count == n)
        for i in 0..<n {
            if i == 0 || i == n - 1 {
                smooth[i].1 = 0
            } else {
                let delta = rs[i + 1] - rs[i - 1]
                let r = rs[i]
                smooth[i].1 = (smooth[i + 1].0 - smooth[i - 1].0) / (delta * r)
            }
            let f1 = smooth[i].0
            let f2 = (i + 1 >= n) ? 0 : smooth[i + 1].0
            fast[i] = (f2 + f1) / 2
        }
    }

    func min_smooth_fst() -> sz {
        var tmp = sz(0)
        for i_inv in 0..<smooth.count {
            let i = smooth.count - i_inv - 1
            if i_inv == 0 || smooth[i].0 < smooth[tmp].0 {
                tmp = sz(i)
            }
        }
        return tmp
    }

    mutating func widen_smooth_fst(rs: flv, left: fl, right: fl) {
        var tmp = flv(repeating: 0, count: smooth.count)
        let min_index = min_smooth_fst()
        assert(min_index < rs.count)
        assert(rs.count == smooth.count)
        let optimal_r = rs[min_index]
        for i in 0..<smooth.count {
            var r = rs[i]
            if r < optimal_r - left {
                r += left
            } else if r > optimal_r + right {
                r -= right
            } else {
                r = optimal_r
            }
            if r < 0 {
                r = 0
            }
            if r > rs.last! {
                r = rs.last!
            }
            tmp[i] = eval_deriv(r2: sqr(r)).0
        }
        for i in 0..<smooth.count {
            smooth[i].0 = tmp[i]
        }
    }

    mutating func widen(rs: flv, left: fl, right: fl) {
        widen_smooth_fst(rs: rs, left: left, right: right)
        init_from_smooth_fst(rs: rs)
    }

}

final class precalculate {
    
    private var m_cutoffSqr: fl
    private var m_maxCutoffSqr: fl
    private var m_n: sz
    private var m_factor: fl

    private var m_data: triangular_matrix<precalculate_element>
    
    private var verbosityLevel: VerbosityLevel = .critical

    init(sf: ScoringFunction, v: fl = fl.greatestFiniteMagnitude, factor: fl = 32.0, verbosityLabel: VerbosityLevel = .critical) {
        self.m_factor = factor
        self.m_cutoffSqr = sqr(sf.get_cutoff())
        self.m_maxCutoffSqr = sqr(sf.get_max_cutoff())
        self.m_n = sz(m_factor * m_maxCutoffSqr) + 3 // sz(factor * r^2) + 1 <= sz(factor * max_cutoff_sqr) + 2 <= n-1 < n  // see assert below
        self.verbosityLevel = verbosityLabel
        
        if self.verbosityLevel <= .debug {
            print("-- DEBUG -- sf.cutoff^2 in precalculate = \(self.m_maxCutoffSqr)")
        }
        
        let n: Int = Int(num_atom_types(sf.get_atom_typing()))
        self.m_data = triangular_matrix<precalculate_element>(n: sz(n), filler_val: precalculate_element(n: m_n, factor: m_factor))

        assert(m_factor > epsilon_fl)
        assert(sz(m_maxCutoffSqr * m_factor) + 1 < m_n) // max_cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
        assert(m_maxCutoffSqr * m_factor + 1 < fl(m_n))

        let rs = calculate_rs()
                
        for t1 in 0..<self.m_data.dim() {
            for t2 in t1..<self.m_data.dim() {
                self.m_data[t1, t2].smooth.withUnsafeMutableBufferPointer { bufferPointer in
                    for i in 0..<self.m_data[t1, t2].smooth.count {
                        let eval = min(v, sf.eval(t1, t2, rs[i]))
//                        self.m_data[t1, t2].smooth[i].0 = eval
                        bufferPointer[i].0 = eval
                    }
                }
                // init the rest
                self.m_data[t1, t2].init_from_smooth_fst(rs: rs)
            }
        }
        
    }

    @inlinable
    func eval_fast(type_pair_index: sz, r2: fl) -> fl {
        assert(r2 <= m_maxCutoffSqr)
        return m_data[type_pair_index].eval_fast(r2: r2)
    }

    func eval_deriv(type_pair_index: sz, r2: fl) -> pr {
        assert(r2 <= m_maxCutoffSqr)
        return m_data[type_pair_index].eval_deriv(r2: r2)
    }

    func index_permissive(t1: sz, t2: sz) -> sz {
        return m_data.index_permissive(i: t1, j: t2)
    }
    
    func cutoff_sqr() -> fl {
        return m_cutoffSqr
    }

    func max_cutoff_sqr() -> fl {
        return m_maxCutoffSqr
    }

    func widen(left: fl, right: fl) { 
        let rs = calculate_rs()
        for t1 in 0..<m_data.dim() {
            for t2 in t1..<m_data.dim() {
                m_data[t1, t2].widen(rs: rs, left: left, right: right)
            }
        }
    }
        
    private func calculate_rs() -> flv {
        var tmp = flv(repeating: 0, count: Int(m_n))
        for i in 0..<m_n {
            tmp[i] = sqrt(fl(i) / m_factor)
        }
        return tmp
    }

}

final class precalculate_byatom {

    private var m_cutoffSqr: fl
    private var m_maxCutoffSqr: fl
    private var m_n: sz
    private var m_factor: fl

    private var m_data: triangular_matrix<precalculate_element>
    private var verbosityLevel: VerbosityLevel = .critical

    init(sf: ScoringFunction, model: Model, v: fl = fl.greatestFiniteMagnitude, factor: fl = 32, verbosityLabel: VerbosityLevel = .critical) {
        // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
        self.m_factor = factor
        self.m_cutoffSqr = sqr(sf.get_cutoff())
        self.m_maxCutoffSqr = sqr(sf.get_max_cutoff())
        self.m_n = sz(m_factor * m_maxCutoffSqr) + 3 // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
        self.verbosityLevel = verbosityLabel
        
        if self.verbosityLevel <= .debug {
            print("-- DEBUG -- sf.cutoff^2 in precalculate = \(self.m_maxCutoffSqr)")
        }
        
        let n_atoms = model.num_atoms
        let atoms = model.get_atoms()
        self.m_data = triangular_matrix<precalculate_element>(n: n_atoms, filler_val: precalculate_element(n: m_n, factor: m_factor))

        assert(m_factor > epsilon_fl)
        assert(sz(m_maxCutoffSqr * m_factor) + 1 < m_n) // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
        assert(m_maxCutoffSqr * m_factor + 1 < fl(m_n))

        let rs = calculate_rs()

        for i in 0..<m_data.dim() {
            for j in i..<m_data.dim() {
                self.m_data[i, j].smooth.withUnsafeMutableBufferPointer { bufferPointer in
                    for k in 0..<self.m_data[i, j].smooth.count {
                        let eval = sf.eval(atoms[i], atoms[j], rs[k])
//                        self.m_data[i, j].smooth[k].0 = eval
                        bufferPointer[k].0 = eval
                    }
                }
                // init the rest
                m_data[i, j].init_from_smooth_fst(rs: rs)
            }
        }
        
    }

    func eval_fast(i: sz, j: sz, r2: fl) -> fl {
        assert(r2 <= m_maxCutoffSqr)
        return m_data[i, j].eval_fast(r2: r2)
    }

    func eval_deriv(i: sz, j: sz, r2: fl) -> pr {
        assert(r2 <= m_maxCutoffSqr)
        return m_data[i, j].eval_deriv(r2: r2)
    }

    func index_permissive(i: sz, j: sz) -> sz {
        return m_data.index_permissive(i: i, j: j)
    }
    
    func cutoff_sqr() -> fl {
        return m_cutoffSqr
    }

    func max_cutoff_sqr() -> fl {
        return m_maxCutoffSqr
    }

    func get_factor() -> sz {
        return sz(m_factor)
    }

    func widen(left: fl, right: fl) { 
        let rs = calculate_rs()
        for t1 in 0..<m_data.dim() {
            for t2 in t1..<m_data.dim() {
                m_data[t1, t2].widen(rs: rs, left: left, right: right)
            }
        }
    }

    @inline(__always)
    private func calculate_rs() -> flv {
        var tmp = flv(repeating: 0, count: Int(m_n))
        for i in 0..<m_n {
            tmp[i] = sqrt(fl(i) / m_factor)
        }
        return tmp
    }

}
