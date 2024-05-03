//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

@usableFromInline
typealias flmat = triangular_matrix<fl>

func minus_mat_vec_product(_ m: flmat, in_: Change, out: inout Change) {
    let n = m.dim()
    for i in 0..<n {
        var sum: fl = 0
        for j in 0..<n {
            sum += m[m.index_permissive(i: i, j: j)] * in_[j]
        }
        out[i] = -sum
    }
}

@inline(__always)
func scalar_product(a: Change, b: Change, n: sz) -> fl {
    var tmp: fl = 0
    for i in 0..<n {
        tmp += a[i] * b[i]
    }
    return tmp
}

@discardableResult
func bfgs_update(h: inout flmat, p: Change, y: Change, alpha: fl) -> Bool {
    let yp: fl = scalar_product(a: y, b: p, n: h.dim())
    if alpha * yp < epsilon_fl { return false } // fixme?? from original impl
    var minus_hy: Change = y
    minus_mat_vec_product(h, in_: y, out: &minus_hy)
    let yhy: fl = -scalar_product(a: y, b: minus_hy, n: h.dim())
    let r: fl = 1 / (alpha * yp) // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
    let n: sz = p.num_floats()
    for i in 0..<n {
        for j in i..<n { // includes i
            h[i, j] +=  alpha * r * (minus_hy[i] * p[j]
                                   + minus_hy[j] * p[i])
                      + alpha * alpha * (r * r * yhy + r) * p[i] * p[j] // s * s == alpha * alpha * p * p
        }
    }
    return true
}

func line_search(f: QuasiNewtonAux, n: sz, x: Conf, g: Change, f0: fl, p: Change, x_new: inout Conf, g_new: inout Change, f1: inout fl, evalcount: inout Int) -> fl { // returns alpha
    let c0: fl = 0.0001
    let max_trials: UInt = 10
    let multiplier: fl = 0.5
    var alpha: fl = 1
    
    let pg: fl = scalar_product(a: p, b: g, n: n)
    
    for _ in 0..<max_trials {
        x_new = x
        x_new.increment(p, alpha)
        f1 = f.operate(&x_new, &g_new)
        evalcount += 1
        if f1 - f0 < c0 * alpha * pg { break } // FIXME check - div by norm(p) ? no?
        alpha *= multiplier
    }
    return alpha
}

@inlinable
func set_diagonal(_ m: inout flmat, x: fl) {
    for i in 0..<m.dim() {
        m[i, i] = x
    }
}

@inlinable
func subtract_change(_ b: inout Change, _ a: inout Change, n: sz) { // b -= a
    for i in 0..<n {
        b[i] -= a[i]
    }
}


func bfgs(_ f: inout QuasiNewtonAux, x: inout Conf, g: inout Change, max_steps: UInt, average_required_improvement: fl, over: sz, evalcount: inout Int) -> fl {
    // x is I/O, final value is returned
    
    let n: sz = g.num_floats()
    var h: flmat = flmat(n: n, filler_val: 0)
    set_diagonal(&h, x: 1)
    
    var g_new: Change = g
    var x_new: Conf = x
    
    var f0: fl = f.operate(&x, &g)
    evalcount += 1
    
    let f_orig: fl = f0
    let g_orig: Change = g
    let x_orig: Conf = x
    
    var p: Change = g
    
    var f_values: flv = flv()
    f_values.append(f0)
    
    for step in 0..<max_steps {
        // print("BFGS Step: \(step)")
        minus_mat_vec_product(h, in_: g, out: &p)
        var f1: fl = 0
        let alpha: fl = line_search(f: f, n: n, x: x, g: g, f0: f0, p: p, x_new: &x_new, g_new: &g_new, f1: &f1, evalcount: &evalcount)
        var y: Change = g_new
        subtract_change(&y, &g, n: n)
        
        f_values.append(f1)
        f0 = f1
        x = x_new
        
        let scalarValue = sqrt(scalar_product(a: g, b: g, n: n))
        
        if !(scalarValue >= 1e-5) || scalarValue.isNaN { 
//            print("scalar product value in break: \(scalarValue)")
            break
        }
        
        g = g_new // ? original impl
        
        if step == 0 {
            let yy: fl = scalar_product(a: y, b: y, n: n)
            if abs(yy) > epsilon_fl {
                set_diagonal(&h, x: alpha * scalar_product(a: y, b: p, n: n) / yy)
            }
        }
        
        bfgs_update(h: &h, p: p, y: y, alpha: alpha)
    }
    
    if(!(f0 <= f_orig)) { // maybe succeeds for nans too
        f0 = f_orig
        x = x_orig
        g = g_orig
    }
    
    return f0
}
