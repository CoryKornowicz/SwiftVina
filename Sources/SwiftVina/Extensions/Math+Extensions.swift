//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation
import simd

@_transparent
@_specialize(where T == Double)
@_specialize(where T == Float)
public func sqr<T: FloatingPoint>(_ x: T) -> T {
    return x*x
}

@_transparent
public func sqr(_ x: vec) -> fl {
    return sqr(x[0]) + sqr(x[1]) + sqr(x[2])
}

@_transparent
public func vec_distance_sqr(_ x: vec, _ y: vec) -> fl {
    return (sqr(x[0] - y[0]) +
            sqr(x[1] - y[1]) +
            sqr(x[2] - y[2]))
}

@_transparent
func curl(e: inout fl, deriv: inout fl, v: fl) {
    if e > 0 && not_max(v) { // FIXME authentic_v can be gotten rid of everywhere now
        let tmp: fl = (v < epsilon_fl) ? 0 : (v / (v + e))
        e *= tmp
        deriv *= sqr(tmp)
    }
}

@_transparent
func curl(e: inout fl, deriv: inout vec, v: fl) {
    if e > 0 && not_max(v) { // FIXME authentic_v can be gotten rid of everywhere now
        let tmp: fl = (v < epsilon_fl) ? 0 : (v / (v + e))
        e *= tmp
        deriv *= sqr(tmp)
    }
}

@_transparent
func curl(e: inout fl, v: fl) {
    if e > 0 && not_max(v) {
        e *= (v < epsilon_fl) ? 0 : (v / (v + e))
    }
}

