//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation
import simd

@inline(__always)
public func random_fl(_ a: fl, _ b: fl) -> fl { // expects a < b, returns rand in [a, b]
    assert(a < b)
    return fl.random(in: a...b)
}

// Uses the Box-Muller transform to generate a random normal distribution
func random_normal(mean: fl, sigma: fl) -> fl {
    assert(sigma >= 0, "Sigma must be non-negative")

    let u1 = fl.random(in: fl.ulpOfOne...1) // Avoid 0 for log
    let u2 = fl.random(in: 0...1)
    
    let z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)
    return z0 * sigma + mean
}

@inline(__always)
func random_int(_ a: Int, _ b: Int) -> Int { // expects a <= b, returns rand in [a, b]
    assert(a <= b)
    return Int.random(in: a...b)
}

@inline(__always)
func random_sz(_ a: sz, _ b: sz) -> sz { // expects a <= b, returns rand in [a, b]
    assert(a <= b)
    let i = sz.random(in: a...b)
    return i
}

//func random_inside_sphere() -> vec {
//    let radius = pow(fl.random(in: 0...1), 1.0 / 3.0) // Cube root to ensure uniform distribution
//    let theta = fl.random(in: 0...2 * .pi) // Azimuth angle
//    let phi = acos(1 - 2 * fl.random(in: 0...1)) // Inclination angle
//
//    let x = radius * sin(phi) * cos(theta)
//    let y = radius * sin(phi) * sin(theta)
//    let z = radius * cos(phi)
//
//    let tmp = vec(x, y, z)
//    
//    if eq(simd_length(tmp), 1) {
//        return tmp
//    } else {
//        return random_inside_sphere()
//    }
//}

func random_inside_sphere() -> vec {
    while true { // on average will need to be run twice
        let tmp = vec(random_fl(-1, 1),
                      random_fl(-1, 1),
                      random_fl(-1, 1))
        
        if sqr(tmp) < 1 {
            return tmp
        }
    }
}

@inline(__always)
func random_in_box(_ corner1: vec, _ corner2: vec) -> vec { // expects corner1[i] < corner2[i]
    return vec(random_fl(corner1[0], corner2[0]),
               random_fl(corner1[1], corner2[1]),
               random_fl(corner1[2], corner2[2]))
}
