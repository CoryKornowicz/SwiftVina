//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation
import simd

public typealias qt = simd_quatd
let qt_identity: qt = qt(ix: 0, iy: 0, iz: 0, r: 1)

@inline(__always)
func abs(_ q: qt) -> fl {
    let maxim = sup(q)

    if maxim == 0 {
        return maxim
    } else {
        let mixam = 1 / maxim

        let a = q.vector.x * mixam
        let b = q.vector.y * mixam
        let c = q.vector.z * mixam
        let d = q.vector.w * mixam

        return maxim * sqrt(a*a + b*b + c*c + d*d)
    }
}

@inline(__always)
func sup(_ q: qt) -> fl {
    return max(max(abs(q.vector.x), abs(q.vector.y)),
               max(abs(q.vector.z), abs(q.vector.w)))
}


@inline(__always)
func quaternion_is_normalized(_ q: qt) -> Bool {
    // Check if the length is close enough to 1 (considering floating-point precision)
    return eq(quaternion_norm_sqr(q), 1.0) && eq(abs(q), 1.0)
}

@inline(__always)
func eq(_ a: qt, _ b: qt) -> Bool {
    return (eq(a.vector.x, b.vector.x) &&
            eq(a.vector.y, b.vector.y) &&
            eq(a.vector.z, b.vector.z) &&
            eq(a.vector.w, b.vector.w))
}

@_transparent
func quaternion_norm_sqr(_ q: qt) -> fl {
    return sqr(q.vector.x) + sqr(q.vector.y) + sqr(q.vector.z) + sqr(q.vector.w)
}

@inline(__always)
func quaternion_normalize(_ q: inout qt) { 
    let s: fl = quaternion_norm_sqr(q)
    assert(eq(s, sqr(abs(q))))
    let a: fl = sqrt(s)
    assert(a > epsilon_fl)
    q *= 1 / a 
    assert(quaternion_is_normalized(q))
} 

@inline(__always)
func quaternion_normalize_approx(_ q: inout qt, tolerance: fl = 1e-6) { 
    let s: fl = quaternion_norm_sqr(q)
    assert(eq(s, sqr(abs(q))))
    if abs(s - 1) < tolerance {
        return 
    } else { 
        let a: fl = sqrt(s)
        assert(a > epsilon_fl)
        q *= 1 / a 
        assert(quaternion_is_normalized(q))
    }
}

func angle_to_quaternion(axis: vec, angle: fl) -> qt { // axis is assumed to be a unit vector
    assert(eq(axis.norm(), 1.0)) // axis is assumed to be a unit vector
    let norm_angle = normalized_angle(x: angle) // this is probably only necessary if angles can be very big
    let c: fl = cos(norm_angle / 2)
    let s: fl = sin(norm_angle / 2)
    return qt(ix: s*axis[0], iy: s*axis[1], iz: s*axis[2], r: c)
}

func angle_to_quaternion(rotation: vec) -> qt {
    let angle = rotation.norm()
    if angle > epsilon_fl {
        let axis: vec = (1/angle) * rotation
        return angle_to_quaternion(axis: axis, angle: angle)
    }
    return qt_identity
}

func quaternion_to_angle(q: qt) -> vec { 
    assert(quaternion_is_normalized(q))
    let c = q.real
    if c >= -1 && c <= 1 { // c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
        var angle = 2 * acos(c) // acos is in [0, pi]
        if angle > pi {
            angle -= two_pi // now angle is in [-pi, pi]
        }
        let axis = vec(x: q.imag.x, y: q.imag.y, z: q.imag.z)
        let s = sin(angle/2)
        if abs(s) < epsilon_fl {
            return zero_vec
        }
        return axis * (angle / s)
    } else { // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
        return zero_vec
    }
}

func quaternion_to_r3(q: qt) -> mat {
    assert(quaternion_is_normalized(q));

    var result = mat()
    
    let a: fl = q.vector.x
    let b: fl = q.vector.y
    let c: fl = q.vector.z
    let d: fl = q.vector.w

    let aa: fl = a*a
    let ab: fl = a*b
    let ac: fl = a*c
    let ad: fl = a*d
    let bb: fl = b*b
    let bc: fl = b*c
    let bd: fl = b*d
    let cc: fl = c*c
    let cd: fl = c*d
    let dd: fl = d*d

    assert(eq(aa+bb+cc+dd, 1))

    // from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
    result[0, 0] = (aa+bb-cc-dd)
    result[0, 1] = 2*(-ad+bc)
    result[0, 2] = 2*(ac+bd)

    result[1, 0] = 2*(ad+bc)
    result[1, 1] = (aa-bb+cc-dd)
    result[1, 2] = 2*(-ab+cd)

    result[2, 0] = 2*(-ac+bd)
    result[2, 1] = 2*(ab+cd)
    result[2, 2] = (aa-bb-cc+dd)
    
    return result
}

func random_orientation() -> qt {
    var q = qt(vector:[random_normal(mean: 0, sigma: 1),
                       random_normal(mean: 0, sigma: 1),
                       random_normal(mean: 0, sigma: 1),
                       random_normal(mean: 0, sigma: 1)])
    let nrm: fl = abs(q)
    
    if nrm > epsilon_fl {
        q /= nrm
        assert(quaternion_is_normalized(q))
        return q
    } else {
        return random_orientation() // this call should almost never happen
    }
}

@inline(__always)
func quaternion_increment(q: inout qt, rotation: vec) {
    assert(quaternion_is_normalized(q))
    q = angle_to_quaternion(rotation: rotation) * q
    quaternion_normalize_approx(&q)
}

@inline(__always)
func quaternion_difference(a: qt, b: qt) -> vec { // rotation that needs to be applied to convert a to b
    assert(quaternion_is_normalized(a))
    assert(quaternion_is_normalized(b))
    var tmp = b
    tmp /= a // b = tmp * a => b * inv(a) = tmp
    return quaternion_to_angle(q: tmp) // already assert normalization
}
