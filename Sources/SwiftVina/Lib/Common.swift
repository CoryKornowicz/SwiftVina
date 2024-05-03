//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/11/23.
//

import Foundation
import simd

public typealias fl = Double
public typealias sz = UInt64

public typealias mat = simd_double3x3
public typealias vec = simd_double3

public let zero_vec = SIMD3<fl>.zero

typealias Pair<T: Any, U: Any> = (T, U)
typealias Triple<T: Any, U: Any, V: Any> = (T, U, V)
typealias Quad<T: Any, U: Any, V: Any, W: Any> = (T, U, V, W)

typealias pr = Pair<fl, fl>
typealias vecv = ContiguousArray<vec>
typealias vecp = Pair<vec, vec>
public typealias flv = ContiguousArray<fl>
typealias prv = ContiguousArray<pr>
typealias szv = ContiguousArray<sz>

public let pi: fl = 3.1415926535897931
public let two_pi = 2 * pi
public let three_pi = 3 * pi
public let fl_tolerance: fl = 0.001
public let epsilon_fl: fl = fl.ulpOfOne

public let max_fl: fl = fl.greatestFiniteMagnitude
public let max_vec = SIMD3<fl>(max_fl, max_fl, max_fl)

func normalize_angle(x: inout fl) { // subtract or add enough 2*pi's to make x be in [-pi, pi]
    if x > three_pi { // very large
        let n = (x - pi) / two_pi // how many 2*pi's do you want to subtract?
        x -= two_pi*ceil(n) // ceil can be very slow, but this should not be called often
        normalize_angle(x: &x)
    } else if x < -three_pi { // very small
        let n = (-x - pi) / two_pi // how many 2*pi's do you want to add?
        x += two_pi*ceil(n) // ceil can be very slow, but this should not be called often
        normalize_angle(x: &x)
    } else if x > pi { // in (   pi, 3*pi]
        x -= two_pi
    } else if x < -pi { // in [-3*pi,  -pi)
        x += two_pi
    }
    assert(x >= -pi && x <= pi)
}

@inline(__always)
func normalized_angle(x: fl) -> fl {
    var tmp_x = x
    normalize_angle(x: &tmp_x)
    return tmp_x
}

@_transparent
func not_max(_ x: fl) -> Bool {
    return x < 0.1 * fl.greatestFiniteMagnitude
}

func fl_to_sz(_ x: fl, max_sz: sz) -> sz {
    if(x <= 0) { return 0 }
    if(x >= fl(max_sz)) { return max_sz }
    let tmp: sz = sz(x)
    return min(tmp, max_sz) // sz -> fl cast loses precision. 'min'
}

@_transparent
func eq(_ a: fl, _ b: fl) -> Bool {
    return abs(a - b) < fl_tolerance
}

@_transparent
func eq(_ a: vec, _ b: vec) -> Bool {
    return eq(a.x, b.x) && eq(a.y, b.y) && eq(a.z, b.z)
}

@_transparent
func eq(_ a: Array<fl>, _ b: Array<fl>) -> Bool { 
    return a.count == b.count && zip(a, b).allSatisfy(eq)
}

@_transparent
func eq(_ a: Array<vec>, _ b: Array<vec>) -> Bool {
    return a.count == b.count && zip(a, b).allSatisfy(eq)
}

@_transparent
func cross_product(_ a: vec, _ b: vec) -> vec {
    return vec(a[1]*b[2] - a[2]*b[1],
               a[2]*b[0] - a[0]*b[2],
               a[0]*b[1] - a[1]*b[0])
}

//public func *(_ lhs: mat, _ v: vec) -> vec {
//    let x1 = lhs[0, 0] * v[0] + lhs[1, 0] * v[1] + lhs[2, 0] * v[2]
//    let x2 = lhs[0, 1] * v[0] + lhs[1, 1] * v[1] + lhs[2, 1] * v[2]
//    let x3 = lhs[0, 2] * v[0] + lhs[1, 2] * v[1] + lhs[2, 2] * v[2]
//    return vec(x1, x2, x3)
//}

// multiply pK by this to get free energy in kcal/mol:
// K = exp(E/RT)  -- lower K and lower E == better binder
// pK = -log10(K)   => K = 10^(-pK)
// E = RT ln(K) = RT ln (10^(-pK)) = - RT * ln(10) * pK
let pK_to_energy_factor: fl = -8.31 /* RT in J/K/mol */ * 0.001 /* kilo */ * 300 /* K */ / 4.184 /* J/cal */ * log(10.0)
//  -0.6 kcal/mol * log(10) = -1.38

@inline(__always)
func pK_to_energy(pK: fl) -> fl {
    return pK_to_energy_factor * pK
}

struct VinaError: Error {
    
    private let message: String
    
    private let line: Int
    private let file: String
    private let function: String
    private let column: Int
    
    init(_ message: String, file: String = #file, function: String = #function, line: Int = #line, column: Int = #column) {
        self.message = message
        self.line = line
        self.file = file
        self.function = function
        self.column = column
    }
}

// Extensions for vec & mat

extension vec {
    
    @_transparent
    func norm_sqr() -> fl {
        return sqr(x) + sqr(y) + sqr(z)
    }
    
    @_transparent
    func norm() -> fl {
        return sqrt(norm_sqr())
    }
    
    static func * (_ lhs: Self, _ rhs: Self) -> fl {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z
    }
    
}
