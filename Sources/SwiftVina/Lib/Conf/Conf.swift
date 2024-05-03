//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation
import simd

struct Scale {
    let position: fl
    let orientation: fl
    let torsion: fl
}

class ConfSize {
    var ligands: szv = szv()
    var flex: szv = szv() 
    func num_degrees_of_freedom() -> sz {
        let lig_size = ligands.reduce(0, +)
        let flex_size = flex.reduce(0, +)
        return lig_size + flex_size + sz(6 * ligands.count)
    }
}

@inline(__always)
func torsions_set_to_null(_ torsions: inout flv) {
    torsions = flv(repeating: 0, count: torsions.count)
}

@inline(__always)
func torsions_increment(_ torsions: inout flv, _ c: flv, _ factor: fl) {
    for i in 0..<torsions.count {
        torsions[i] += normalized_angle(x: factor * c[i])
        normalize_angle(x: &torsions[i])
    }
}

@_transparent
func torsions_randomize(_ torsions: inout flv) {
    for i in 0..<torsions.count {
        torsions[i] = random_fl(-pi, pi)
    }
}

// This function checks if all torsions are too close and returns true if they are
// if even one torsion is not too close, it returns false
// which might be troublesome...
@inline(__always)
func torsions_too_close(_ torsions1: flv, _ torsions2: flv, _ cutoff: fl) -> Bool {
    assert(torsions1.count == torsions2.count)
    for i in 0..<torsions1.count {
        if abs(normalized_angle(x: torsions1[i] - torsions2[i])) > cutoff {
            return false
        }
    }
    return true
}

func torsions_generate(_ torsions: inout flv, _ spread: fl, _ rp: fl, _ rs: flv? = nil) {
    assert(rs == nil || rs!.count == torsions.count)
    for i in 0..<torsions.count {
        if rs != nil && random_fl(0, 1) < rp {
            torsions[i] = rs![i]
        } else {
            torsions[i] += random_fl(-spread, spread)
        }
    }
}

public struct RigidChange {
    public var position: vec
    public var orientation: vec
    
    init() {
        self.position = zero_vec
        self.orientation = zero_vec
    }
    
    func print() {
        Swift.print(self.position)
        Swift.print(self.orientation)
    }
}

struct RigidConf {
    var position: vec
    var orientation: qt
    
    init() {
        self.position = zero_vec
        self.orientation = qt_identity
    }
    
    mutating func set_to_null() {
        position = zero_vec
        orientation = qt_identity
    }
    
    mutating func increment(_ c: RigidChange, _ factor: fl) {
        position += factor * c.position
        let rotation: vec = factor * c.orientation
        quaternion_increment(q: &orientation, rotation: rotation)
    }
    
    mutating func randomize(_ corner1: vec, _ corner2: vec) {
        position = random_in_box(corner1, corner2)
        orientation = random_orientation()
    }
    
    func too_close(_ c: RigidConf, _ position_cutoff: fl, _ orientation_cutoff: fl) -> Bool {
        if vec_distance_sqr(position, c.position) > sqr(position_cutoff) {
            return false
        }
        if sqr(quaternion_difference(a: orientation, b: c.orientation)) > sqr(orientation_cutoff) {
            return false
        }
        return true
    }
    
    mutating func mutate_position(_ spread: fl) {
        position += spread * random_inside_sphere()
    }
    
    mutating func mutate_orientation(_ spread: fl) {
        let tmp: vec = spread * random_inside_sphere()
        quaternion_increment(q: &orientation, rotation: tmp)
    }
    
    mutating func generate(_ position_spread: fl, _ orientation_spread: fl, _ rp: fl, _ rs: RigidConf? = nil) {
        if rs != nil && random_fl(0, 1) < rp {
            position = rs!.position
        } else {
            mutate_position(position_spread)
        }
        if rs != nil && random_fl(0, 1) < rp {
            orientation = rs!.orientation
        } else {
            mutate_orientation(orientation_spread)
        }
    }
    
    func apply(_ in_: vecv, _ out: inout vecv, _ begin: sz, _ end: sz) {
        assert(in_.count == out.count)
        let m = quaternion_to_r3(q: orientation)
        for i in begin..<end {
            out[i] = m * in_[i] + position
        }
    }

    func print() { 
        Swift.print(position)
        Swift.print(orientation)
    }
}

public struct LigandChange {
    public var rigid: RigidChange
    public var torsions: flv = flv()
    
    init(rigid: RigidChange) {
        self.rigid = rigid
    }
    
    func print() {
        self.rigid.print()
        Swift.print(self.torsions)
    }
}

struct LigandConf {
    var rigid: RigidConf
    var torsions: flv = flv()

    init(rigid: RigidConf) {
        self.rigid = rigid
    }
    
    mutating func set_to_null() {
        rigid.set_to_null()
        torsions_set_to_null(&torsions)
    }
    
    mutating func increment(_ c: LigandChange, _ factor: fl) {
        rigid.increment(c.rigid, factor)
        torsions_increment(&torsions, c.torsions, factor)
    }
    
    mutating func randomize(_ corner1: vec, _ corner2: vec) {
        rigid.randomize(corner1, corner2)
        torsions_randomize(&torsions)
    }

    func print() { 
        rigid.print()
        Swift.print(torsions)
    }
}

public struct ResidueChange {
    public var torsions: flv = flv()
    func print() {
        Swift.print(self.torsions)
    }
}

public struct ResidueConf {
    public var torsions: flv = flv()
    
    mutating func set_to_null() {
        torsions_set_to_null(&torsions)
    }
    
    mutating func increment(_ c: ResidueChange, _ factor: fl) {
        torsions_increment(&torsions, c.torsions, factor)
    }
    
    mutating func randomize() {
        torsions_randomize(&torsions)
    }
    
    func print() {
        Swift.print(torsions)
    }
}


public class Change {

    public var ligands: [LigandChange]
    public var flex: [ResidueChange]

    init(_ s: ConfSize) {
        ligands = [LigandChange](repeating: LigandChange(rigid: RigidChange()), count: s.ligands.count)
        flex = [ResidueChange](repeating: ResidueChange(), count: s.flex.count)
        for i in 0..<ligands.count {
            ligands[i].torsions = flv(repeating: 0, count: Int(s.ligands[i]))
        }
        for i in 0..<flex.count {
            flex[i].torsions = flv(repeating: 0, count: Int(s.flex[i]))
        }
    }
    
    @inlinable
    subscript(index: sz) -> fl {
        get {
            var index = index
            for i in 0..<ligands.count {
                let lig = ligands[i]
                if index < 3 {
                    return lig.rigid.position[Int(index)]
                }
                index -= 3
                if index < 3 {
                    return lig.rigid.orientation[Int(index)]
                }
                index -= 3
                if index < lig.torsions.count {
                    return lig.torsions[Int(index)]
                }
                index -= sz(lig.torsions.count)
            }
            for i in 0..<flex.count {
                let res = flex[i]
                if index < res.torsions.count {
                    return res.torsions[Int(index)]
                }
                index -= sz(res.torsions.count)
            }
            assert(false)
            return 0
        }
        
        set {
            var index = index
            for i in 0..<ligands.count {
                if index < 3 {
                    ligands[i].rigid.position[Int(index)] = newValue
                    return
                }
                index -= 3
                if index < 3 {
                    ligands[i].rigid.orientation[Int(index)] = newValue
                    return
                }
                index -= 3
                if index < ligands[i].torsions.count {
                    ligands[i].torsions[Int(index)] = newValue
                    return
                }
                index -= sz(ligands[i].torsions.count)
            }
            for i in 0..<flex.count {
                if index < flex[i].torsions.count {
                    flex[i].torsions[Int(index)] = newValue
                    return
                }
                index -= sz(flex[i].torsions.count)
            }
        }
    }

    func num_floats() -> sz { 
        var tmp: sz = 0
        for i in 0..<ligands.count {
            tmp += 6 + sz(ligands[i].torsions.count)
        }
        for i in 0..<flex.count {
            tmp += sz(flex[i].torsions.count)
        }
        return tmp
    }

    func print() {
        for i in 0..<ligands.count {
            ligands[i].print()
        }
        for i in 0..<flex.count {
            flex[i].print()
        }
    }
}

struct Conf {
    
    var ligands: [LigandConf]
    var flex: [ResidueConf]

    init(_ s: ConfSize) {
        ligands = [LigandConf](repeating: LigandConf(rigid: RigidConf()), count: s.ligands.count)
        flex = [ResidueConf](repeating: ResidueConf(), count: s.flex.count)
        for i in 0..<ligands.count {
            ligands[i].torsions = flv(repeating: 0, count: Int(s.ligands[i]))
        }
        for i in 0..<flex.count {
            flex[i].torsions = flv(repeating: 0, count: Int(s.flex[i]))
        }
    }

    mutating func set_to_null() {
        for i in 0..<ligands.count {
            ligands[i].set_to_null()
        }
        for i in 0..<flex.count {
            flex[i].set_to_null()
        }
    }

    mutating func increment(_ c: Change, _ factor: fl) { // torsions get normalized, orientations do not
        for i in 0..<ligands.count {
            ligands[i].increment(c.ligands[i], factor)
        }
        for i in 0..<flex.count {
            flex[i].increment(c.flex[i], factor)
        }
    }

    func internal_too_close(_ c: Conf, _ torsions_cutoff: fl) -> Bool {
        assert(ligands.count == c.ligands.count)
        for i in 0..<ligands.count {
            if !torsions_too_close(ligands[i].torsions, c.ligands[i].torsions, torsions_cutoff) {
                return false
            }
        }
        return true
    }

    func external_too_close(_ c: Conf, _ cutoff: Scale) -> Bool {
        assert(ligands.count == c.ligands.count)
        for i in 0..<ligands.count {
            if !ligands[i].rigid.too_close(c.ligands[i].rigid, cutoff.position, cutoff.orientation) {
                return false
            }
        }
        assert(flex.count == c.flex.count)
        for i in 0..<flex.count {
            if !torsions_too_close(flex[i].torsions, c.flex[i].torsions, cutoff.torsion) {
                return false
            }
        }
        return true
    }

    func too_close(_ c: Conf, _ cutoff: Scale) -> Bool {
        return internal_too_close(c, cutoff.torsion) && 
               external_too_close(c, cutoff) // a more efficient implementation is possible, probably
    }

    mutating func generate_internal(_ torsion_spread: fl, _ rp: fl, _ rs: Conf? = nil) { // torsions are not normalized after this
        for i in 0..<ligands.count {
            ligands[i].rigid.position = zero_vec
            ligands[i].rigid.orientation = qt_identity
            let torsions_rs = rs != nil ? rs!.ligands[i].torsions : nil
            torsions_generate(&ligands[i].torsions, torsion_spread, rp, torsions_rs)
        }
    }

    mutating func generate_external(_ spread: Scale, _ rp: fl, _ rs: Conf? = nil) { // torsions are not normalized after this
        for i in 0..<ligands.count {
            let rigid_conf_rs = rs != nil ? rs!.ligands[i].rigid : nil
            ligands[i].rigid.generate(spread.position, spread.orientation, rp, rigid_conf_rs)
        }
        for i in 0..<flex.count {
            let torsions_rs = rs != nil ? rs!.flex[i].torsions : nil
            torsions_generate(&flex[i].torsions, spread.torsion, rp, torsions_rs)
        }
    }

    mutating func randomize(_ corner1: vec, _ corner2: vec) {
        for i in 0..<ligands.count {
            ligands[i].randomize(corner1, corner2)
        }
        for i in 0..<flex.count {
            flex[i].randomize()
        }
    }

    func print() {
        for i in 0..<ligands.count {
            ligands[i].print()
        }
        for i in 0..<flex.count {
            flex[i].print()
        }
    }

}

typealias output_container = ContiguousArray<OutputType>

struct OutputType: Comparable {
    var c: Conf
    var e: fl
    var lb: fl
    var ub: fl
    var intra: fl
    var inter: fl
    var conf_independent: fl
    var unbound: fl
    var total: fl
    var coords: vecv
    
    init(_ c_: Conf, _ e_: fl) {
        self.c = c_
        self.e = e_
        self.lb = 0
        self.ub = 0
        self.intra = 0
        self.inter = 0
        self.conf_independent = 0
        self.unbound = 0
        self.total = 0
        self.coords = vecv()
    }
    
    static func < (a: OutputType, b: OutputType) -> Bool {
        return a.e < b.e
    }
    
    static func == (lhs: OutputType, rhs: OutputType) -> Bool {
        lhs.e == rhs.e && lhs.total == rhs.total
    }

}
