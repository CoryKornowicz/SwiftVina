//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation
import simd

class Frame: NSCopying {
    
    var origin: vec
    var orientation_q: qt
    var orientation_m: mat
   
    init(origin_: vec) {
        self.origin = origin_
        self.orientation_q = qt_identity
        self.orientation_m = quaternion_to_r3(q: qt_identity)
    }
    
    @inline(__always)
    func set_orientation(_ q: qt) {
        orientation_q = q
        orientation_m = quaternion_to_r3(q: orientation_q)
    }

    @inline(__always)
    func local_to_lab(_ local_coords: vec) -> vec {
        return origin + orientation_m * local_coords
    }

    @inline(__always)
    func local_to_lab_direction(_ local_direction: vec) -> vec {
        return (orientation_m * local_direction)
    }

    @inline(__always)
    func orientation() -> qt {
        return orientation_q
    }

    @inline(__always)
    func get_origin() -> vec {
        return origin
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let newFrame = Frame(origin_: self.origin)
        newFrame.orientation_q = self.orientation_q
        newFrame.orientation_m = self.orientation_m
        return newFrame
    }
    
}

protocol Transformable {
    typealias F = (sz) -> sz
    mutating func transform(_ f: F)
}

protocol AtomRange: Transformable {
    var begin: sz { get set }
    var end: sz { get set }
}

extension AtomRange {
    mutating func transform(_ f: F) {
        let diff = end - begin
        begin = f(begin)
        end = begin + diff
    }
}

class AtomFrame: Frame, AtomRange {
    var begin: sz = 0
    var end: sz = 0
    
    init(origin_: vec, begin: sz, end: sz) {
        super.init(origin_: origin_)
        self.begin = begin
        self.end = end
    }
    
    func setCoords(_ atoms: atomv, _ coords: inout vecv) {
        for i in begin..<end {
            coords[i] = local_to_lab(atoms[i].coords)
        }
    }
    
    func sumForceAndTorque(coords: vecv, forces: vecv) -> vecp {
        var tmp: vecp = Pair<vec, vec>(zero_vec, zero_vec)
        for i in begin..<end {
            tmp.0 += forces[i]
            tmp.1 += cross_product(coords[i] - origin, forces[i])
        }
        return tmp
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newAtomFrame = AtomFrame(origin_: self.origin, begin: self.begin, end: self.end)
        newAtomFrame.orientation_q = self.orientation_q
        newAtomFrame.orientation_m = self.orientation_m
        return newAtomFrame
    }
}

protocol TorsionCountable {
    func count_torsions(_ s: inout sz)
}

class RigidBody: AtomFrame {
        
    init(origin_: vec, begin_: sz, end_: sz) {
        super.init(origin_: origin_, begin: begin_, end: end_)
        self.orientation_q = qt_identity
        self.orientation_m = quaternion_to_r3(q: self.orientation_q)
    }
    
    func set_conf(_ atoms: atomv, _ coords: inout vecv, _ c: RigidConf) {
        origin = c.position
        set_orientation(c.orientation)
        setCoords(atoms, &coords)
    }

    func set_derivative(_ force_torque: vecp, _ c: inout RigidChange) {
        c.position = force_torque.0
        c.orientation = force_torque.1
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newAtomFrame = RigidBody(origin_: self.origin, begin_: self.begin, end_: self.end)
        newAtomFrame.orientation_q = self.orientation_q
        newAtomFrame.orientation_m = self.orientation_m
        return newAtomFrame
    }
}

extension RigidBody: TorsionCountable {
    func count_torsions(_ s: inout sz) { } // do nothing
}

class AxisFrame: AtomFrame {
    
    var axis: vec = zero_vec
    
    init(origin_: vec, begin_: sz, end_: sz, axis_root: vec) {
        
        let diff = origin_ - axis_root
        let nrm = diff.norm()
        assert(nrm >= epsilon_fl)
        self.axis = (1 / nrm) * diff
        
        super.init(origin_: origin_, begin: begin_, end: end_)
    }
    
    private override init(origin_: vec, begin: sz, end: sz) {
        super.init(origin_: origin_, begin: begin, end: end)
    }

    func set_derivative(_ force_torque: vecp, _ c: inout fl) {
        c = force_torque.1 * axis
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newAxisFrame = AxisFrame(origin_: self.origin, begin: self.begin, end: self.end)
        newAxisFrame.axis = self.axis
        newAxisFrame.orientation_q = self.orientation_q
        newAxisFrame.orientation_m = self.orientation_m
        return newAxisFrame
    }
    
}

class Segment: AxisFrame, TorsionCountable {
    
    var relative_axis: vec = zero_vec
    var relative_origin: vec = zero_vec
    
    init(origin_: vec, begin_: sz, end_: sz, axis_root: vec, parent: Frame) {
        assert(eq(parent.orientation(), qt_identity))
        super.init(origin_: origin_, begin_: begin_, end_: end_, axis_root: axis_root) // will set self.axis in the process!
        relative_axis = axis
        relative_origin = origin - parent.get_origin()
    }
    
    // This function creates a local instance for copying.
    // Considering the axis_root is a zero_vec, the axis will have to be set in the copying function
    private init(origin_: vec, begin: sz, end: sz) {
        super.init(origin_: origin_, begin_: begin, end_: end, axis_root: zero_vec)
    }

    func set_conf(_ parent: Frame, _ atoms: atomv, _ coords: inout vecv, _ c: inout Iterator<fl>) {
        guard let torsion = c.next() else {
            fatalError("Unable to get torsion value")
        }
        origin = parent.local_to_lab(relative_origin)
        axis = parent.local_to_lab_direction(relative_axis)
        var tmp = angle_to_quaternion(axis: axis, angle: torsion) * parent.orientation()
        quaternion_normalize_approx(&tmp)
        set_orientation(tmp)
        setCoords(atoms, &coords)
    }

    func count_torsions(_ s: inout sz) {
        s += 1
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newAxisFrame = Segment(origin_: self.origin, begin: self.begin, end: self.end)
        newAxisFrame.relative_axis = self.relative_axis
        newAxisFrame.relative_origin = self.relative_origin
        newAxisFrame.axis = self.axis
        newAxisFrame.orientation_q = self.orientation_q
        newAxisFrame.orientation_m = self.orientation_m
        return newAxisFrame
    }
    
}

class FirstSegment: AxisFrame, TorsionCountable {
        
    func set_conf(_ atoms: atomv, _ coords: inout vecv, _ torsion: fl) {
        set_orientation(angle_to_quaternion(axis: self.axis, angle: torsion))
        setCoords(atoms, &coords)
    }
    
    func count_torsions(_ s: inout sz) {
        s += 1
    }
}

typealias Branch = Tree<Segment>
typealias Branches = [Branch]

func branches_set_conf(_ b: inout Branches, _ parent: Frame, _ atoms: atomv, _ coords: inout vecv, _ c: inout Iterator<fl>) {
    for i in 0..<b.count {
        b[i].set_conf(parent, atoms, &coords, &c)
    }
}

func branches_derivative(_ b: Branches, _ origin: vec, _ coords: vecv, _ forces: vecv, _ out: inout vecp, _ d: inout Iterator<fl>) {
    for i in 0..<b.count {
        let force_torque = b[i].derivative(coords, forces, &d)
        out.0 += force_torque.0
        let r = b[i].node.get_origin() - origin
        out.1 += cross_product(r, force_torque.0) + force_torque.1
    }
}

protocol NodeBearing<Node> {
    associatedtype Node where Node: TorsionCountable, Node: Transformable, Node: AtomFrame
    var node: Node { get set }
    var children: Branches { get set }
}

class Tree<T: AtomFrame & TorsionCountable & NSCopying>: NodeBearing, NSCopying {
    
    var node: T
    var children: Branches
    
    init(node: AtomFrame & TorsionCountable, children: Branches) {
        self.node = node as! T
        self.children = children
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let nodeCopy = node.copy() as! T
        let childrenCopies = self.children.map { $0.copy() as! Tree<Segment> }
        let newTree = Tree<T>(node: nodeCopy, children: childrenCopies)
        return newTree
    }
}

extension Tree where T == Segment {
    
    func set_conf (_ parent: Frame, _ atoms: atomv, _ coords: inout vecv, _ c: inout Iterator<fl>) {
        node.set_conf(parent, atoms, &coords, &c)
        branches_set_conf(&children, node, atoms, &coords, &c)
    }
    
    // MARK: Not entirely sure why we are mutating `d` given that it is not used further in the scope??
    func derivative(_ coords: vecv, _ forces: vecv, _ p: inout Iterator<fl>) -> vecp {
        var force_torque: vecp = node.sumForceAndTorque(coords: coords, forces: forces)
        guard var d = p.next() else {
            fatalError("No more torsion values")
        }
        branches_derivative(children, node.get_origin(), coords, forces, &force_torque, &p)
        node.set_derivative(force_torque, &d)
        return force_torque
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let nodeCopy = node.copy() as! T
        let childrenCopies = self.children.map { $0.copy() as! Tree<Segment> }
        let newTree = Tree<T>(node: nodeCopy, children: childrenCopies)
        return newTree
    }
}


class Heterotree <T: AtomFrame & TorsionCountable & NSCopying>: NodeBearing {
    
    var node: T
    var children: Branches
    
    init(node: AtomFrame & TorsionCountable, children: Branches) {
        self.node = node as! T
        self.children = children
    }

}

typealias MainBranch = Heterotree<FirstSegment>
typealias FlexibleBody = Heterotree<RigidBody>

extension Heterotree where T == FirstSegment {

    func set_conf(_ atoms: atomv, _ coords: inout vecv, _ c: ResidueConf) {
        var p = Iterator<fl>(c.torsions)
        guard let startValue = p.next() else {
            fatalError("Cannot get first torsion value")
        }
        node.set_conf(atoms, &coords, startValue)
        branches_set_conf(&children, node, atoms, &coords, &p)
    }

    func derivative(_ coords: vecv, _ forces: vecv, _ c: inout ResidueChange) {
        var force_torque = node.sumForceAndTorque(coords: coords, forces: forces)
        var p = Iterator<fl>(c.torsions)
        guard var d = p.next() else {
            fatalError("Cannot get first torsion value")
        }
        branches_derivative(children, node.get_origin(), coords, forces, &force_torque, &p)
        node.set_derivative(force_torque, &d)
    }
    
}

extension Heterotree where T == RigidBody {
    
    func set_conf(_ atoms: atomv, _ coords: inout vecv, _ c: LigandConf) {
        node.set_conf(atoms, &coords, c.rigid)
        var p = Iterator<fl>(c.torsions)
        branches_set_conf(&children, node, atoms, &coords, &p)
    }

    func derivative(_ coords: vecv, _ forces: vecv, _ c: inout LigandChange) {
        var force_torque: vecp = node.sumForceAndTorque(coords: coords, forces: forces)
        var p = Iterator<fl>(c.torsions)
        branches_derivative(children, node.get_origin(), coords, forces, &force_torque, &p)
        node.set_derivative(force_torque, &c.rigid)
    }
    
}

func count_torsions(_ t: some NodeBearing, _ s: inout sz) {
    t.node.count_torsions(&s)
    for i in 0..<t.children.count {
        count_torsions(t.children[i], &s)
    }
}


final class VectorMutable<T, C, D>: NSCopying where T: NSCopying {
    
    var elements: [T] = [T]()
    
    var count: Int {
        return elements.count
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let newVM = VectorMutable<T, C, D>()
        newVM.elements = self.elements.map { $0.copy() as! T }
        return newVM
    }
    
}

extension VectorMutable where T == Ligand, C == LigandConf, D == LigandChange {
    
    func count_torsions() -> szv {
        var tmp = ContiguousArray<sz>(repeating: 0, count: elements.count)
        for i in 0..<elements.count {
            SwiftVina.count_torsions(elements[i], &tmp[i])
        }
        return tmp
    }
    
    func derivative(coords: vecv, forces: vecv, c: inout [D]) { // D == ligand_change || residue_change
        for i in 0..<elements.count {
            elements[i].derivative(coords, forces, &c[i])
        }
    }
    
    func set_conf(atoms: atomv, coords: inout vecv, c: [C]) { // C == ligand_conf || residue_conf
        for i in 0..<elements.count {
            elements[i].set_conf(atoms, &coords, c[i])
        }
    }
}

extension VectorMutable where T == Residue, C == ResidueConf, D == ResidueChange {
    
    func count_torsions() -> szv {
        var tmp = ContiguousArray<sz>(repeating: 0, count: elements.count)
        for i in 0..<elements.count {
            SwiftVina.count_torsions(elements[i], &tmp[i])
        }
        return tmp
    }
    
    func derivative(coords: vecv, forces: vecv, c: inout [D]) { // D == ligand_change || residue_change
        for i in 0..<elements.count {
            elements[i].derivative(coords, forces, &c[i])
        }
    }
    
    func set_conf(atoms: atomv, coords: inout vecv, c: [C]) { // C == ligand_conf || residue_conf
        for i in 0..<elements.count {
            elements[i].set_conf(atoms, &coords, c[i])
        }
    }
}



func transform_ranges(_ t: inout some NodeBearing, _ f: (sz) -> sz) {
    t.node.transform(f)
    for i in 0..<t.children.count {
        transform_ranges(&t.children[i], f)
    }
}


