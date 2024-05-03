//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation
import simd

struct AtomIndex: Equatable {
    var i: sz
    var in_grid: Bool
    
    init(i: sz, in_grid: Bool) {
        self.i = i
        self.in_grid = in_grid
    }
    
    init() {
        self.i = sz.max
        self.in_grid = false
    }
}

class Bond: NSCopying {
    var connected_atom_index: AtomIndex
    var length: fl
    var rotatable: Bool
    
    init(connected_atom_index: AtomIndex, length: fl, rotatable: Bool) {
        self.connected_atom_index = connected_atom_index
        self.length = length
        self.rotatable = rotatable
    }
        
    func copy(with zone: NSZone? = nil) -> Any {
        let newBond = Bond(connected_atom_index: self.connected_atom_index,
                           length: self.length, rotatable: self.rotatable)
        return newBond
    }
    
}

typealias atomv = ContiguousArray<Atom>

class Atom: AtomBase {
    var coords: vec
    var bonds: ContiguousArray<Bond>
    
    init(coords: vec, bonds: ContiguousArray<Bond>) {
        self.coords = coords
        self.bonds = bonds
        super.init()
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newBonds = ContiguousArray<Bond>(self.bonds.map { $0.copy() as! Bond })
        let newAtom = Atom(coords: self.coords, bonds: newBonds)
        newAtom.charge = self.charge
        newAtom.xs = self.xs
        newAtom.el = self.el
        newAtom.ad = self.ad
        newAtom.sy = self.sy
        return newAtom
    }
}
