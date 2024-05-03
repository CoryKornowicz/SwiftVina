//
//  Appender.swift
//  
//
//  Created by Cory Kornowicz on 11/14/23.
//

import Foundation

class AppenderInfo {
    var grid_atoms_size: sz
    var m_num_movable_atoms: sz
    var atoms_size: sz
    
    init(_ m: Model) {
        self.grid_atoms_size = sz(m.grid_atoms.count)
        self.m_num_movable_atoms = m.m_num_movable_atoms
        self.atoms_size = sz(m.atoms.count)
    }
}

class Appender {
    var a_info: AppenderInfo
    var b_info: AppenderInfo 
    var is_a: Bool

    init(_ a: Model, _ b: Model) {
        self.a_info = AppenderInfo(a)
        self.b_info = AppenderInfo(b)
        self.is_a = true
    }

    func new_grid_index(_ x: sz) -> sz {
        return (is_a ? x : (a_info.grid_atoms_size + x)) // a-grid_atoms spliced before b-grid_atoms
    }

    func transform_coord_index(_ x: sz) -> sz {
        if(is_a) {
            if(x < a_info.m_num_movable_atoms)  { return x } // a-movable unchanged
            else                                { return x + b_info.m_num_movable_atoms } // b-movable spliced before a-inflex
        }
        else {
            if(x < b_info.m_num_movable_atoms)  { return x + a_info.m_num_movable_atoms } // a-movable spliced before b-movable
            else                                { return x + a_info.atoms_size } // all a's spliced before b-inflex
        }
    }

    func transform_atom_index(_ x: AtomIndex) -> AtomIndex {
        var tmp = x
        if(tmp.in_grid) { tmp.i = new_grid_index(tmp.i) }
        else            { tmp.i = transform_coord_index(tmp.i) }
        return tmp
    }
    
    @inlinable
    func update(_ ip: inout InteractingPair) {
        ip.a = transform_coord_index(ip.a)
        ip.b = transform_coord_index(ip.b)
    }
    
    func append(_ a: inout [InteractingPair], _ b: [InteractingPair]) {
        let a_sz = a.count
        a.append(contentsOf: b)
        
        is_a = true
        for i in 0..<a_sz {
            update(&a[i])
        }
        
        is_a = false
        for i in a_sz..<a.count {
            update(&a[i])
        }
    }
    
    func update(_ vec: inout vec) { return } // coordinates & forces - do nothing
    
    func update(_ lig: inout Ligand) {
        lig.transform(self.transform_coord_index)
        transform_ranges(&lig, self.transform_coord_index)
        
        for var i in lig.pairs {
            self.update(&i)
        }
        
        for var i in lig.cont {
            self.update(&i)
        }
    }
    
    func append(_ a: inout [Ligand], _ b: [Ligand]) {
        let a_sz = a.count
        a.append(contentsOf: b)
        
        is_a = true
        for i in 0..<a_sz {
            update(&a[i])
        }
        
        is_a = false
        for i in a_sz..<a.count {
            update(&a[i])
        }
    }

    func update(_ res: inout Residue) {
        transform_ranges(&res, self.transform_coord_index)
    }
    
    func append(_ a: inout [Residue], _ b: [Residue]) {
        let a_sz = a.count
        a.append(contentsOf: b)
        
        is_a = true
        for i in 0..<a_sz {
            update(&a[i])
        }
        
        is_a = false
        for i in a_sz..<a.count {
            update(&a[i])
        }
    }
    
    func update(_ p: inout ParsedLine) {
        if(p.1 != nil) { p.1 = transform_coord_index(p.1!) }
    }
    
    func append(_ a: inout [ParsedLine], _ b: [ParsedLine]) {
        let a_sz = a.count
        a.append(contentsOf: b)
        
        is_a = true
        for i in 0..<a_sz {
            update(&a[i])
        }
        
        is_a = false
        for i in a_sz..<a.count {
            update(&a[i])
        }
    }
    
    func update(_ a: inout Atom) {
        for i in 0..<a.bonds.count {
            a.bonds[i].connected_atom_index = transform_atom_index(a.bonds[i].connected_atom_index)
        }
    }
    
    func append(_ a: inout ContiguousArray<Atom>, _ b: ContiguousArray<Atom>) {
        let a_sz = a.count
        a.append(contentsOf: b)
        
        is_a = true
        for i in 0..<a_sz {
            update(&a[i])
        }
        
        is_a = false
        for i in a_sz..<a.count {
            update(&a[i])
        }
    }
    
    //  coords, minus_forces, atoms
    func coords_append(_ a: inout vecv, _ b: vecv) {
        var b_copy = b // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

        is_a = true
        for i in 0..<a_info.atoms_size {
            update(&a[i])
        }

        is_a = false
        for i in 0..<b_info.atoms_size {
            update(&b_copy[i])
        }

        // interleave
        let b1 = b_copy[..<Int(b_info.m_num_movable_atoms)]
        let b2 = b_copy[Int(b_info.m_num_movable_atoms)...]

        a.insert(contentsOf: b1, at: Int(a_info.m_num_movable_atoms))
        a.append(contentsOf: b2)
    }   
    
    func coords_append(_ a: inout atomv, _ b: atomv) {
        var b_copy = b // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

        is_a = true
        for i in 0..<a_info.atoms_size {
            update(&a[i])
        }

        is_a = false
        for i in 0..<b_info.atoms_size {
            update(&b_copy[i])
        }

        // interleave
        let b1 = b_copy[..<Int(b_info.m_num_movable_atoms)]
        let b2 = b_copy[Int(b_info.m_num_movable_atoms)...]

        a.insert(contentsOf: b1, at: Int(a_info.m_num_movable_atoms))
        a.append(contentsOf: b2)
    }
    
    
}

