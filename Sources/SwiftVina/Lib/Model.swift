//
//  File.swift
//
//
//  Created by Cory Kornowicz on 11/11/23.
//

import Foundation
import simd

typealias InteractingPairs = Array<InteractingPair>
typealias ParsedLine = Pair<String, sz?>
typealias Context = Array<ParsedLine>

typealias DistanceTypeMatrix = strictly_triangular_matrix<DistanceType>

enum DistanceType {
    case DISTANCE_FIXED
    case DISTANCE_ROTOR
    case DISTANCE_VARIABLE
}

struct InteractingPair {
    var type_pair_index: sz
    var a: sz
    var b: sz
}

struct BranchMetrics {
    var length: sz
    var corner2corner: sz
    init() {
        self.length = 0
        self.corner2corner = 0
    }
}

func get_atom_range(_ t: some NodeBearing) -> AtomRange {
    var tmp: AtomRange = t.node
    for i in t.children {
        let r: AtomRange = get_atom_range(i)
        if(tmp.begin > r.begin) {
            tmp.begin = r.begin
        }
        if(tmp.end < r.end  ) {
            tmp.end = r.end
        }
    }
    return tmp
}

func get_branch_metrics(_ t: some NodeBearing) -> BranchMetrics {
    var tmp: BranchMetrics = BranchMetrics()
    
    if (!t.children.isEmpty) {
        var corner2corner_max: sz = 0
        var lengths: szv = []
        
        for i in t.children {
            let res: BranchMetrics = get_branch_metrics(i)
            if corner2corner_max < res.corner2corner {
                corner2corner_max = res.corner2corner
            }
            lengths.append(res.length + 1)
        }
        
        lengths.sort()
        tmp.length = lengths.last!
        tmp.corner2corner = tmp.length
        
        if lengths.count >= 2 {
            tmp.corner2corner += lengths[lengths.count - 1]
        }
        
        if tmp.corner2corner < corner2corner_max {
            tmp.corner2corner = corner2corner_max
        }
    }
    
    return tmp
}

final class Ligand: FlexibleBody, AtomRange, NSCopying {
    
    var begin: sz = 0
    var end: sz = 0
    var degrees_of_freedom: UInt = 0
    var pairs: InteractingPairs = []
    var cont: Context = []
    
    init(f: FlexibleBody, degrees_of_freedom: UInt) {
        super.init(node: f.node, children: f.children)
        self.degrees_of_freedom = degrees_of_freedom
    }
    
    func set_range() {
        let tmp: AtomRange = get_atom_range(self)
        begin = tmp.begin
        end   = tmp.end
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let nodeCopy = self.node.copy() as! RigidBody
        let childrenCopies = self.children.map { $0.copy() } as! Branches
        let newLigand = Ligand(f: FlexibleBody(node: nodeCopy, children: childrenCopies), degrees_of_freedom: self.degrees_of_freedom)
        newLigand.begin = self.begin
        newLigand.end = self.end
        newLigand.pairs = self.pairs
        newLigand.cont = self.cont
        return newLigand
    }
    
}

final class Residue: MainBranch, NSCopying {
    
    func copy(with zone: NSZone? = nil) -> Any {
        let nodeCopy = self.node.copy() as! FirstSegment
        let childrenCopies = self.children.map { $0.copy() } as! Branches
        let newResidue = Residue(node: nodeCopy, children: childrenCopies)
        return newResidue
    }
    
}

final class Model: NSCopying {
    
    var atoms: atomv = ContiguousArray<Atom>() // movable, inflex
    var coords: vecv = []
    var minus_forces: vecv = []
    var grid_atoms: atomv = ContiguousArray<Atom>()
    var m_num_movable_atoms: sz = 0
    var ligands: VectorMutable<Ligand, LigandConf, LigandChange> = VectorMutable<Ligand, LigandConf, LigandChange>()
    var flex: VectorMutable<Residue, ResidueConf, ResidueChange> = VectorMutable<Residue, ResidueConf, ResidueChange>()
    var flex_context: Context = Context()
    var other_pairs: InteractingPairs = InteractingPairs() // INTRAmolecular interactions: flex_i - flex_j and flex_i - flex_i
    var inter_pairs: InteractingPairs = InteractingPairs() // INTERmolecular interactions: ligand - flex and ligand_i - ligand_j
    var glue_pairs: InteractingPairs = InteractingPairs() // INTRAmolecular interactions: glue_i - glue_i
    
    let m_atom_typing_used: AtomType.T
    
    var atom_typing_used: AtomType.T {
        self.m_atom_typing_used
    }
    
    var num_movable_atoms: sz {
        return m_num_movable_atoms
    }
    
    var num_internal_pairs: sz {
        var tmp = 0
        for i in 0..<ligands.count {
            tmp += ligands.elements[i].pairs.count
        }
        return sz(tmp)
    }
    
    var num_other_pairs: sz {
        return sz(other_pairs.count)
    }
    
    var num_ligands: sz {
        return sz(ligands.count)
    }
    
    var num_flex: sz {
        return sz(flex.count)
    }
    
    var num_atoms: sz {
        return sz(atoms.count)
    }
    
    init() {
        self.m_atom_typing_used = .XS
    }
    
    init(_ atype: AtomType.T) {
        self.m_atom_typing_used = atype
    }

    func copy(with zone: NSZone? = nil) -> Any {
        let copy = Model(self.m_atom_typing_used)
        
        copy.atoms = ContiguousArray(self.atoms.map({ atom in
            atom.copy() as! Atom
        }))
        
        copy.coords = self.coords
        copy.minus_forces = self.minus_forces
        
        copy.grid_atoms = ContiguousArray(self.grid_atoms.map({ grid_atom in
            grid_atom.copy() as! Atom
        }))
        
        copy.m_num_movable_atoms = self.m_num_movable_atoms
        copy.ligands = self.ligands.copy() as! VectorMutable<Ligand, LigandConf, LigandChange>
        copy.flex = self.flex.copy() as! VectorMutable<Residue, ResidueConf, ResidueChange>        
        copy.flex_context = self.flex_context
        copy.other_pairs = self.other_pairs
        copy.inter_pairs = self.inter_pairs
        copy.glue_pairs = self.glue_pairs
        return copy
    }
    
    func initialize(mobility: DistanceTypeMatrix) throws {
        for i in 0..<ligands.count {
            ligands.elements[i].set_range()
        }
        try assign_bonds(mobility: mobility)
        assign_types()
        initialize_pairs(mobility)
    }
    
    func initialize_pairs(_ mobility: DistanceTypeMatrix) {
        /* Interactions:
         - ligand_i - ligand_i : YES (1-4 only) (ligand.pairs)
         - flex_i   - flex_i   : YES (1-4 only) (other_pairs)
         - flex_i   - flex_j   : YES (other_pairs)
         - intra macrocycle closure interactions: NO (1-2, 1-3, 1-4)
         */
        for i in 0..<sz(atoms.count) {
            let i_lig: sz = find_ligand(i)
            let bonded_atoms: szv = bonded_to(i, 3) // up to 1-4
            
            for j in i + 1..<sz(atoms.count) {
                if mobility[i, j] == .DISTANCE_VARIABLE && !bonded_atoms.contains(j) {
                    if is_closure_clash(i, j) || is_unmatched_closure_dummy(i, j) { continue }
                    
                    let t1: sz = atoms[i].get(atom_typing_used)
                    let t2: sz = atoms[j].get(atom_typing_used)
                    let n: sz = num_atom_types(atom_typing_used)
                    
                    if t1 < n && t2 < n { // exclude, say, Hydrogens
                        let type_pair_index: sz = triangular_matrix_index_permissive(n, t1, t2)
                        let ip: InteractingPair = InteractingPair(type_pair_index: type_pair_index, a: i, b: j)
                        
                        if is_glue_pair(i, j) {
                            // Add glue_i - glue_i interaction pair
                            glue_pairs.append(ip)
                        } else if i_lig < ligands.count && find_ligand(j) == i_lig {
                            // Add INTRAmolecular ligand_i - ligand_i
                            ligands.elements[i_lig].pairs.append(ip)
                        } else if !is_atom_in_ligand(i) && !is_atom_in_ligand(j) {
                            // Add INTRAmolecular flex_i  - flex_i but also flex_i - flex_j
                            other_pairs.append(ip)
                        }
                    }
                }
            }
        }
    }
    
    func get_atoms() -> atomv {
        return atoms
    }
    
    func get_atom(_ i: sz) -> Atom {
        return atoms[i]
    }
    
    func get_atom(_ i: AtomIndex) -> Atom {
        return i.in_grid ? grid_atoms[i.i] : atoms[i.i]
    }
    
    func get_ligand(_ i: sz) -> Ligand {
        return self.ligands.elements[i]
    }
    
    func get_coords(_ i: sz) -> vec {
        return self.coords[i]
    }
    
    func ligand_degrees_of_freedom(_ ligand_number: sz) -> sz {
        return sz(ligands.elements[ligand_number].degrees_of_freedom)
    }
    
    func ligand_longest_branch(_ ligand_number: sz) -> sz {
        return get_branch_metrics(ligands.elements[ligand_number]).length
    }
    
    func ligand_length(_ ligand_number: sz) -> sz {
        return get_branch_metrics(ligands.elements[ligand_number]).corner2corner
    }
    
    func find_ligand(_ a: sz) -> sz {
        for i in 0..<ligands.count {
            if a >= ligands.elements[i].begin && a < ligands.elements[i].end {
                return sz(i)
            }
        }
        return sz(ligands.count)
    }
    
    func is_atom_in_ligand(_ a: sz) -> Bool {
        for i in 0..<ligands.count {
            if a >= ligands.elements[i].begin && a < ligands.elements[i].end {
                return true
            }
        }
        return false
    }
    
    func get_movable_atom_types(_ atom_typing_used_: AtomType.T) -> szv {
        var tmp: szv = []
        let n: sz = num_atom_types(atom_typing_used_)
        
        for i in 0..<m_num_movable_atoms {
            let a: Atom = atoms[i]
            let t: sz = a.get(atom_typing_used_)
            if t < n && !tmp.contains(t) {
                tmp.append(t)
            }
        }
        return tmp
    }
    
    func get_size() -> ConfSize {
        let tmp: ConfSize = ConfSize()
        tmp.ligands = ligands.count_torsions()
        tmp.flex    = flex.count_torsions()
        return tmp
    }
    
    func get_initial_conf() -> Conf { // torsions = 0, orientations = identity, ligand positions = current
        let cs: ConfSize = get_size()
        var tmp: Conf = Conf(cs)
        tmp.set_to_null()
        for i in 0..<ligands.count {
            tmp.ligands[i].rigid.position = ligands.elements[i].node.get_origin()
        }
        return tmp
    }
    
//    func get_ligand_coords() -> vecv { // FIXME?
//        var tmp: vecv = []
//        let lig: Ligand = ligands[0]
//        for i in lig.begin..<lig.end {
//            tmp.append(coords[i])
//        }
//        return tmp
//    }
    
    func get_ligand_coords() -> Array<fl> {
        var tmp: Array<fl> = []
        let lig: Ligand = ligands.elements[0]
        for i in lig.begin..<lig.end {
            tmp.append(coords[i][0])
            tmp.append(coords[i][1])
            tmp.append(coords[i][2])
        }
        return tmp
    }
    
    func get_heavy_atom_movable_coords() -> vecv { // FIXME mv
        var tmp: vecv = []
        for i in 0..<num_movable_atoms {
            if atoms[i].el != EL_TYPE_H {
                tmp.append(coords[i])
            }
        }
        return tmp
    }
    
    func is_movable_atom(_ a: sz) -> Bool {
        if (a < num_movable_atoms) {
            return true
        } else {
            return false
        }
    }
    
    func append(_ model: Model) {
        
        assert(atom_typing_used == model.atom_typing_used)
        
        let t: Appender = Appender(self, model)
        
        t.append(&other_pairs, model.other_pairs)
        t.append(&inter_pairs, model.inter_pairs)
        t.append(&glue_pairs, model.glue_pairs)
        
        assert(minus_forces.count == coords.count)
        assert(model.minus_forces.count == model.coords.count)
        
        t.coords_append(&coords, model.coords)
        t.coords_append(&minus_forces, model.minus_forces)
        
        t.append(&ligands.elements, model.ligands.elements)
        t.append(&flex.elements, model.flex.elements)
        t.append(&flex_context, model.flex_context)
        
        t.append(&grid_atoms, model.grid_atoms)
        t.coords_append(&atoms, model.atoms)
        
        // Add interaction pairs between previously added atoms and (very likely) ligand atoms
        /* Interactions:
         - flex     - ligand   : YES (inter_pairs)
         - flex_i   - flex_j   : YES (other_pairs) but append is used mostly for adding ligand
         - ligand_i - ligand_j : YES (inter_pairs)
         - macrocycle closure interactions: NO (1-2, 1-3, 1-4)
         */
        
        for i in 0..<m_num_movable_atoms {
            for j in m_num_movable_atoms..<m_num_movable_atoms + model.m_num_movable_atoms {
                if is_closure_clash(i, j) || is_unmatched_closure_dummy(i, j) {
                    continue
                }
                
                let a: Atom = atoms[i]
                let b: Atom = atoms[j]
                
                let t1: sz = a.get(atom_typing_used)
                let t2: sz = b.get(atom_typing_used)
                let n: sz = num_atom_types(atom_typing_used)
                
                if t1 < n && t2 < n {
                    let type_pair_index: sz = triangular_matrix_index_permissive(n, t1, t2)
                    
                    if is_glue_pair(i, j) {
                        glue_pairs.append(InteractingPair(type_pair_index: type_pair_index, a: i, b: j)) // glue_i - glue_j
                    } else if is_atom_in_ligand(i) && is_atom_in_ligand(j) {
                        inter_pairs.append(InteractingPair(type_pair_index: type_pair_index, a: i, b: j)) // INTER: ligand_i - ligand_j
                    } else if is_atom_in_ligand(i) || is_atom_in_ligand(j) {
                        inter_pairs.append(InteractingPair(type_pair_index: type_pair_index, a: i, b: j)) // INTER: flex - ligand
                    } else {
                        other_pairs.append(InteractingPair(type_pair_index: type_pair_index, a: i, b: j)) // INTRA: flex_i - flex_j
                    }
                }
            }
        }
        
        m_num_movable_atoms += model.m_num_movable_atoms
    }
    
    func sz_to_atom_index(_ i: sz) -> AtomIndex {
        if i < grid_atoms.count { return AtomIndex(i: i                       , in_grid:  true)}
        else                    { return AtomIndex(i: i - sz(grid_atoms.count), in_grid: false)}
    }
    
    func distance_type_between(_ mobility: DistanceTypeMatrix, _ i: AtomIndex, _ j: AtomIndex) -> DistanceType {
        if i.in_grid && j.in_grid { return .DISTANCE_FIXED }
        if i.in_grid { return (j.i < m_num_movable_atoms) ? .DISTANCE_VARIABLE : .DISTANCE_FIXED }
        if j.in_grid { return (i.i < m_num_movable_atoms) ? .DISTANCE_VARIABLE : .DISTANCE_FIXED }
        assert(!i.in_grid)
        assert(!j.in_grid)
        assert(i.i < atoms.count)
        assert(j.i < atoms.count)
        let a: sz = i.i
        let b: sz = j.i
        if a == b { return .DISTANCE_FIXED }
        return (a < b) ? mobility[a, b] : mobility[b, a]
    }
    
    @inlinable
    func atom_coords(_ i: AtomIndex) -> vec {
        return i.in_grid ? grid_atoms[i.i].coords : coords[i.i]
    }
    
    @inline(__always)
    func distance_sqr_between(_ a: AtomIndex, _ b: AtomIndex) -> fl {
        return vec_distance_sqr(atom_coords(a), atom_coords(b))
    }
    
    func atom_exists_between(_ mobility: DistanceTypeMatrix, _ a: AtomIndex, _ b: AtomIndex, _ relevant_atoms: szv) -> Bool {
        // there is an atom closer to both a and b then they are to each other and immobile relative to them
        let r2: fl = distance_sqr_between(a, b)
        for i in 0..<relevant_atoms.count {
            let c: AtomIndex = sz_to_atom_index(relevant_atoms[i])
            if a == c || b == c { continue }
            let ac: DistanceType = distance_type_between(mobility, a, c)
            let bc: DistanceType = distance_type_between(mobility, b, c)
            if ac != .DISTANCE_VARIABLE && bc != .DISTANCE_VARIABLE && distance_sqr_between(a, c) < r2 && distance_sqr_between(b, c) < r2 {
                return true
            }
        }
        return false
    }
    
    final class Beads {
        
        var radius_sqr: fl
        var data: Array<Pair<vec, szv>> = [Pair<vec, szv>]()
        
        init(radius_sqr: fl) {
            self.radius_sqr = radius_sqr
        }
        
        func add(_ index: sz, coords: vec) {
            for i in 0..<data.count {
                if vec_distance_sqr(coords, data[i].0) < radius_sqr {
                    data[i].1.append(index)
                    return
                }
            }
            // not found
            let tmp: Pair<vec, szv> = Pair<vec, szv>(coords, [index])
            data.append(tmp)
        }
    }
    
    
    func assign_bonds(mobility: DistanceTypeMatrix) throws {
        let bond_length_allowance_factor: fl = 1.1
        let n: sz = sz(grid_atoms.count + atoms.count)
        // construct beads
        let bead_radius: fl = 15
        let beads_instance: Beads = Beads(radius_sqr: sqr(bead_radius))
        for i in 0..<n {
            let i_atom_index: AtomIndex = sz_to_atom_index(i)
            beads_instance.add(i, coords: atom_coords(i_atom_index))
        }

        // assign bonds
        for i in 0..<n {
            let i_atom_index: AtomIndex = sz_to_atom_index(i)
            let i_atom_coords: vec = atom_coords(i_atom_index)
            let i_atom: Atom = get_atom(i_atom_index)
            let max_covalent_r: fl = max_covalent_radius() // FIXME mv to atom_constants
            var i_atom_covalent_radius: fl = max_covalent_r
            if i_atom.ad < AD_TYPE_SIZE {
                i_atom_covalent_radius = ad_type_property(i: i_atom.ad).covalent_radius
            }
            
            //find relevant atoms
            var relevant_atoms: szv = []
            let bead_cutoff_sqr: fl = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r))
            for b in 0..<beads_instance.data.count {
                if vec_distance_sqr(beads_instance.data[b].0, i_atom_coords) > bead_cutoff_sqr { continue }
                let bead_elements: szv = beads_instance.data[b].1
                for bead_elements_i in 0..<bead_elements.count {
                    let j = bead_elements[bead_elements_i]
                    let j_atom_index: AtomIndex = sz_to_atom_index(j)
                    // TODO: Fix why bond_length was removed
                    // var j_atom: Atom = atoms[unsafe: j_atom_index.i]
                    // let bond_length: fl = i_atom.optimal_covalent_bond_length(j_atom)
                    let dt: DistanceType = distance_type_between(mobility, i_atom_index, j_atom_index)
                    if dt != .DISTANCE_VARIABLE && i != j {
                        let r2: fl = distance_sqr_between(i_atom_index, j_atom_index)
                        //if(r2 < sqr(bond_length_allowance_factor * bond_length))
                        if r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)) {
                            relevant_atoms.append(j)
                        }
                    }
                }
            }
            
            // find bonded atoms
            for relevant_atoms_i in 0..<relevant_atoms.count {
                let j = relevant_atoms[relevant_atoms_i]
                if j <= i { continue } // already considered
                let j_atom_index: AtomIndex = sz_to_atom_index(j)
                let j_atom: Atom = get_atom(j_atom_index)
                let bond_length: fl = i_atom.optimal_covalent_bond_length(j_atom)
                let dt: DistanceType = distance_type_between(mobility, i_atom_index, j_atom_index)
                let r2: fl = distance_sqr_between(i_atom_index, j_atom_index)
                if r2 < sqr(bond_length_allowance_factor * bond_length) && !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms) {
                    let rotatable: Bool = (dt == .DISTANCE_ROTOR)
                    let length: fl = sqrt(r2)
                    i_atom.bonds.append(Bond(connected_atom_index: j_atom_index, length: length, rotatable: rotatable))
                    j_atom.bonds.append(Bond(connected_atom_index: i_atom_index, length: length, rotatable: rotatable))
                }
            }
        }
    }
    
    
    func bonded_to_HD(_ a: Atom) -> Bool {
        for b in a.bonds {
            if get_atom(b.connected_atom_index).ad == AD_TYPE_HD {
                return true
            }
        }
        return false
    }
    
    func bonded_to_heteroatom(_ a: Atom) -> Bool {
        for b in a.bonds {
            if get_atom(b.connected_atom_index).is_heteroatom() {
                return true
            }
        }
        return false
    }
    
    func assign_types() {
        for i in 0..<sz(grid_atoms.count + atoms.count) {
            let ai: AtomIndex = sz_to_atom_index(i)
            let a: Atom = get_atom(ai)
            a.assign_el()
            
            let acceptor: Bool = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA) // X-Score forumaltion apparently ignores SA
            let donor_NorO: Bool = (a.el == EL_TYPE_Met || bonded_to_HD(a))
            
            switch a.el {
            case EL_TYPE_H    : break
            case EL_TYPE_C    :
                if     (a.ad == AD_TYPE_CG0) { a.xs = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG0 : XS_TYPE_C_H_CG0 }
                else if(a.ad == AD_TYPE_CG1) { a.xs = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG1 : XS_TYPE_C_H_CG1 }
                else if(a.ad == AD_TYPE_CG2) { a.xs = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG2 : XS_TYPE_C_H_CG2 }
                else if(a.ad == AD_TYPE_CG3) { a.xs = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG3 : XS_TYPE_C_H_CG3 }
                else                         { a.xs = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H }
                break
            case EL_TYPE_N    : a.xs = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break
            case EL_TYPE_O    : a.xs = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break
            case EL_TYPE_S    : a.xs = XS_TYPE_S_P; break
            case EL_TYPE_P    : a.xs = XS_TYPE_P_P; break
            case EL_TYPE_F    : a.xs = XS_TYPE_F_H; break
            case EL_TYPE_Cl   : a.xs = XS_TYPE_Cl_H; break
            case EL_TYPE_Br   : a.xs = XS_TYPE_Br_H; break
            case EL_TYPE_I    : a.xs = XS_TYPE_I_H; break
            case EL_TYPE_Si   : a.xs = XS_TYPE_Si; break
            case EL_TYPE_At   : a.xs = XS_TYPE_At; break
            case EL_TYPE_Met  : a.xs = XS_TYPE_Met_D; break
            case EL_TYPE_Dummy:
                if      (a.ad == AD_TYPE_G0) { a.xs = XS_TYPE_G0 }
                else if (a.ad == AD_TYPE_G1) { a.xs = XS_TYPE_G1 }
                else if (a.ad == AD_TYPE_G2) { a.xs = XS_TYPE_G2 }
                else if (a.ad == AD_TYPE_G3) { a.xs = XS_TYPE_G3 }
                else if (a.ad == AD_TYPE_W)  { a.xs = XS_TYPE_SIZE } // no W atoms in XS types
                else { assert(false) }
                break
            case EL_TYPE_SIZE : break
            default: assert(false)
            }
        }
    }
    
    
    // remove 1-2, 1-3 and 1-4 interactions around CG-CG bond
    func is_closure_clash(_ i: sz, _ j: sz) -> Bool {
        let t1: sz = atoms[i].get(.AD)
        let t2: sz = atoms[j].get(.AD)
        
        if ((t1 == AD_TYPE_CG0 && t2 == AD_TYPE_G0) || (t2 == AD_TYPE_CG0 && t1 == AD_TYPE_G0) ||
            (t1 == AD_TYPE_CG1 && t2 == AD_TYPE_G1) || (t2 == AD_TYPE_CG1 && t1 == AD_TYPE_G1) ||
            (t1 == AD_TYPE_CG2 && t2 == AD_TYPE_G2) || (t2 == AD_TYPE_CG2 && t1 == AD_TYPE_G2) ||
            (t1 == AD_TYPE_CG3 && t2 == AD_TYPE_G3) || (t2 == AD_TYPE_CG3 && t1 == AD_TYPE_G3))
        { return false } // it's a G-CG pair, not a clash
        
        let neighbors_of_i: szv = bonded_to(i, 1) // up to 1-2 interactions
        let neighbors_of_j: szv = bonded_to(j, 1) // up to 1-2 interactions
        var i_has_CG0: Bool = false
        var i_has_CG1: Bool = false
        var i_has_CG2: Bool = false
        var i_has_CG3: Bool = false
        
        for index in neighbors_of_i {
            let type_i: sz = atoms[index].get(.AD)
            if      (type_i == AD_TYPE_CG0) { i_has_CG0 = true }
            else if (type_i == AD_TYPE_CG1) { i_has_CG1 = true }
            else if (type_i == AD_TYPE_CG2) { i_has_CG2 = true }
            else if (type_i == AD_TYPE_CG3) { i_has_CG3 = true }
        }
        
        for index in neighbors_of_j {
            let type_j: sz = atoms[index].get(.AD)
            if ((type_j == AD_TYPE_CG0 && i_has_CG0) ||
                (type_j == AD_TYPE_CG1 && i_has_CG1) ||
                (type_j == AD_TYPE_CG2 && i_has_CG2) ||
                (type_j == AD_TYPE_CG3 && i_has_CG3))
            { return true }
        }
        return false
    }
    
    
    func is_glue_pair(_ i: sz, _ j: sz) -> Bool {
        let t1: sz = atoms[i].get(.AD)
        let t2: sz = atoms[j].get(.AD)
        
        if ((t1 == AD_TYPE_CG0 && t2 == AD_TYPE_G0) || (t2 == AD_TYPE_CG0 && t1 == AD_TYPE_G0) ||
            (t1 == AD_TYPE_CG1 && t2 == AD_TYPE_G1) || (t2 == AD_TYPE_CG1 && t1 == AD_TYPE_G1) ||
            (t1 == AD_TYPE_CG2 && t2 == AD_TYPE_G2) || (t2 == AD_TYPE_CG2 && t1 == AD_TYPE_G2) ||
            (t1 == AD_TYPE_CG3 && t2 == AD_TYPE_G3) || (t2 == AD_TYPE_CG3 && t1 == AD_TYPE_G3))
        { return true }
        else
        { return false }
    }
    
    func is_unmatched_closure_dummy(_ i: sz, _ j: sz) -> Bool {
        let t1: sz = atoms[i].get(.AD)
        let t2: sz = atoms[j].get(.AD)
        
        if ((t1 == AD_TYPE_G0 && t2 != AD_TYPE_CG0) || (t2 == AD_TYPE_G0 && t1 != AD_TYPE_CG0) ||
            (t1 == AD_TYPE_G1 && t2 != AD_TYPE_CG1) || (t2 == AD_TYPE_G1 && t1 != AD_TYPE_CG1) ||
            (t1 == AD_TYPE_G2 && t2 != AD_TYPE_CG2) || (t2 == AD_TYPE_G2 && t1 != AD_TYPE_CG2) ||
            (t1 == AD_TYPE_G3 && t2 != AD_TYPE_CG3) || (t2 == AD_TYPE_G3 && t1 != AD_TYPE_CG3))
        { return true }
        else
        { return false }
    }
    
    func bonded_to(_ a: sz, _ n: sz, _ out: inout szv) {
        if !out.contains(a) {
            out.append(a)
            if n > 0 {
                for i in atoms[a].bonds {
                    let b: Bond = i
                    if !b.connected_atom_index.in_grid {
                        bonded_to(b.connected_atom_index.i, n-1, &out)
                    }
                }
            }
        }
    }
    
    func bonded_to(_ a: sz, _ n: sz) -> szv {
        var tmp: szv = []
        bonded_to(a, n, &tmp)
        return tmp
    }
    
    // Sums the movable atoms and returns the center by dividing by the number of movable atoms (uses SIMD)
    func center() -> [fl] {
        var center: vec = zero_vec
        for i in 0..<num_movable_atoms {
            center += coords[i]
        }
        center /= fl(num_movable_atoms)
        return [center[0], center[1], center[2]]
    }
    
    func string_write_coord(_ i: sz, _ x: fl, _ str: inout String) {
        assert(i > 0)
        let out: String = String(format: "%8.3f", x)
        let _i = Int(i - 1)
        let start = str.index(str.startIndex, offsetBy: _i)
        let end = str.index(str.startIndex, offsetBy: _i + 8)
        
        str = str.replacingCharacters(in: start..<end, with: out)
    }
    
    func coords_to_pdbqt_string(_ coords: vec, str: String) -> String {
        var tmp: String = str
        string_write_coord(31, coords[0], &tmp)
        string_write_coord(39, coords[1], &tmp)
        string_write_coord(47, coords[2], &tmp)
        return tmp
    }
    
    func write_context(_ c: Context, _ out: inout TextStreamer) {
        verify_bond_lengths()
        for i in c {
            let str: String = i.0
            if let second = i.1 {
                print(coords_to_pdbqt_string(coords[second], str: str), to: &out)
            } else {
                print(str, to: &out)
            }
        }
    }
    
    func write_context(_ c: Context, _ out: inout String) {
        verify_bond_lengths()
        for i in c {
            let str: String = i.0
            if str.isEmpty { continue } // skip empty lines
            if let second = i.1 {
                out += coords_to_pdbqt_string(coords[second], str: str) + "\n"
            } else {
                out += str + "\n"
            }
        }
    }
    
    func write_model(model_number: sz, remark: String) -> String {
        var out: String = ""
        out += "MODEL \(model_number)\n"
        out += remark
        
        for i in 0..<ligands.count {
            write_context(ligands.elements[i].cont, &out)
        }
        
        if num_flex > 0 {
            write_context(flex_context, &out)
        }
        
        out += "ENDMDL\n"
        return out
    }
    
    func set(_ c: inout Conf) {
        self.ligands.set_conf(atoms: atoms, coords: &coords, c: c.ligands)
        self.flex   .set_conf(atoms: atoms, coords: &coords, c: c.flex)
    }
    
    func gyration_radius(_ ligand_number: sz) -> fl {
        assert(ligand_number < ligands.count)
        let lig: Ligand = ligands.elements[ligand_number]
        var acc: fl = 0
        var counter: UInt = 0
        for i in lig.begin..<lig.end {
            if atoms[i].el != EL_TYPE_H { // only heavy atoms are used
                acc += vec_distance_sqr(coords[i], lig.node.get_origin())
                counter += 1
            }
        }
        return (counter > 0) ? sqrt(acc / fl(counter)) : 0
    }
    
    func evalo(_ p: precalculate_byatom, _ v: vec) -> fl {
        return eval_interacting_pairs(p, v[2], other_pairs, coords)
    }
    
    func eval_inter(_ p: precalculate_byatom, _ v: vec) -> fl {
        return eval_interacting_pairs(p, v[2], inter_pairs, coords)
    }
    
    func evali(_ p: precalculate_byatom, v: vec) -> fl {
        var e: fl = 0
        for i in 0..<ligands.count {
            e += eval_interacting_pairs(p, v[0], ligands.elements[i].pairs, coords) // probably might was well use coords here
        }
        return e
    }
    
    func eval_deriv(_ p: precalculate_byatom, _ ig: igrid, _ v: vec, _ g: inout Change) -> fl {
        // INTER ligand - grid
        var e: fl = ig.eval_deriv(m: self, v: v[1]) // sets minus_forces, except inflex
        
        // INTRA ligand_i - ligand_i
        for i in 0..<ligands.count {
            e += eval_interacting_pairs_deriv(p, v[0], ligands.elements[i].pairs, coords, &minus_forces) // adds to minus_forces
        }
        
        // INTER ligand_i - ligand_j and ligand_i - flex_i
        if !inter_pairs.isEmpty {
            e += eval_interacting_pairs_deriv(p, v[2], inter_pairs, coords, &minus_forces) // adds to minus_forces
        }
        
        // INTRA flex_i - flex_i and flex_i - flex_j
        if !other_pairs.isEmpty {
            e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, &minus_forces) // adds to minus_forces
        }
        
        // glue_i - glue_i and glue_i - glue_j
        if !glue_pairs.isEmpty {
            e += eval_interacting_pairs_deriv(p, v[2], glue_pairs, coords, &minus_forces, with_max_cutoff: true) // adds to minus_forces
        }
        
        // calculate derivatives
        ligands.derivative(coords: coords, forces: minus_forces, c: &g.ligands)
        flex.derivative(coords: coords, forces: minus_forces, c: &g.flex) // inflex forces are ignored
        
        return e
    }
    
    func eval_intramolecular(_ p: precalculate_byatom, ig: igrid, v: vec) -> fl {
        
        var e: fl = 0
        let cutoff_sqr: fl = p.cutoff_sqr()
        
        // internal for each ligand
        for i in 0..<ligands.count {
            e += eval_interacting_pairs(p, v[0], ligands.elements[i].pairs, coords)
        }
        
        // flex - rigid
        e += ig.eval_intra(m: self, v: v[1])
        
        // flex_i - flex_i and flex_i - flex_j
        for i in 0..<other_pairs.count {
            let pair: InteractingPair = other_pairs[i]
            let r2: fl = vec_distance_sqr(coords[pair.a], coords[pair.b])
            if r2 < cutoff_sqr {
                var this_e: fl = p.eval_fast(i: pair.a, j: pair.b, r2: r2)
                curl(e: &this_e, v: v[2])
                e += this_e
            }
        }
        
        return e
    }
    
    func rmsd_lower_bound_asymmetric(_ x: Model, _ y: Model) -> fl {
        let n: sz = x.m_num_movable_atoms
        assert(n == y.m_num_movable_atoms)
        var sum: fl = 0
        var counter: UInt = 0
        for i in 0..<n {
            let a: Atom = x.atoms[i]
            if a.el != EL_TYPE_H {
                var r2: fl = fl.greatestFiniteMagnitude
                for j in 0..<n {
                    let b: Atom = y.atoms[j]
                    if a.same_element(b) && !b.is_hydrogen() {
                        let this_r2: fl = vec_distance_sqr(x.coords[i], y.coords[j])
                        if this_r2 < r2 {
                            r2 = this_r2
                        }
                    }
                }
                assert(not_max(r2))
                sum += r2
                counter += 1
            }
        }
        return (counter == 0) ? 0 : sqrt(sum / fl(counter))
    }
    
    func rmsd_lower_bound(_ m: Model) -> fl {
        return max(rmsd_lower_bound_asymmetric(self, m), rmsd_lower_bound_asymmetric(m, self))
    }
    
    func rmsd_upper_bound(_ m: Model) -> fl {
        assert(m_num_movable_atoms == m.m_num_movable_atoms)
        var sum: fl = 0
        var counter: UInt = 0
        for i in 0..<m_num_movable_atoms {
            let a: Atom = atoms[i]
            let b: Atom = m.atoms[i]
            assert(a.ad == b.ad)
            assert(a.xs == b.xs)
            if a.el != EL_TYPE_H {
                sum += vec_distance_sqr(coords[i], m.coords[i])
                counter += 1
            }
        }
        return (counter == 0) ? 0 : sqrt(sum / fl(counter))
    }
    
    func rmsd_ligands_upper_bound(_ m: Model) -> fl {
        assert(ligands.count == m.ligands.count)
        var sum: fl = 0
        var counter: UInt = 0
        for ligand_i in 0..<ligands.count {
            let lig: Ligand = ligands.elements[ligand_i]
            let m_lig: Ligand = m.ligands.elements[ligand_i]
            assert(lig.begin == m_lig.begin)
            assert(lig.end == m_lig.end)
            for i in lig.begin..<lig.end {
                let a: Atom = atoms[i]
                let b: Atom = m.atoms[i]
                assert(a.ad == b.ad)
                assert(a.xs == b.xs)
                if a.el != EL_TYPE_H {
                    sum += vec_distance_sqr(coords[i], m.coords[i])
                    counter += 1
                }
            }
        }
        return (counter == 0) ? 0 : sqrt(sum / fl(counter))
    }
    
    func verify_bond_lengths() {
        for i in 0..<grid_atoms.count + atoms.count {
            let ai: AtomIndex = sz_to_atom_index(sz(i))
            let a: Atom = get_atom(ai)
            for j in 0..<a.bonds.count {
                let b: Bond = a.bonds[j]
                let d: fl = sqrt(distance_sqr_between(ai, b.connected_atom_index))
                assert(eq(d, b.length), "d = \(d), b.length = \(b.length)")
            }
        }
    }
    
    
    func check_ligand_internal_pairs() {
        for i in 0..<ligands.count {
            let lig: Ligand = ligands.elements[i]
            for j in lig.pairs {
                let ip: InteractingPair = j
                assert(ip.a >= lig.begin)
                assert(ip.b < lig.end)
            }
        }
    }
    
    public func about() {
        show_variable(name: "atom_typing_used", value: atom_typing_used)
        show_variable(name: "num_movable_atoms", value: num_movable_atoms)
        show_variable(name: "num_internal_pairs", value: num_internal_pairs)
        show_variable(name: "num_other_pairs", value: num_other_pairs)
        show_variable(name: "num_ligands", value: num_ligands)
        show_variable(name: "num_flex", value: num_flex)
    }
    
    public func show_pairs() {
        print("INTER PAIRS")
        for i in inter_pairs {
            let ip: InteractingPair = i
            if is_atom_in_ligand(ip.a) {
                let lig_i: sz = find_ligand(ip.a)
                print("LIGAND (\(lig_i)) : ", separator: "", terminator: "")
            } else {
                print("  FLEX     : ", separator: "", terminator: "")
            }
            if is_atom_in_ligand(ip.b) {
                let lig_i: sz = find_ligand(ip.b)
                print("LIGAND (\(lig_i)) ", separator: "", terminator: "")
            } else {
                print("  FLEX     ", separator: "", terminator: "")
            }
            print(" - ", ip.a, " : ", ip.b, " - ", get_coords(ip.a)[0], " ", get_coords(ip.a)[1], " ", get_coords(ip.a)[2], " - ", get_coords(ip.b)[0], " ", get_coords(ip.b)[1], " ", get_coords(ip.b)[2])
        }
        print("INTRA LIG PAIRS")
        for i in 0..<num_ligands {
            let lig: Ligand = get_ligand(i)
            for j in lig.pairs {
                let ip: InteractingPair = j
                let lig_i: sz = find_ligand(ip.a)
                print("LIGAND (\(lig_i)) ", separator: "", terminator: "")
                print(" - ", ip.a, " : ", ip.b, " - ", get_coords(ip.a)[0], " ", get_coords(ip.a)[1], " ", get_coords(ip.a)[2], " - ", get_coords(ip.b)[0], " ", get_coords(ip.b)[1], " ", get_coords(ip.b)[2])
            }
        }
        print("INTRA FLEX PAIRS")
        for i in 0..<num_other_pairs {
            let ip: InteractingPair = other_pairs[i]
            print("FLEX       ", separator: "", terminator: "")
            print(" - ", ip.a, " : ", ip.b, " - ", get_coords(ip.a)[0], " ", get_coords(ip.a)[1], " ", get_coords(ip.a)[2], " - ", get_coords(ip.b)[0], " ", get_coords(ip.b)[1], " ", get_coords(ip.b)[2])
        }
        print("GLUE - GLUE PAIRS")
        for i in 0..<glue_pairs.count {
            let ip: InteractingPair = glue_pairs[i]
            print("FLEX       ", separator: "", terminator: "")
            print(" - ", ip.a, " : ", ip.b, " - ", get_coords(ip.a)[0], " ", get_coords(ip.a)[1], " ", get_coords(ip.a)[2], " - ", get_coords(ip.b)[0], " ", get_coords(ip.b)[1], " ", get_coords(ip.b)[2])
        }
    }
    
    public func show_atoms() {
        print("ATOM INFORMATION")
        for i in 0..<atoms.count {
            let a: Atom = atoms[i]
            if i < num_movable_atoms {
                print("     MOVABLE: ", separator: "", terminator: "")
            } else {
                print(" NOT MOVABLE: ", separator: "", terminator: "")
            }
            print(i, " - ", coords[i][0], " ", coords[i][1], " ", coords[i][2], " - ", a.ad, " - ", a.xs, " - ", a.charge)
        }
    }
    
    func show_forces() {
        print("ATOM FORCES")
        for i in 0..<atoms.count {
            print(i, " ", minus_forces[i][0], " ", minus_forces[i][1], " ", minus_forces[i][2])
        }
    }
    
    public func print_stuff(_ show_coords: Bool = true,
                            _ show_internal: Bool = true,
                            _ show_atoms: Bool = true,
                            _ show_grid: Bool = true,
                            _ show_about: Bool = true) {
        if show_coords {
            print("coords:\n")
            for i in 0..<coords.count {
                print(coords[i], separator: " ")
            }
        }
        if show_atoms {
            print("atoms:\n")
            for i in 0..<atoms.count {
                let a: Atom = atoms[i]
                print(a.el, a.ad, a.xs, a.sy, a.charge, separator: "\t")
                print(a.bonds.count, a.coords, separator: " ")
            }
        }
        if show_grid {
            print("grid_atoms:\n")
            for i in 0..<grid_atoms.count {
                let a: Atom = grid_atoms[i]
                print(a.el, a.ad, a.xs, a.sy, a.charge, separator: "\t")
                print(a.bonds.count, a.coords, separator: " ")
            }
        }
        if show_about {
            about()
        }
    }
    
    func clash_penalty_aux(_ pairs: InteractingPairs) -> fl {
        var e: fl = 0
        for i in pairs {
            let r: fl = sqrt(vec_distance_sqr(coords[i.a], coords[i.b]))
            let covalent_r: fl = atoms[i.a].covalent_radius() + atoms[i.b].covalent_radius()
            e += pairwise_clash_penalty(r: r, covalent_r: covalent_r)
        }
        return e
    }
    
    func clash_penalty() -> fl {
        var e: fl = 0
        for i in 0..<ligands.count {
            e += clash_penalty_aux(ligands.elements[i].pairs)
        }
        e += clash_penalty_aux(other_pairs)
        return e
    }
        
    func write_structure(_ out: inout TextStreamer) {
        for i in 0..<ligands.count {
            write_context(ligands.elements[i].cont, &out)
        }
        if num_flex > 0 {
            write_context(flex_context, &out)
        }
    }
    
    func write_structure(_ out: inout TextStreamer, _ remark: String) {
        print(remark, to: &out)
        write_structure(&out)
    }
    
    func write_structure(_ out: inout TextStreamer, remarks: [String]) { 
        for i in remarks {
            print(i, to: &out)
        }
        write_structure(&out)
    }
    
    func write_structure(_ name: String) {
        var out: TextStreamer = TextStreamer(filepath: name)
        write_structure(&out)
    }
    
    func write_model(_ out: inout TextStreamer, model_number: sz, remark: String) {
        print("MODEL \(model_number)", to: &out)
        write_structure(&out, remark)
        print("ENDMDL", to: &out)
    }
    
}
