//
//  PDBQTParser.swift
//  
//
//  Created by Cory Kornowicz on 11/14/23.
//

import Foundation
import simd

// can throw pdbqt_parse_error
func parse_receptor_pdbqt(rigid: String, flex: String? = nil, atype: AtomType.T = .XS) throws -> Model {
    // Parse PDBQT receptor with flex residues

    var r: Rigid = Rigid(atoms: [])
    var nrp: NonRigidParsed = NonRigidParsed(ligands: VectorMutable<Ligand, LigandConf, LigandChange>(), 
                                             flex: VectorMutable<Residue, ResidueConf, ResidueChange>(),
                                             atoms: [], inflex: [], atoms_atoms_bonds: DistanceTypeMatrix(),
                                             atoms_inflex_bonds: matrix<DistanceType>(), inflex_inflex_bonds: DistanceTypeMatrix())
    var c: Context = Context()

    let tmp: PDBQTInitializer = PDBQTInitializer(atom_typing_used: atype)

    if !rigid.isEmpty {
        let rigidParser: RowIterator = RowIterator(filePath: rigid)
        try parse_pdbqt_rigid(in_: rigidParser, r: &r)
    }

    if let flex = flex, !flex.isEmpty {
        let flexParser: RowIterator = RowIterator(filePath: flex)
        try parse_pdbqt_flex(in_: flexParser, nr: &nrp, c: &c)
    }
    
    if !rigid.isEmpty {
        try tmp.initialize_from_rigid(r)
        if flex == nil {
            let mobility_matrix: DistanceTypeMatrix = DistanceTypeMatrix()
            try tmp.initialize(mobility: mobility_matrix)
        }
    }

    if let flex = flex, !flex.isEmpty {
        tmp.initialize_from_nrp(nrp: nrp, c: c, is_ligand: false)
        try tmp.initialize(mobility: nrp.mobility_matrix())
    }

    return tmp.m
}

// can throw pdbqt_parse_error
func parse_ligand_pdbqt_from_file(filePath: String, atype: AtomType.T) throws -> Model {
    
    var nrp: NonRigidParsed = NonRigidParsed(ligands: VectorMutable<Ligand, LigandConf, LigandChange>(), flex: VectorMutable<Residue, ResidueConf, ResidueChange>(), 
                                             atoms: [], inflex: [], atoms_atoms_bonds: DistanceTypeMatrix(),
                                             atoms_inflex_bonds: matrix<DistanceType>(), inflex_inflex_bonds: DistanceTypeMatrix())
    var c: Context = Context()
    
    let rowIterator: RowIterator = RowIterator(filePath: filePath)

    try parse_pdbqt_ligand(rowIterator: rowIterator, nr: &nrp, c: &c)

    let tmp: PDBQTInitializer = PDBQTInitializer(atom_typing_used: atype)
    
    tmp.initialize_from_nrp(nrp: nrp, c: c, is_ligand: true)
    
    try tmp.initialize(mobility: nrp.mobility_matrix())
    
    return tmp.m
}

func parse_ligand_pdbqt_from_string(stringContents: String, atype: AtomType.T) throws -> Model {
    // can throw parse_error
    var nrp: NonRigidParsed = NonRigidParsed(ligands: VectorMutable<Ligand, LigandConf, LigandChange>(), flex: VectorMutable<Residue, ResidueConf, ResidueChange>(), atoms: [], inflex: [], atoms_atoms_bonds: DistanceTypeMatrix(), atoms_inflex_bonds: matrix<DistanceType>(), inflex_inflex_bonds: DistanceTypeMatrix())
    var c: Context = Context()
    
    let rowIterator: RowIterator = RowIterator(stringContents: stringContents)

    try parse_pdbqt_ligand(in_: rowIterator, nr: &nrp, c: &c)

    let tmp: PDBQTInitializer = PDBQTInitializer(atom_typing_used: atype)
    tmp.initialize_from_nrp(nrp: nrp, c: c, is_ligand: true)
    try tmp.initialize(mobility: nrp.mobility_matrix())
    return tmp.m
}

final class ParsedAtom: Atom {
    var number: sz
    init(ad_: sz, charge_: fl, coords_: vec, number_: sz) {
        self.number = number_
        super.init(coords: coords_, bonds: [])
        self.ad = ad_
        self.charge = charge_
    }
}

func add_context(c: inout Context, str: String) {
    c.append(ParsedLine(str, nil)) // in original it was a boost::optional<sz> value?
}

func checked_convert_substring<T>(_ str: String, i: Int, j: Int, dest_nature: String, converter: @escaping (String) -> T) throws -> T {
    
    assert(i >= 1)
    assert(i <= j+1)
    
    if j > str.count {
        throw VinaError("This line is too short: " + str)
    }
    
    // cast range to substring and then remove all whitespace
    let tmp = str[slow: i..<j]
    let substring = tmp.trimmingCharacters(in: .whitespaces)
    
    // convert string to T type
    return converter(substring)
}

func parse_pdbqt_atom_string(str: String) throws -> ParsedAtom {
    
    let number: sz = try checked_convert_substring(str, i: 7, j: 11, dest_nature: "Atom number") { str in
        return sz((str as NSString).integerValue)
    }
    
    var coords: vec = zero_vec
    coords[0] = try checked_convert_substring(str, i: 31, j: 38, dest_nature: "Coordinate") { str in
        return fl(str)!
    }
    coords[1] = try checked_convert_substring(str, i: 39, j: 46, dest_nature: "Coordinate") { str in
        return fl(str)!
    }
    coords[2] = try checked_convert_substring(str, i: 47, j: 54, dest_nature: "Coordinate") { str in
        return fl(str)!
    }

    var charge: fl = 0
    if !str[slow: 69..<76].isEmpty {
        charge = try checked_convert_substring(str, i: 69, j: 76, dest_nature: "Charge") { str in
            return fl(str)!
        }
    }

    var name: String = ""
    if !str[slow: 77...78].isEmpty {
        name = str[slow: 77...78].trimmingCharacters(in: .whitespaces)
    }

    let ad: sz = string_to_ad_type(name: name)
    
    let tmp = ParsedAtom(ad_: ad, charge_: charge, coords_: coords, number_: number)

    if is_non_ad_metal_name(name: name) {
        tmp.xs = XS_TYPE_Met_D
    }

    if tmp.acceptable_type() {
        return tmp
    } else {
        throw VinaError("Atom type " + name + " is not a valid AutoDock type (atom types are case-sensitive).")
    }
}

class AtomReference : NSCopying {
    var index: sz
    var inflex: Bool
    
    init(index: sz, inflex: Bool) {
        self.index = index
        self.inflex = inflex
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        return AtomReference(index: self.index, inflex: self.inflex)
    }
}

final class MoveableAtom: Atom {
    
    var relative_coords: vec
    
    init(_ a: Atom, relative_coords: vec) {
        self.relative_coords = relative_coords
        super.init(coords: a.coords, bonds: a.bonds)
        self.charge = a.charge
        self.ad = a.ad
        self.el = a.el
        self.sy = a.sy
        self.xs = a.xs
    }
    
    private init() {
        self.relative_coords = zero_vec
        super.init(coords: zero_vec, bonds: [])
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newAtom = MoveableAtom()
        let newBonds = ContiguousArray<Bond>(self.bonds.map { $0.copy() as! Bond })
        newAtom.relative_coords = self.relative_coords
        newAtom.coords = self.coords
        newAtom.bonds = newBonds
        newAtom.charge = self.charge
        newAtom.xs = self.xs
        newAtom.el = self.el
        newAtom.ad = self.ad
        newAtom.sy = self.sy
        return newAtom
    }
    
}

final class Rigid: NSCopying {
    var atoms: atomv
    
    init(atoms: atomv) {
        self.atoms = atoms
    }
    
    func copy(with zone: NSZone? = nil) -> Any {
        let copyAtoms = ContiguousArray(self.atoms.map({ $0.copy() as! Atom }))
        let newRigid = Rigid(atoms: copyAtoms)
        return newRigid
    }
}

typealias mav = ContiguousArray<MoveableAtom>

class NonRigidParsed {
    let ligands: VectorMutable<Ligand, LigandConf, LigandChange>
    let flex: VectorMutable<Residue, ResidueConf, ResidueChange>
    
    var atoms: mav = []
    var inflex: atomv = []
    
    var atoms_atoms_bonds: DistanceTypeMatrix
    var atoms_inflex_bonds: matrix<DistanceType>
    var inflex_inflex_bonds: DistanceTypeMatrix
    
    init(ligands: VectorMutable<Ligand, LigandConf, LigandChange>, flex: VectorMutable<Residue, ResidueConf, ResidueChange>, 
         atoms: mav, inflex: atomv,
         atoms_atoms_bonds: DistanceTypeMatrix, atoms_inflex_bonds: matrix<DistanceType>, inflex_inflex_bonds: DistanceTypeMatrix) {
        self.ligands = ligands
        self.flex = flex
        self.atoms = atoms
        self.inflex = inflex
        self.atoms_atoms_bonds = atoms_atoms_bonds
        self.atoms_inflex_bonds = atoms_inflex_bonds
        self.inflex_inflex_bonds = inflex_inflex_bonds
    }
    
    func mobility_matrix() -> DistanceTypeMatrix {
        let tmp: DistanceTypeMatrix = atoms_atoms_bonds
        tmp.append(rectangular: atoms_inflex_bonds, triangular: inflex_inflex_bonds)
        return tmp
    }
}


protocol ParsingStructProtocol { 
    var axis_begin: AtomReference? { get set }
    var axis_end: AtomReference? { get set }

    func insert_immobile_inflex(nr: inout NonRigidParsed)
    func insert_immobile(nr: inout NonRigidParsed, c: inout Context, frame_origin: vec)
}

final class ParsingStruct: ParsingStructProtocol {
    
    final class node_t<T: ParsingStructProtocol> {
        var context_index: sz
        var a: ParsedAtom
        var ps: ContiguousArray<T> = ContiguousArray<T>()
        
        init(context_index: sz, a: ParsedAtom) {
            self.context_index = context_index
            self.a = a
        }

        // inflex atom insertion
        func insert_inflex(nr: inout NonRigidParsed) {
            for i in 0..<ps.count {
                ps[i].axis_begin = AtomReference(index: sz(nr.inflex.count), inflex: true)
            }
            nr.inflex.append(a)
        }

        func insert_immobiles_inflex(nr: inout NonRigidParsed) {
            for i in 0..<ps.count {
                ps[i].insert_immobile_inflex(nr: &nr)
            }
        }
        // insertion into non_rigid_parsed
        func insert(nr: inout NonRigidParsed, c: inout Context, frame_origin: vec) {
            for i in 0..<ps.count {
                ps[i].axis_begin = AtomReference(index: sz(nr.atoms.count), inflex: false)
            }
            let relative_coords: vec = a.coords - frame_origin
            c[context_index].1 = sz(nr.atoms.count)
            nr.atoms.append(MoveableAtom(a, relative_coords: relative_coords))
        }
        func insert_immobiles(nr: inout NonRigidParsed, c: inout Context, frame_origin: vec) {
            for i in 0..<ps.count {
                ps[i].insert_immobile(nr: &nr, c: &c, frame_origin: frame_origin)
            }
        }

    }

    typealias node = ParsingStruct.node_t<ParsingStruct>
    
    var immobile_atom: sz? = nil // which of `atoms' is immobile, if any
    var axis_begin: AtomReference? = nil // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
    var axis_end: AtomReference? = nil // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
    var atoms: ContiguousArray<node> = []

    func add(a: ParsedAtom, c: Context) {
        assert(c.count > 0)
        atoms.append(node(context_index: sz(c.count - 1), a: a))
    }
    
    func countTotalAtoms() -> Int {
        var sum = 0
        for node in atoms {
            sum += 1
            for sub_node in node.ps {
                sum += sub_node.countTotalAtoms()
            }
        }
        return sum
    }

    func immobile_atom_coords() -> vec {
        assert(immobile_atom != nil)
        assert(immobile_atom! < atoms.count)
        return atoms[immobile_atom!].a.coords
    }
    // inflex insertion
    func insert_immobile_inflex(nr: inout NonRigidParsed) {
        if(!atoms.isEmpty) {
            assert(immobile_atom != nil)
            assert(immobile_atom! < atoms.count)
            axis_end = AtomReference(index: sz(nr.inflex.count), inflex: true)
            atoms[immobile_atom!].insert_inflex(nr: &nr)
        }
    }

    // insertion into non_rigid_parsed
    func insert_immobile(nr: inout NonRigidParsed, c: inout Context, frame_origin: vec) {
        if(!atoms.isEmpty) {
            assert(immobile_atom != nil)
            assert(immobile_atom! < atoms.count)
            axis_end = AtomReference(index: sz(nr.atoms.count), inflex: false)
            atoms[immobile_atom!].insert(nr: &nr, c: &c, frame_origin: frame_origin)
        }
    }
    
    func essentially_empty() -> Bool { // no sub-branches besides immobile atom, including sub-sub-branches, etc
        for i in 0..<atoms.count {
            if let imAtom = immobile_atom, imAtom != i { return false }
            if(!atoms[i].ps.isEmpty) {
                return false // FIXME : iffy
            }
        }
        return true
    }
}

func parse_one_unsigned(str: String, start: String) throws -> UInt {    
    let substr = str[slow: start.count..<str.count]
    let components = substr.components(separatedBy: " ").filter({ !$0.isEmpty })
    let first = UInt((components[0] as NSString).integerValue)
    return first
}

func parse_two_unsigneds(str: String, start: String) throws -> (UInt, UInt) {
    let substr = str[slow: start.count..<str.count]
    let components = substr.components(separatedBy: " ").filter({ !$0.isEmpty })
    let first = UInt((components[0] as NSString).integerValue)
    let second = UInt((components[1] as NSString).integerValue)
    return (first, second)
}

func parse_pdbqt_rigid(in_: some RowIterable<String>, r: inout Rigid) throws {

    while let line = in_.next() {
        if line.isEmpty { return }
        if line.hasPrefix("TER") { return } // ignore
        if line.hasPrefix("END") { return } // ignore 
        if line.hasPrefix("WARNING") { return } // ignore - AutoDockTools bug workaround
        if line.hasPrefix("REMARK") { return } // ignore
        if line.hasPrefix("ATOM  ") || line.hasPrefix("HETATM") {
            let tempPAtom = try parse_pdbqt_atom_string(str: line)
            r.atoms.append(tempPAtom)
        }
        else if line.hasPrefix("MODEL") {
            throw VinaError("Unexpected multi-MODEL tag found in rigid receptor.\nOnly one model can be used for the rigid receptor.")
        }
        else {
            throw VinaError("Unknown or inappropriate tag found in rigid receptor: \(line).")
        }
    }

}

//void parse_pdbqt_root_aux(std::istream& in, parsing_struct& p, context& c) {
//    std::string str;
//
//    while(std::getline(in, str)) {
//        add_context(c, str);
//
//        if(str.empty()) {} // ignore ""
//        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
//        else if(starts_with(str, "REMARK")) {} // ignore
//        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM"))
//                p.add(parse_pdbqt_atom_string(str), c);
//        else if(starts_with(str, "ENDROOT"))
//                return;
//        else if(starts_with(str, "MODEL"))
//                throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
//                                        "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
//        else
//            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
//    }
//}

func parse_pdbqt_root_aux(in_: some RowIterable<String>, p: inout ParsingStruct, c: inout Context) throws {
    while let line = in_.next() {
        if line.isEmpty { return } // ignore ""
        add_context(c: &c, str: line)
        if line.hasPrefix("WARNING") { return } // ignore - AutoDockTools bug workaround
        else if line.hasPrefix("REMARK") { return } // ignore
        else if line.hasPrefix("ATOM") || line.hasPrefix("HETATM") {
            p.add(a: try! parse_pdbqt_atom_string(str: line), c: c)
        }
        else if line.hasPrefix("ENDROOT") {
            return
        }
        else if line.hasPrefix("MODEL") {
            throw VinaError("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.")
        }
        else {
            throw VinaError("Unknown or inappropriate tag found in flex residue or ligand: \(line)")
        }
    }
}

//void parse_pdbqt_root(std::istream& in, parsing_struct& p, context& c) {
//    std::string str;
//
//    while(std::getline(in, str)) {
//        add_context(c, str);
//
//        if(str.empty()) {} // ignore
//        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
//        else if(starts_with(str, "REMARK")) {} // ignore
//        else if(starts_with(str, "ROOT")) {
//            parse_pdbqt_root_aux(in, p, c);
//            break;
//        }
//        else if(starts_with(str, "MODEL"))
//                throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
//                                        "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
//        else
//            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
//    }
//}

func parse_pdbqt_root(in_: some RowIterable<String>, p: inout ParsingStruct, c: inout Context) throws {
    while let line = in_.next() {
        if line.isEmpty { continue } // ignore
        add_context(c: &c, str: line)
        if line.hasPrefix("WARNING") { } // ignore - AutoDockTools bug workaround
        else if line.hasPrefix("REMARK") { } // ignore
        else if line.hasPrefix("ROOT") {
            try parse_pdbqt_root_aux(in_: in_, p: &p, c: &c)
            return
        }
        else if line.hasPrefix("MODEL") {
            throw VinaError("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.")
        }
        else {
            throw VinaError("Unknown or inappropriate tag found in flex residue or ligand: \(line).")
        }
    }
}

//void parse_pdbqt_branch_aux(std::istream& in, const std::string& str, parsing_struct& p, context& c) {
//    unsigned first, second;
//    parse_two_unsigneds(str, "BRANCH", first, second);
//    sz i = 0;
//
//    for(; i < p.atoms.size(); ++i)
//        if(p.atoms[i].a.number == first) {
//        p.atoms[i].ps.push_back(parsing_struct());
//        parse_pdbqt_branch(in, p.atoms[i].ps.back(), c, first, second);
//        break;
//    }
//
//    if(i == p.atoms.size())
//        throw pdbqt_parse_error("Atom number " + std::to_string(first) + " is missing in this branch.", str);
//}

func parse_pdbqt_branch_aux(in_: some RowIterable<String>, str: String, p: inout ParsingStruct, c: inout Context) throws {
    let components = str.components(separatedBy: " ").filter({ !$0.isEmpty })
    let first = UInt((components[1] as NSString).integerValue)
    let second = UInt((components[2] as NSString).integerValue)
    var i: sz = 0
    while i < p.atoms.count {
        if p.atoms[i].a.number == first {
            p.atoms[i].ps.append(ParsingStruct())
            try p.atoms[i].ps.withLastElement { new_p in
                try parse_pdbqt_branch(in_: in_, p: &new_p, c: &c, from: first, to: second)
            }
            break
        }
        i+=1
    }
    
    if i == p.atoms.count {
        throw VinaError("Atom number \(first) is missing in this branch.")
    }
}

//void parse_pdbqt_aux(std::istream& in, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue) {
//    parse_pdbqt_root(in, p, c);
//
//    std::string str;
//
//    while(std::getline(in, str)) {
//        add_context(c, str);
//
//        if(str.empty()) {} // ignore ""
//        if(str[0] == '\0') {} // ignore a different kind of emptiness (potential issues on Windows)
//        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
//        else if(starts_with(str, "REMARK")) {} // ignore
//        else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, str, p, c);
//        else if(!residue && starts_with(str, "TORSDOF")) {
//            if(torsdof)
//                throw pdbqt_parse_error("TORSDOF keyword can be defined only once.");
//            torsdof = parse_one_unsigned(str, "TORSDOF");
//        }
//        else if(residue && starts_with(str, "END_RES"))
//                return;
//        else if(starts_with(str, "MODEL"))
//                throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
//                                        "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
//        else
//            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
//    }
//}

func parse_pdbqt_aux(in_: some RowIterable<String>, p: inout ParsingStruct, c: inout Context, torsdof: inout UInt?, residue: Bool) throws {
    try parse_pdbqt_root(in_: in_, p: &p, c: &c)

    while let line = in_.next() {
        if line.isEmpty { continue } // ignore ""
        add_context(c: &c, str: line)
        if line.hasPrefix("WARNING") { } // ignore - AutoDockTools bug workaround
        else if line.hasPrefix("REMARK") { } // ignore
        else if line.hasPrefix("BRANCH") { try parse_pdbqt_branch_aux(in_: in_, str: line, p: &p, c: &c) }
        else if !residue && line.hasPrefix("TORSDOF") {
            if torsdof != nil {
                throw VinaError("TORSDOF keyword can be defined only once.")
            }
            torsdof = try parse_one_unsigned(str: line, start: "TORSDOF")
        }
        else if residue && line.hasPrefix("END_RES") {
            return
        }
        else if line.hasPrefix("MODEL") {
            throw VinaError("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.")
        }
        else {
            throw VinaError("Unknown or inappropriate tag found in flex residue or ligand.")
        }
    }
}

//void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
//    if(atm)
//        VINA_RANGE(i, r.begin, r.end) {
//        atom_reference& ar = atm.get();
//        if(ar.inflex)
//            nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED; //(max_unsigned); // first index - atoms, second index - inflex
//        else
//            nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED; // (max_unsigned);
//    }
//}

func add_bonds(nr: inout NonRigidParsed, atm: AtomReference?, r: AtomRange) {
    if let ar = atm {
        for i in r.begin..<r.end {
            if ar.inflex {
                nr.atoms_inflex_bonds[i, ar.index] = .DISTANCE_FIXED // first index - atoms, second index - inflex
            } else {
                nr.atoms_atoms_bonds[ar.index, i] = .DISTANCE_FIXED
            }
        }
    }
}

//void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
//    if(axis_begin && axis_end) {
//        atom_reference& r1 = axis_begin.get();
//        atom_reference& r2 = axis_end  .get();
//        if(r2.inflex) {
//            VINA_CHECK(r1.inflex); // no atom-inflex rotors
//            nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
//        }
//        else
//            if(r1.inflex)
//                nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (atoms, inflex)
//        else
//            nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
//    }
//}

func set_rotor(nr: inout NonRigidParsed, axis_begin: AtomReference?, axis_end: AtomReference?) {
    if let r1 = axis_begin, let r2 = axis_end {
        if r2.inflex {
            assert(r1.inflex) // no atom-inflex rotors
            nr.inflex_inflex_bonds[r1.index, r2.index] = .DISTANCE_ROTOR
        } else if r1.inflex {
            nr.atoms_inflex_bonds[r2.index, r1.index] = .DISTANCE_ROTOR // (atoms, inflex)
        } else {
            nr.atoms_atoms_bonds[r1.index, r2.index] = .DISTANCE_ROTOR
        }
    }
}

typealias axis_numbers = Pair<sz, sz>
typealias axis_numbers_option = axis_numbers?

//void nr_update_matrixes(non_rigid_parsed& nr) {
//    // atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)
//
//    nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);
//    nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // first index - inflex, second index - atoms
//    nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
//}

func nr_update_matrixes(nr: inout NonRigidParsed) {
    // atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)
    nr.atoms_atoms_bonds.resize(n: sz(nr.atoms.count), filler_val: .DISTANCE_VARIABLE)
    nr.atoms_inflex_bonds.resize(m: sz(nr.atoms.count), n: sz(nr.inflex.count), filler_val: .DISTANCE_VARIABLE) // first index - inflex, second index - atoms
    nr.inflex_inflex_bonds.resize(n: sz(nr.inflex.count), filler_val: .DISTANCE_FIXED) // FIXME?
}

//template<typename B> // B == branch or main_branch or flexible_body
//void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
//    b.node.begin = nr.atoms.size();
//    VINA_FOR_IN(i, p.atoms) {  // postprocess atoms into 'b.node'
//        parsing_struct::node& p_node = p.atoms[i];
//        if(p.immobile_atom && i == p.immobile_atom.get()) {} // skip immobile_atom - it's already inserted in "THERE"
//        else p_node.insert(nr, c, b.node.get_origin());
//        p_node.insert_immobiles(nr, c, b.node.get_origin());
//    }
//    b.node.end = nr.atoms.size();
//
//    nr_update_matrixes(nr);
//    add_bonds(nr, p.axis_begin, b.node); // b.node is used as atom_range
//    add_bonds(nr, p.axis_end  , b.node); // b.node is used as atom_range
//    set_rotor(nr, p.axis_begin, p.axis_end);
//
//    VINA_RANGE(i, b.node.begin, b.node.end)
//    VINA_RANGE(j, i+1, b.node.end)
//    nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME
//
//
//    VINA_FOR_IN(i, p.atoms) {   // postprocess children
//        parsing_struct::node& p_node = p.atoms[i];
//        VINA_FOR_IN(j, p_node.ps) {
//            parsing_struct& ps = p_node.ps[j];
//            if(!ps.essentially_empty()) { // immobile already inserted // FIXME ?!
//                b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch will assign begin and end
//                postprocess_branch(nr, ps, c, b.children.back());
//            }
//        }
//    }
//    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
//    VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
//    VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
//}

func postprocess_branch<B: NodeBearing>(nr: inout NonRigidParsed, p: inout ParsingStruct, c: inout Context, b: inout B) {
    
    b.node.begin = sz(nr.atoms.count)
    
    for i in 0..<p.atoms.count { // postprocess atoms into 'b.node'
        if let pAtom = p.immobile_atom, pAtom == i { } // skip immobile_atom - it's already inserted in "THERE"
        else { p.atoms[i].insert(nr: &nr, c: &c, frame_origin: b.node.get_origin()) }
        
        p.atoms[i].insert_immobiles(nr: &nr, c: &c, frame_origin: b.node.get_origin())
    }
    
    b.node.end = sz(nr.atoms.count)
    
    nr_update_matrixes(nr: &nr)
    
    add_bonds(nr: &nr, atm: p.axis_begin, r: b.node) // b.node is used as atom_range
    add_bonds(nr: &nr, atm: p.axis_end, r: b.node) // b.node is used as atom_range
    set_rotor(nr: &nr, axis_begin: p.axis_begin, axis_end: p.axis_end)

    for i in b.node.begin..<b.node.end {
        for j in i+1..<b.node.end {
            nr.atoms_atoms_bonds[i, j] = .DISTANCE_FIXED // FIXME
        }
    }

    for i in 0..<p.atoms.count { // postprocess children
        for j in 0..<p.atoms[i].ps.count {
            if !p.atoms[i].ps[j].essentially_empty() { // immobile already inserted // FIXME ?!
                // postprocess_branch will assign begin and end
                var new_segment = Branch(node: Segment(origin_: p.atoms[i].ps[j].immobile_atom_coords(), begin_: 0, end_: 0, axis_root: p.atoms[i].a.coords, parent: b.node), children: [])
                
                postprocess_branch(nr: &nr, p: &p.atoms[i].ps[j], c: &c, b: &new_segment)
                
                b.children.append(new_segment)
            }
        }
    }

    assert(nr.atoms_atoms_bonds.dim() == nr.atoms.count)
    assert(nr.atoms_inflex_bonds.dim_1() == nr.atoms.count)
    assert(nr.atoms_inflex_bonds.dim_2() == nr.inflex.count)
}

//void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
//    VINA_CHECK(!p.atoms.empty());
//    nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch will assign begin and end
//    postprocess_branch(nr, p, c, nr.ligands.back());
//    nr_update_matrixes(nr); // FIXME ?
//}

func postprocess_ligand(nr: inout NonRigidParsed, p: inout ParsingStruct, c: inout Context, torsdof: UInt) {
    assert(!p.atoms.isEmpty)
    var new_ligand = Ligand(f: FlexibleBody(node: RigidBody(origin_: p.atoms[0].a.coords, begin_: 0, end_: 0), children: []), degrees_of_freedom: torsdof)
    // postprocess_branch will assign begin and end
    postprocess_branch(nr: &nr, p: &p, c: &c, b: &new_ligand)
    nr.ligands.elements.append(new_ligand)
    nr_update_matrixes(nr: &nr) // FIXME ?
}

//void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
//    VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
//        parsing_struct::node& p_node = p.atoms[i];
//        p_node.insert_inflex(nr);
//        p_node.insert_immobiles_inflex(nr);
//    }
//    VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
//        parsing_struct::node& p_node = p.atoms[i];
//        VINA_FOR_IN(j, p_node.ps) {
//            parsing_struct& ps = p_node.ps[j];
//            if(!ps.essentially_empty()) { // immobile atom already inserted // FIXME ?!
//                nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords))); // postprocess_branch will assign begin and end
//                postprocess_branch(nr, ps, c, nr.flex.back());
//            }
//        }
//    }
//    nr_update_matrixes(nr); // FIXME ?
//    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
//    VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
//    VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
//}

func postprocess_residue(nr: inout NonRigidParsed, p: inout ParsingStruct, c: inout Context) {
    for i in 0..<p.atoms.count { // iterate over "root" of a "residue"
        p.atoms[i].insert_inflex(nr: &nr)
        p.atoms[i].insert_immobiles_inflex(nr: &nr)
    }
    for i in 0..<p.atoms.count { // iterate over "root" of a "residue"
        for j in 0..<p.atoms[i].ps.count {
            if !p.atoms[i].ps[j].essentially_empty() { // immobile atom already inserted // FIXME ?!
                var new_main_branch = MainBranch(node: FirstSegment(origin_: p.atoms[i].ps[j].immobile_atom_coords(), begin_: 0, end_: 0, axis_root: p.atoms[i].a.coords), children: [])
                // postprocess_branch will assign begin and end
                postprocess_branch(nr: &nr, p: &p.atoms[i].ps[j], c: &c, b: &new_main_branch)
                nr.flex.elements.append(new_main_branch as! Residue)
            }
        }
    }
    nr_update_matrixes(nr: &nr) // FIXME ?
    assert(nr.atoms_atoms_bonds.dim() == nr.atoms.count)
    assert(nr.atoms_inflex_bonds.dim_1() == nr.atoms.count)
    assert(nr.atoms_inflex_bonds.dim_2() == nr.inflex.count)
} 

//dkoes, stream version
//void parse_pdbqt_ligand(std::istream& in, non_rigid_parsed& nr, context& c) {
//    parsing_struct p;
//    boost::optional<unsigned> torsdof;
//
//    parse_pdbqt_aux(in, p, c, torsdof, false);
//
//    if(p.atoms.empty())
//        throw pdbqt_parse_error("No atoms in this ligand.");
//    if(!torsdof)
//        throw pdbqt_parse_error("Missing TORSDOF keyword.");
//
//    postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
//
//    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
//}

func parse_pdbqt_ligand(in_: some RowIterable<String>, nr: inout NonRigidParsed, c: inout Context) throws {
    var p: ParsingStruct = ParsingStruct()
    var torsdof: UInt? = nil
    
    try parse_pdbqt_aux(in_: in_, p: &p, c: &c, torsdof: &torsdof, residue: false)
    
    if p.atoms.isEmpty {
        throw VinaError("No atoms in this ligand.")
    }
    if torsdof == nil {
        throw VinaError("Missing TORSDOF keyword.")
    }
    
    postprocess_ligand(nr: &nr, p: &p, c: &c, torsdof: torsdof!) 
    
    assert(nr.atoms_atoms_bonds.dim() == nr.atoms.count)
}

//void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
//    ifile in(name);
//    parsing_struct p;
//    boost::optional<unsigned> torsdof;
//
//    parse_pdbqt_aux(in, p, c, torsdof, false);
//
//    if(p.atoms.empty())
//        throw pdbqt_parse_error("No atoms in this ligand.");
//    if(!torsdof)
//        throw pdbqt_parse_error("Missing TORSDOF keyword in this ligand.");
//
//    postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
//
//    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
//}

func parse_pdbqt_ligand(rowIterator: some RowIterable<String>, nr: inout NonRigidParsed, c: inout Context) throws {
    
    var p: ParsingStruct = ParsingStruct()
    var torsdof: UInt? = nil

    try parse_pdbqt_aux(in_: rowIterator, p: &p, c: &c, torsdof: &torsdof, residue: false)

    if p.atoms.isEmpty {
        throw VinaError("No atoms in this ligand.")
    }

    if torsdof == nil {
        throw VinaError("Missing TORSDOF keyword in this ligand.")
    }

    postprocess_ligand(nr: &nr, p: &p, c: &c, torsdof: torsdof!)

    assert(nr.atoms_atoms_bonds.dim() == nr.atoms.count)
}

//void parse_pdbqt_residue(std::istream& in, parsing_struct& p, context& c) {
//    boost::optional<unsigned> dummy;
//    parse_pdbqt_aux(in, p, c, dummy, true);
//}

func parse_pdbqt_residue(in_: some RowIterable<String>, p: inout ParsingStruct, c: inout Context) throws {
    var dummy: UInt? = nil
    try parse_pdbqt_aux(in_: in_, p: &p, c: &c, torsdof: &dummy, residue: true)
}

//void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
//    ifile in(name);
//    std::string str;
//    
//    while(std::getline(in, str)) {
//        add_context(c, str);
//        
//        if(str.empty()) {} // ignore ""
//        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
//        else if(starts_with(str, "REMARK")) {} // ignore
//        else if(starts_with(str, "BEGIN_RES")) {
//            parsing_struct p;
//            parse_pdbqt_residue(in, p, c);
//            postprocess_residue(nr, p, c);
//        }
//        else if(starts_with(str, "MODEL"))
//                throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue PDBQT file. "
//                                        "Use \"vina_split\" to split flex residues in multiple PDBQT files.");
//        else
//            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue.", str);
//    }
//    
//    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
//}

func parse_pdbqt_flex(in_: some RowIterable<String>, nr: inout NonRigidParsed, c: inout Context) throws {
    while let line = in_.next() {
        if line.isEmpty { continue } // ignore ""
        add_context(c: &c, str: line)
        if line.hasPrefix("WARNING") { } // ignore - AutoDockTools bug workaround
        else if line.hasPrefix("REMARK") { } // ignore
        else if line.hasPrefix("BEGIN_RES") {
            var p: ParsingStruct = ParsingStruct()
            try parse_pdbqt_residue(in_: in_, p: &p, c: &c)
            postprocess_residue(nr: &nr, p: &p, c: &c)
        }
        else if line.hasPrefix("MODEL") {
            throw VinaError("Unexpected multi-MODEL tag found in flex residue PDBQT file. Use \"vina_split\" to split flex residues in multiple PDBQT files.")
        }
        else {
            throw VinaError("Unknown or inappropriate tag found in flex residue.")
        }
    }

    assert(nr.atoms_atoms_bonds.dim() == nr.atoms.count)
}
 
func parse_pdbqt_branch(in_: some RowIterable<String>, p: inout ParsingStruct, c: inout Context, from: UInt, to: UInt) throws {
    while let line = in_.next() {
        if line.isEmpty { continue } //ignore ""
        add_context(c: &c, str: line)
        if line.hasPrefix("WARNING") { } // ignore - AutoDockTools bug workaround
        else if line.hasPrefix("REMARK") { } // ignore
        else if line.hasPrefix("BRANCH") { try parse_pdbqt_branch_aux(in_: in_, str: line, p: &p, c: &c) }
        else if line.hasPrefix("ENDBRANCH") {
            let components = line.components(separatedBy: " ").filter({ !$0.isEmpty })
            let first = UInt((components[1] as NSString).integerValue)
            let second = UInt((components[2] as NSString).integerValue)
            if first != from || second != to {
                throw VinaError("Inconsistent branch numbers.")
            }
            if p.immobile_atom == nil {
                throw VinaError("Atom \(to) has not been found in this branch.")
            }
            return
        }
        else if line.hasPrefix("ATOM") || line.hasPrefix("HETATM") {
            let a = try parse_pdbqt_atom_string(str: line)
            if a.number == to {
                p.immobile_atom = sz(p.atoms.count)
            }
            p.add(a: a, c: c)
        }
        else if line.hasPrefix("MODEL") {
            throw VinaError("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.")
        }
        else {
            throw VinaError("Unknown or inappropriate tag found in flex residue or ligand.")
        }
    }
}

class PDBQTInitializer {
    
    let m: Model
    let atom_typing_used: AtomType.T
    
    init(atom_typing_used: AtomType.T) {
        self.atom_typing_used = atom_typing_used
        self.m = Model(atom_typing_used)
    }
    
    func initialize_from_rigid(_ r: Rigid) throws {
        guard m.grid_atoms.isEmpty else {
            throw VinaError("Model grid atoms are not emtpy when attempting initialization")
        }
//        m.grid_atoms = ContiguousArray(r.atoms.map { $0.copy() as! Atom })
        m.grid_atoms = r.atoms
    }
    
    func initialize_from_nrp(nrp: NonRigidParsed, c: Context, is_ligand: Bool) {
        assert(m.ligands.elements.isEmpty)
        assert(m.flex.elements.isEmpty)
        
//        m.ligands = nrp.ligands //.copy() as! VectorMutable<Ligand, LigandConf, LigandChange>
        m.ligands = nrp.ligands.copy() as! VectorMutable<Ligand, LigandConf, LigandChange>
//        m.flex = nrp.flex //.copy() as! VectorMutable<Residue, ResidueConf, ResidueChange>
        m.flex = nrp.flex.copy() as! VectorMutable<Residue, ResidueConf, ResidueChange>

        assert(m.atoms.isEmpty)

        let n = nrp.atoms.count + nrp.inflex.count

        for i in 0..<nrp.atoms.count {
            let a: MoveableAtom = nrp.atoms[i] as MoveableAtom
            let b = a.copy() as! Atom
            b.coords = a.relative_coords
            m.atoms.append(b)
            m.coords.append(a.coords)
        }

        for i in 0..<nrp.inflex.count {
            let a: Atom = nrp.inflex[i]
            let b = a.copy() as! Atom
            b.coords = zero_vec // to avoid any confusion; presumably these will never be looked at
            m.atoms.append(b)
            m.coords.append(a.coords)
        }

        assert(m.coords.count == n)

        m.minus_forces = m.coords
        m.m_num_movable_atoms = sz(nrp.atoms.count)

        if is_ligand {
            assert(m.ligands.elements.count == 1)
            m.ligands.elements[0].cont = c
        } else {
            m.flex_context = c
        }
        
    }
    
    func initialize(mobility: DistanceTypeMatrix) throws {
        try m.initialize(mobility: mobility)
    }
}
