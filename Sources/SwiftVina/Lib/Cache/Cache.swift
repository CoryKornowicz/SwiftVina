//
//  Cache.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

func convert_XS_to_string(_ t: sz) -> String {
    switch(t) {
        case XS_TYPE_C_H     : return "C_H"
        case XS_TYPE_C_P     : return "C_P"
        case XS_TYPE_N_P     : return "N_P"
        case XS_TYPE_N_D     : return "N_D"
        case XS_TYPE_N_A     : return "N_A"
        case XS_TYPE_N_DA    : return "N_DA"
        case XS_TYPE_O_P     : return "O_P"
        case XS_TYPE_O_D     : return "O_D"
        case XS_TYPE_O_A     : return "O_A"
        case XS_TYPE_O_DA    : return "O_DA"
        case XS_TYPE_S_P     : return "S_P"
        case XS_TYPE_P_P     : return "P_P"
        case XS_TYPE_F_H     : return "F_H"
        case XS_TYPE_Cl_H    : return "Cl_H"
        case XS_TYPE_Br_H    : return "Br_H"
        case XS_TYPE_I_H     : return "I_H"
        case XS_TYPE_Si      : return "Si"
        case XS_TYPE_At      : return "At"
        case XS_TYPE_Met_D   : return "Met_D"
        case XS_TYPE_W       : return "W"
        default: assert(false)
    }
    // Should not hit this branch
    assert(false)
    return ""
}

func read_vina_map(filename: URL, gds: inout Array<grid_dims>, g: inout Grid) throws {
    var pt_counter: Int = 0
    var x: sz = 0
    var y: sz = 0
    var z: sz = 0
    var gd: grid_dims = grid_dims(repeating: grid_dim(), count: 3)
    var spacing: fl = 0
    var center: fl = 0
    var halfspan: fl = 0
    
    // open a read filename URL line by line
    try filename.foreachRow { line, line_counter in // offset by one from original impl.
        
        if line_counter == 3 { // spacing
            let fields = line.split(separator: " ")
            spacing = fl((fields[1] as NSString).doubleValue)
        }
        
        if line_counter == 4  { // n_voxels
            let fields = line.split(separator: " ")
            for i in 0..<3 {
                // n_voxels must be EVEN
                // because the number of sampled points in the grid is always ODD
                // (number of sampled points == n_voxels + 1)
                gd[i].n_voxels = sz((fields[i + 1] as NSString).integerValue)
                if gd[i].n_voxels % 2 == 1 {
                    throw VinaError("ERROR: number of voxels (NELEMENTS) must be even")
                }
            }
        }
        
        if line_counter == 5 { // center
            let fields = line.split(separator: " ")
            for i in 0..<3 {
                center = fl((fields[i + 1] as NSString).doubleValue)
                halfspan = fl((gd[i].n_voxels)) * spacing / 2.0
                gd[i].begin = center - halfspan
                gd[i].end = center + halfspan
                // print("center: \(center) halfspan: \(halfspan) begin: \(gd[i].begin) end: \(gd[i].end)")
            }
            gds.append(gd)
            g = .init(gd: gd)
        }
        
        if line_counter > 5 { // data_points
            print("\(pt_counter) \(x) \(y) \(z) \(String(describing: fl(line)))")
            g.m_data[x, y, z] = fl(line)!
            y += x == (gd[0].n_voxels + 1) ? 1 : 0
            z += y == (gd[1].n_voxels + 1) ? 1 : 0
            x = x % (gd[0].n_voxels + 1)
            y = y % (gd[1].n_voxels + 1)
            pt_counter += 1
            x += 1
        }
    }
}

final class Cache: igrid {
    
    private var m_gd: grid_dims
    private var m_slope: fl = 0
    private var m_grids: ContiguousArray<Grid> = ContiguousArray<Grid>.init(repeating: Grid(), count: Int(XS_TYPE_SIZE))
    
    var gd: grid_dims {
        return m_gd
    }
    
    init(_ slope: fl = 1e6) {
        self.m_slope = slope
        self.m_gd = grid_dims(repeating: grid_dim(), count: 3)
    }

    init(_ gd: grid_dims, _ slope: fl = 1e6) {
        self.m_gd = gd
        self.m_slope = slope
    }
    
    func corner1() -> vec {
        return vec(m_gd[0].begin, m_gd[1].begin, m_gd[2].begin)
    }
    
    func corner2() -> vec {
        return vec(m_gd[0].end, m_gd[1].end, m_gd[2].end)
    }
    
    func is_in_grid(_ m: Model, _ margin: fl = 0.0001) -> Bool {
        
        for i in 0..<m.num_movable_atoms {
            if(m.atoms[i].is_hydrogen()) { continue }
            let a_coords: vec = m.coords[i]
            for j in 0..<m_gd.count {
                if(m_gd[j].n_voxels > 0) {
                    if(a_coords[j] < m_gd[j].begin - margin || a_coords[j] > m_gd[j].end + margin) {
                        return false
                    }
                }
            }
        }
        
        return true
    }
    
    func is_atom_type_grid_initialized(_ t: sz) -> Bool {
        return m_grids[t].initialized
    }
    
    func are_atom_types_grid_initialized(_ atom_types: szv) -> Bool {
        let nat: sz = num_atom_types(.XS)
        for i in 0..<atom_types.count {
            var t: sz = atom_types[i]
            if(t >= nat) { continue }
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
            XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
            XS_TYPE_C_H_CG3:
                t = XS_TYPE_C_H
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
            XS_TYPE_C_P_CG3:
                t = XS_TYPE_C_P
                break
            default:
                break
            }
            
            if(!is_atom_type_grid_initialized(t)) {
                print("ERROR: Affinity map for atom type \(convert_XS_to_string(t)) is not present.")
                return false
            }
        }
        return true
    }
    
    func eval(m: Model, v: fl) -> fl {
        var e: fl = 0
        let nat: sz = num_atom_types(.XS)
        
        for i in 0..<m.num_movable_atoms {
            var t: sz = m.atoms[i].xs
            
            if(t >= nat) { continue }
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
            XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
            XS_TYPE_C_H_CG3:
                t = XS_TYPE_C_H
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
            XS_TYPE_C_P_CG3:
                t = XS_TYPE_C_P
                break
            default:
                break
            }
            
            e += m_grids[t].evaluate(location: m.coords[i], slope: m_slope, c: v)
        }
        return e
    }
    
    func eval_intra(m: Model, v: fl) -> fl {
        var e: fl = 0
        let nat: sz = num_atom_types(.XS)
        
        for i in 0..<m.num_movable_atoms {
            if(m.is_atom_in_ligand(i)) { continue } // we only want flex-rigid interaction
            var t: sz = m.atoms[i].xs
            
            if(t >= nat) { continue }
            
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
            XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
            XS_TYPE_C_H_CG3:
                t = XS_TYPE_C_H
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
            XS_TYPE_C_P_CG3:
                t = XS_TYPE_C_P
                break
            default:
                break
            }
            
            e += m_grids[t].evaluate(location: m.coords[i], slope: m_slope, c: v)
        }
        return e
    }
    
    func eval_deriv(m: Model, v: fl) -> fl {
        var e: fl = 0
        let nat: sz = num_atom_types(.XS)
        
        for i in 0..<m.num_movable_atoms {
            var t: sz = m.atoms[i].xs
            
            if(t >= nat) { continue }
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
            XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
            XS_TYPE_C_H_CG3:
                t = XS_TYPE_C_H
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
            XS_TYPE_C_P_CG3:
                t = XS_TYPE_C_P
                break
            default:
                break
            }
            
            e += m_grids[t].evaluate(location: m.coords[i], slope: m_slope, c: v, deriv: &m.minus_forces[i])
        }
        
        return e
    }
    
    
    func read(_ map_prefix: String) throws {
        let nat = num_atom_types(.XS)
        var type: String, filename: String
        var gds: Array<grid_dims> = []  // to check all maps have same dims (TODO)
        
        var got_C_H_already: Bool = false
        var got_C_P_already: Bool = false
        var found_at_least_1_map: Bool = false
        
        for atom_type in 0..<XS_TYPE_SIZE {
            var t: sz = atom_type
            
            if(t >= nat) { continue }
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
            XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
            XS_TYPE_C_H_CG3:
                if(got_C_H_already) { continue }
                t = XS_TYPE_C_H
                got_C_H_already = true
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
            XS_TYPE_C_P_CG3:
                if(got_C_P_already) { continue }
                t = XS_TYPE_C_P
                got_C_P_already = true
                break
            default:
                break
            }
            
            type = convert_XS_to_string(t)
            //MARK: this might be bad because what if the file is not in the local scope?
            filename = map_prefix + "." + type + ".map"
            
            if FileManager.default.fileExists(atPath: filename) {
                // convert to URL?
                let p = URL(filePath: filename)
                try read_vina_map(filename: p, gds: &gds, g: &m_grids[t])
                found_at_least_1_map = true
            } else {
                throw VinaError("File not found at \(filename)")
            }
        } // map loop
        
        if !found_at_least_1_map {
            throw VinaError("\nERROR: No *.map files with prefix \"\(map_prefix)\"")
        }
        
        // Store in Cache object
        m_gd = gds[0]
    }
    
    
    func write(out_prefix: String, atom_types: szv, gpf_filename: String?,
               fld_filename: String?, receptor_filename: String?) throws {
        
        let nat = num_atom_types(.XS)
        var atom_type: String
        var filename: String
        var got_C_H_already: Bool = false
        var got_C_P_already: Bool = false
        
        for i in 0..<atom_types.count {
            var t: sz = atom_types[i]
            
            if(t >= nat) { continue }
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
                XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
                XS_TYPE_C_H_CG3:
                if(got_C_H_already) { continue }
                t = XS_TYPE_C_H
                got_C_H_already = true
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
                XS_TYPE_C_P_CG3:
                if(got_C_P_already) { continue }
                t = XS_TYPE_C_P
                got_C_P_already = true
                break
            default:
                break
            }
            
            if(m_grids[t].initialized) {
                let nx = m_grids[t].m_data.dim0() - 1
                let ny = m_grids[t].m_data.dim1() - 1
                let nz = m_grids[t].m_data.dim2() - 1
                atom_type = convert_XS_to_string(t)
                
                // For writing a .map, the number of points in the grid must be odd, which means that the number
                // of voxels (NELEMENTS) must be even... n_voxels = n_grid_points - 1
                if(nx % 2 == 1 || ny % 2 == 1 || nz % 2 == 1) {
                    throw VinaError("ERROR: Can't write maps. Number of voxels (NELEMENTS) is odd. Use --force_even_voxels.")
                    
                } else {
                    //TODO: Check if the file exists through FileManager first
                    filename = out_prefix + "." + atom_type + ".map"
                    //                     let p = URL(filePath: filename)
                    var out: TextStreamer = TextStreamer(filepath: filename)
                    
                    // write header
                    print("GRID_PARAMETER_FILE \(String(describing: gpf_filename))", to: &out)
                    print("GRID_DATA_FILE \(String(describing: fld_filename))", to: &out)
                    print("MACROMOLECULE \(String(describing: receptor_filename))", to: &out)
                    print("SPACING \(m_grids[t].m_factor_inv[0])", to: &out)
                    
                    print("NELEMENTS \(nx) \(ny) \(nz)", to: &out)
                    
                    // center
                    let cx = m_grids[t].m_init[0] + m_grids[t].m_range[0] * 0.5
                    let cy = m_grids[t].m_init[1] + m_grids[t].m_range[1] * 0.5
                    let cz = m_grids[t].m_init[2] + m_grids[t].m_range[2] * 0.5
                    print("CENTER \(cx) \(cy) \(cz)", to: &out)
                    
                    // write data
                    for z in 0..<m_grids[t].m_data.dim2() {
                        for y in 0..<m_grids[t].m_data.dim1() {
                            for x in 0..<m_grids[t].m_data.dim0() {
                                print(String(format: "%.4f", m_grids[t].m_data[x, y, z]), to: &out)
                            } // x
                        } // y
                    } // z
                } // even voxels
            } // map initialized
        } // map loop
    } // write
    
    func populate(m: Model, p: precalculate, atom_types_needed: szv) {
        var needed: szv = []
        var got_C_H_already: Bool = false
        var got_C_P_already: Bool = false
    
        for i in 0..<atom_types_needed.count {
            var t: sz = atom_types_needed[i]
            switch(t) {
            case XS_TYPE_G0,
                XS_TYPE_G1,
                XS_TYPE_G2,
                XS_TYPE_G3:
                continue
            case XS_TYPE_C_H_CG0,
                XS_TYPE_C_H_CG1,
                XS_TYPE_C_H_CG2,
                XS_TYPE_C_H_CG3:
                if(got_C_H_already) { continue }
                t = XS_TYPE_C_H
                got_C_H_already = true
                break
            case XS_TYPE_C_P_CG0,
                XS_TYPE_C_P_CG1,
                XS_TYPE_C_P_CG2,
                XS_TYPE_C_P_CG3:
                if(got_C_P_already) { continue }
                t = XS_TYPE_C_P
                got_C_P_already = true
                break
            default:
                break
            }
            
            if(!m_grids[t].initialized) {
                needed.append(t)
                m_grids[t] = .init(gd: m_gd)
            }
        }

        if needed.isEmpty { return }

        var affinities: flv = ContiguousArray<fl>.init(repeating: 0, count: needed.count)
        let nat: sz = num_atom_types(.XS)
        let g: Grid = m_grids[needed[0]]
        let cutoff_sqr: fl = p.cutoff_sqr()
        let gd_reduced: grid_dims = szv_grid_dims(m_gd)
        let ig: SZVGrid = SZVGrid(m, gd_reduced, cutoff_sqr)
        
        for x in 0..<g.m_data.dim0() {
            for y in 0..<g.m_data.dim1() {
                for z in 0..<g.m_data.dim2() {
                    affinities = ContiguousArray<fl>.init(repeating: 0, count: needed.count)
                    let probe_coords: vec = g.index_to_argument(x: x, y: y, z: z)
                    let possibilities: szv = ig.possibilities(probe_coords)
                    for possibilities_i in 0..<possibilities.count {
                        let i: sz = possibilities[possibilities_i]
                        let a: Atom = m.grid_atoms[i]
                        let t1: sz = a.xs
                        if(t1 >= nat) { continue }
                        let r2: fl = vec_distance_sqr(a.coords, probe_coords)
                        if(r2 <= cutoff_sqr) {
                            for j in 0..<needed.count {
                                let t2: sz = needed[j]
                                assert(t2 < nat)
                                let type_pair_index: sz = triangular_matrix_index_permissive(nat, t1, t2)
                                affinities[j] += p.eval_fast(type_pair_index: type_pair_index, r2: r2)
                            }
                        }
                    }
                    for j in 0..<needed.count {
                        let t: sz = needed[j]
                        assert(t < nat)
                        m_grids[t].m_data[x, y, z] = affinities[j]
                    }
                }
            }
        }
     }
}
