//
//  SZVGrid.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

func szv_grid_dims(_ gd: grid_dims) -> grid_dims {
    var tmp: grid_dims = ContiguousArray(repeating: grid_dim(), count: 3)
    for i in 0..<gd.count {
        tmp[i].begin = gd[i].begin
        tmp[i].end = gd[i].end
        let n_fl = (gd[i].end - gd[i].begin) / 3 // 3A preferred size
        let n_int = Int(n_fl)
        tmp[i].n_voxels = (n_int < 1) ? 1 : sz(n_int)
    }
    return tmp
}

final class SZVGrid {
    var m_data: Array3D<szv>
    var m_init: vec = zero_vec
    var m_range: vec = zero_vec

    init(_ m: Model, _ gd: grid_dims, _ cutoff_sqr: fl) {
        self.m_data = Array3D<szv>(i: gd[0].n_voxels, j: gd[1].n_voxels, k: gd[2].n_voxels, szv())
        
        var end: vec = zero_vec
        for i in 0..<gd.count {
            self.m_init[i] = gd[i].begin
            end[i] = gd[i].end
        }
        self.m_range = end - m_init

        let nat = num_atom_types(m.atom_typing_used)

        var relevant_indexes: szv = []
        for i in 0..<m.grid_atoms.count {
            let a = m.grid_atoms[i]
            if a.get(m.atom_typing_used) < nat && brick_distance_sqr(begin: m_init, end: end, v: a.coords) < cutoff_sqr {
                relevant_indexes.append(sz(i))
            }
        }

        for x in 0..<m_data.dim0() {
            for y in 0..<m_data.dim1() {
                for z in 0..<m_data.dim2() {
                    for ri in 0..<relevant_indexes.count {
                        let i = relevant_indexes[ri]
                        let a = m.grid_atoms[i]
                        if brick_distance_sqr(begin: index_to_coord(x, y, z), end: index_to_coord(x+1, y+1, z+1), v: a.coords) < cutoff_sqr {
                            m_data[x, y, z].append(i)
                        }
                    }
                }
            }
        }
    }

    func possibilities(_ coords: vec) -> szv {
        var index: Array<sz> = [0, 0, 0]
        for i in 0..<3 {
            assert(coords[i] + epsilon_fl >= m_init[i])
            assert(coords[i] <= m_init[i] + m_range[i] + epsilon_fl)
            let tmp = (coords[i] - m_init[i]) * fl(m_data.dim(sz(i))) / m_range[i]
            index[i] = fl_to_sz(tmp, max_sz: m_data.dim(sz(i)) - 1)
        }
        return m_data[index[0], index[1], index[2]]
    }

    func average_num_possibilities() -> fl {
        var counter: sz = 0
        for x in 0..<m_data.dim0() {
            for y in 0..<m_data.dim1() {
                for z in 0..<m_data.dim2() {
                    counter += sz(m_data[x, y, z].count)
                }
            }
        }
        return fl(counter) / (fl(m_data.dim0()) * fl(m_data.dim1()) * fl(m_data.dim2()))
    }

    private func index_to_coord(_ i: sz, _ j: sz, _ k: sz) -> vec {
        let index = vec(fl(i), fl(j), fl(k))
        var tmp: vec = zero_vec
        for n in 0..<3 {
            tmp[n] = m_init[n] + m_range[n] * index[n] / fl(m_data.dim(sz(n)))
        }
        return tmp
    }
}
