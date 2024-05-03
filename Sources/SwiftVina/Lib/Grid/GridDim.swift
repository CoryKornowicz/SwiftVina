//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

// MARK: Will always be 3D, maybe make this into a concrete type to always have three dims...
typealias grid_dims = ContiguousArray<grid_dim>

struct grid_dim: Equatable {
    var begin: fl = 0
    var end: fl = 0
    var n_voxels: sz = 0 // number of intervals == number of sample points - 1
    
    func span() -> fl {
        return end - begin
    }
    
    func enabled() -> Bool {
        return n_voxels > 0
    }

    static func ==(lhs: grid_dim, rhs: grid_dim) -> Bool {
        return lhs.n_voxels == rhs.n_voxels && eq(lhs.begin, rhs.begin) && eq(lhs.end, rhs.end)
    }
}

@inline(__always)
func ==(lhs: grid_dims, rhs: grid_dims) -> Bool {
    return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]
}

func print(_ gd: grid_dims) {
    for i in 0..<gd.count {
        print(gd[i].n_voxels, terminator: " ")
        print("[\(gd[i].begin) .. \(gd[i].end)]")
    }
    print("\n")
}

@inline(__always)
func grid_dims_begin(_ gd: grid_dims) -> vec {
//    var tmp: vec = zero_vec
//    for i in 0..<gd.count {
//        tmp[i] = gd[i].begin
//    }
//    return tmp
    return vec(gd[0].begin, gd[1].begin, gd[2].begin)
}

@inline(__always)
func grid_dims_end(_ gd: grid_dims) -> vec {
//    var tmp: vec = zero_vec
//    for i in 0..<gd.count {
//        tmp[i] = gd[i].end
//    }
//    return tmp
    return vec(gd[0].end, gd[1].end, gd[2].end)
}

