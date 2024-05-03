//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation

@inline(__always)
func triangular_matrix_index(_ n: sz, _ i: sz, _ j: sz) -> sz {
//    assert(j < n)
//    assert(i <= j)
    return i + j*(j+1)/2
}

@inline(__always)
func triangular_matrix_index_permissive(_ n: sz, _ i: sz, _ j: sz) -> sz {
    return (i <= j) ? triangular_matrix_index(n, i, j)
                    : triangular_matrix_index(n, j, i)
}
