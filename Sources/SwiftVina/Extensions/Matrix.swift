//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation
import simd

@usableFromInline
final class matrix<T> {
    private var m_data: Array<T>
    private var m_i, m_j: sz
    
    init() {
        self.m_data = []
        self.m_i = 0
        self.m_j = 0
    }
    
    init(i: sz, j: sz, filler_val: T) {
        self.m_data = Array<T>(repeating: filler_val, count: Int(i*j))
        self.m_i = i
        self.m_j = j
    }
    
    @usableFromInline
    func index(i: sz, j: sz) -> sz {
        assert(j < m_j)
        assert(i < m_i)
        return i + m_i*j // column-major
    }
    
    func resize(m: sz, n: sz, filler_val: T) { // new sizes should be the same or greater than the old
        if m == dim_1() && n == dim_2() {
            return // no-op
        }
        assert(m >= dim_1())
        assert(n >= dim_2())
        var tmp = Array<T>(repeating: filler_val, count: Int(m*n))
        for i in 0..<m_i {
            for j in 0..<m_j {
                tmp[i+m*j] = self[i, j]
            }
        }
        m_data = tmp
        m_i = m
        m_j = n
    }
    
    func append(_ x: matrix<T>, _ filler_val: T) {
        let m = dim_1()
        let n = dim_2()
        resize(m: m + x.dim_1(), n: n + x.dim_2(), filler_val: filler_val)
        for i in 0..<x.dim_1() {
            for j in 0..<x.dim_2() {
                self[i+m, j+n] = x[i, j]
            }
        }
    }
    
    @usableFromInline
    subscript(_ i: sz) -> T {
        get {
            return m_data[i]
        }
        
        set {
            m_data[i] = newValue
        }
    }
    
    @usableFromInline
    subscript(_ i: sz, _ j: sz) -> T {
        get {
            return m_data[index(i: i, j: j)]
        }
        
        set {
            m_data[index(i: i, j: j)] = newValue
        }
    }

    func dim_1() -> sz {
        return m_i
    }

    func dim_2() -> sz {
        return m_j
    }
 }

@usableFromInline
final class triangular_matrix<T> {
    private var m_data: ContiguousArray<T>
    private var m_dim: sz
    
    init(n: sz, filler_val: T) {
        self.m_data = ContiguousArray<T>(repeating: filler_val, count: Int(n*(n+1)/2))
        self.m_dim = n
    }
    
    @usableFromInline
    func index(i: sz, j: sz) -> sz {
        return triangular_matrix_index(m_dim, i, j)
    }
    
    @usableFromInline
    func index_permissive(i: sz, j: sz) -> sz {
        return (i < j) ? index(i: i, j: j) : index(i: j, j: i)
    }
    
    @usableFromInline
    func dim() -> sz {
        return m_dim
    }

    @usableFromInline
    subscript(_ i: sz) -> T {
        get {
            m_data[i]
        }
        
        set {
            m_data[i] = newValue
        }
        
    }

    @usableFromInline
    subscript(_ i: sz, _ j: sz) -> T {
        get {
            m_data[index(i: i, j: j)]
        }
        set {
            m_data[index(i: i, j: j)] = newValue
        }
    }
 }

@usableFromInline
class strictly_triangular_matrix<T> {
    
    private var m_data: Array<T>
    private var m_dim: sz
    
    init() {
        self.m_data = []
        self.m_dim = 0
    }
    
    init(n: sz, filler_val: T) {
        self.m_data = Array<T>(repeating: filler_val, count: Int(n*(n-1)/2))
        self.m_dim = n
    }
    
    @usableFromInline
    func index(i: sz, j: sz) -> sz {
        assert(j < m_dim)
        assert(i < j)
        assert(j >= 1)
        return i + j*(j-1)/2
    }
    
    @usableFromInline
    func index_permissive(i: sz, j: sz) -> sz {
        return (i < j) ? index(i: i, j: j) : index(i: j, j: i)
    }
    
    func dim() -> sz {
        return m_dim
    }
    
    @usableFromInline
    subscript(_ i: sz) -> T {
        get {
            m_data[i]
        }
        
        set {
            m_data[i] = newValue
        }
    }

    @usableFromInline
    subscript(_ i: sz, _ j: sz) -> T {
        get {
            m_data[index(i: i, j: j)]
        }
        
        set {
            m_data[index(i: i, j: j)] = newValue
        }
    }

    func resize(n: sz, filler_val: T) {
        if n == m_dim {
            return
        }
        assert(n > m_dim)
        m_dim = n
        m_data.extend(to: Int(n*(n-1)/2), with: filler_val) // preserve original data
    }

    func append(_ m: strictly_triangular_matrix<T>, _ filler_val: T) {
        let n = dim()
        resize(n: n + m.dim(), filler_val: filler_val)
        for i in 0..<m.dim() {
            for j in i+1..<m.dim() {
                self[i+n, j+n] = m[i, j]
            }
        }
    }
    
    func append(rectangular: matrix<T>, triangular: strictly_triangular_matrix<T>) {
        assert(dim() == rectangular.dim_1())
        assert(rectangular.dim_2() == triangular.dim())

        // a filler value is needed by append or resize
        // we will use a value from rectangular as the filler value
        // but it can not be obtained if dim_1 or dim_2 is 0
        // these cases have to be considered separately
        if rectangular.dim_2() == 0 { return }
        if rectangular.dim_1() == 0 {
            // self = triangular can't do this in Swift -_-
            self.m_dim = triangular.m_dim
            self.m_data = triangular.m_data
            return
        }
        let filler_val = rectangular[0, 0]
        let n = dim()
        append(triangular, filler_val)
        for i in 0..<rectangular.dim_1() {
            for j in 0..<rectangular.dim_2() {
                self[i, n+j] = rectangular[i, j]
            }
        }
    }
}
