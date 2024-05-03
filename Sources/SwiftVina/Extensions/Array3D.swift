//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

@usableFromInline
final class Array3D<T> {
    
    var m_data: ContiguousArray<T>
    fileprivate var m_i, m_j, m_k: sz

    init() {
        self.m_i = 0
        self.m_j = 0
        self.m_k = 0
        self.m_data = []
    }

    init(i: sz, j: sz, k: sz, _ initialValue: T) {
        self.m_i = i
        self.m_j = j
        self.m_k = k
        self.m_data = ContiguousArray<T>(repeating: initialValue, count: Int(i*j*k))
    }
    
    @inline(__always)
    func dim0() -> sz {
        return m_i
    }

    @inline(__always)
    func dim1() -> sz {
        return m_j
    }

    @inline(__always)
    func dim2() -> sz {
        return m_k
    }
    
    @inline(__always)
    func dim(_ i: sz) -> sz {
        switch(i) {
            case 0: return m_i
            case 1: return m_j
            case 2: return m_k
            default: assert(false); return 0 // to get rid of the warning
        }
    }

    func resize(_ i: sz, _ j: sz, _ k: sz, fillValue: T) {
        m_i = i
        m_j = j
        m_k = k
        m_data.extend(to: Int(i*j*k), with: fillValue)
    }

    @inline(__always)
    subscript(_ i: sz, _ j: sz, _ k: sz) -> T {
        get {
            return m_data[i + m_i*(j + m_j*k)]
        }
        set {
            m_data[i + m_i*(j + m_j*k)] = newValue
        }
    }

    func print() {
        for i in 0..<m_i {
            for j in 0..<m_j {
                for k in 0..<m_k {
                    Swift.print(self[i, j, k], terminator: " ")
                }
                Swift.print("\n", terminator: "")
            }
            Swift.print("\n", terminator: "")
        }
    }
    
    
}
