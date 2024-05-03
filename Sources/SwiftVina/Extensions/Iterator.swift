//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation

final class Iterator<T>: IteratorProtocol, Sequence {
    
    private var collection: ContiguousArray<T>
    private var index = 0

    init(_ collection: ContiguousArray<T>) {
        self.collection = collection
    }

    var start: T {
        get {
            return collection[0]
        }
        set {
            collection[0] = newValue
        }
    }
    
    var curr: T {
        get {
            assert(index < collection.count)
            assert(index >= 0)
            return collection[index]
        }
        set {
            collection[index] = newValue
        }
    }
    
    subscript(_ index: Int) -> T {
        get {
            guard index < self.collection.count else {
                fatalError("index out of bounds")
            }
            return collection[index]
        }
        set {
            collection[index] = newValue
        }
    }
    
    @discardableResult
    public func next() -> T? {
        defer { index += 1 }
        return index < collection.count ? collection[index] : nil
    }
    
    @discardableResult
    public func next() -> T {
        assert(index < collection.count)
        assert(index + 1 <= collection.count)
        defer { index += 1 }
        return collection[index]
    }
    
    static func += (_ lhs: Iterator<T>, _ rhs: Int) {
        lhs.index += rhs
    }
    
    static func -= (_ lhs: Iterator<T>, _ rhs: Int) {
        lhs.index -= rhs
    }
    
    @discardableResult
    static prefix func ++(_ lhs: Iterator<T>) -> T {
        lhs.index += 1
        return lhs[lhs.index]
    }
    
}
