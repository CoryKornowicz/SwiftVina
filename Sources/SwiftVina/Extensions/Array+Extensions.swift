//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation

// TODO: Optimize these, they seem a bit excessive...

extension Array {
    
//    @inlinable public var count: sz {
//        sz(self.endIndex - self.startIndex)
//    }
    
    @inline(__always)
    subscript(_ i: UInt32) -> Element {
        get {
            self[Int(i)]
        }
        
        set {
            self[Int(i)] = newValue
        }
    }
    
    @inline(__always)
    subscript(_ i: UInt64) -> Element {
        get {
            self[Int(i)]
        }
        
        set {
            self[Int(i)] = newValue
        }
    }
    
    
//    subscript(unsafe index: UInt64) -> Element {
//        // Check if the UInt64 index can be converted to Int
//        get {
//            return withUnsafeBufferPointer { buffer -> Element in
//                // Safely access the buffer pointer and return the element
//                guard let baseAddress = buffer.baseAddress else {
//                    fatalError("Unsafe access buffer address cannot be found")
//                }
//                return baseAddress.advanced(by: Int(index)).pointee
//            }
//        }
//        
//        set {
//            self[unsafe: index] = newValue
//        }
//    }
    
//    subscript(safe index: UInt64) -> Element? {
//        // Check if the UInt64 index can be converted to Int
//        guard let intIndex = Int(exactly: index) else {
//            return nil
//        }
//        // Ensure the index is within the bounds of the array
//        guard intIndex >= 0 && intIndex < count else {
//            return nil
//        }
//        return withUnsafeBufferPointer { buffer -> Element? in
//            // Safely access the buffer pointer and return the element
//            guard let baseAddress = buffer.baseAddress else { return nil }
//            return baseAddress.advanced(by: intIndex).pointee
//        }
//    }
//    
//    subscript(modifiable index: UInt64) -> Element? {
//        get {
//            return self[safe: index]
//        }
//        set {
//            guard newValue != nil else { return }
//            // Check if the UInt64 index can be converted to Int
//            guard let intIndex = Int(exactly: index) else {
//                return
//            }
//            // Ensure the index is within the bounds of the array
//            guard intIndex >= 0 && intIndex < count else {
//                return
//            }
//            self[modifiable: index] = newValue
////            self.modify(at: intIndex) { ele in
////                ele = newValue! // we guard'd this from being nil in this scope ^^^
////            }
//        }
//    }
    
    // MARK: Int extensions below
    
//    subscript(safe index: Int) -> Element? {
//        // Ensure the index is within the bounds of the array
//        guard index >= 0 && index < count else {
//            return nil
//        }
//        return withUnsafeBufferPointer { buffer -> Element? in
//            // Safely access the buffer pointer and return the element
//            guard let baseAddress = buffer.baseAddress else { return nil }
//            return baseAddress.advanced(by: index).pointee
//        }
//    }
//
//    subscript(modifiable index: Int) -> Element? {
//        get {
//            return self[safe: index]
//        }
//        set {
//            guard newValue != nil else { return }
//            // Ensure the index is within the bounds of the array
//            guard index >= 0 && index < count else {
//                return
//            }
//            self[index] = newValue
////            self.modify(at: index) { ele in
////                ele = newValue! // we guard'd this from being nil in this scope ^^^
////            }
//        }
//    }
}

extension Array {
    /// Extends the array to a specified size with a default value.
    /// If the array is already equal to or larger than the specified size, no change is made.
    /// - Parameters:
    ///   - size: The new size of the array.
    ///   - defaultValue: The value to fill new elements with.
    mutating func extend(to size: Int, with defaultValue: Element) {
        guard size > count else { return }
        append(contentsOf: repeatElement(defaultValue, count: size - count))
    }
    
    mutating func withLastElement(_ closure: (inout Element) throws -> Void) rethrows {
        try closure(&self[self.endIndex - 1])
    }
    
}


extension ContiguousArray {
    
//    @inlinable public var count: sz {
//        sz(self.endIndex - self.startIndex)
//    }
    
    @inline(__always)
    subscript(_ i: UInt32) -> Element {
        get {
            self[Int(i)]
        }
        set {
            self[Int(i)] = newValue
        }
    }
    
    @inline(__always)
    subscript(_ i: UInt64) -> Element {
        get {
            self[Int(i)]
        }
        set {
            self[Int(i)] = newValue
        }
    }
    
    mutating func withLastElement(_ closure: (inout Element) throws -> Void) rethrows {
        try closure(&self[self.endIndex - 1])
    }
    
}

extension ContiguousArray {
    /// Extends the array to a specified size with a default value.
    /// If the array is already equal to or larger than the specified size, no change is made.
    /// - Parameters:
    ///   - size: The new size of the array.
    ///   - defaultValue: The value to fill new elements with.
    mutating func extend(to size: Int, with defaultValue: Element) {
        guard size > count else { return }
        append(contentsOf: repeatElement(defaultValue, count: size - count))
    }
}
