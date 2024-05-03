//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/16/23.
//

import Foundation

extension String {
    
    var lines: [String] {
        return self.components(separatedBy: "\n")
    }
    
    /**
     Enables passing in negative indices to access characters
     starting from the end and going backwards.
     if num is negative, then it is added to the
     length of the string to retrieve the true index.
     */
    func negativeIndex(_ num: Int) -> Int {
        return num < 0 ? num + self.count : num
    }
    
    func strOpenRange(index i: Int) -> Range<String.Index> {
        let j = negativeIndex(i)
        return strOpenRange(j..<(j + 1), checkNegative: false)
    }
    
    func strOpenRange(
        _ range: Range<Int>, checkNegative: Bool = true
    ) -> Range<String.Index> {
        
        var lower = range.lowerBound
        var upper = range.upperBound
        
        if checkNegative {
            lower = negativeIndex(lower)
            upper = negativeIndex(upper)
        }
        
        let idx1 = index(self.startIndex, offsetBy: lower)
        let idx2 = index(self.startIndex, offsetBy: upper)
        
        return idx1..<idx2
    }
    
    // MARK: - Subscripts
    
    /**
     Gets and sets a character at a given index.
     Negative indices are added to the length so that
     characters can be accessed from the end backwards
     
     Usage: `string[n]`
     */
    subscript(index i: Int) -> String {
        get {
            return String(self[strOpenRange(index: i)])
        }
        set {
            let range = strOpenRange(index: i)
            replaceSubrange(range, with: newValue)
        }
    }
    
    
    subscript<T: RangeExpression>(slow bounds: T) -> String where T.Bound == Int {
        let range = bounds.relative(to: 0 ..< .max)
        let left = index(startIndex, offsetBy: range.lowerBound)
        let right = range.upperBound >= .max ? endIndex : index(startIndex, offsetBy: range.upperBound)
        return String(self[left ..< right])
    }
}

//extension String: RowIterable {
//    
//    func foreachRow(skipLines: Int, _ rowParcer: ((String, Int) -> Void)) {
//        for (i, v) in lines[skipLines...].enumerated() {
//            rowParcer(v, i)
//        }
//    }
//    
//    func foreachRow(_ rowParcer: ((String, Int) -> Void)) {
//        for (i, v) in lines.enumerated() {
//            rowParcer(v, i)
//        }
//    }
//    
//    func foreachRow(_ rowParcer: ((String, Int) throws -> Void)) rethrows {
//        for (i, v) in lines.enumerated() {
//            try rowParcer(v, i)
//        }
//    }
//
//}
