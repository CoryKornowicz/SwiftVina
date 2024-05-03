//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/18/23.
//

import Foundation

class RowIterator: RowIterable {
    
    typealias T = String
    var lineOffset: Int = 0
    var lines: [String]
    
    init(url: URL) {
        do {
            let urlContentsString = try NSString(contentsOf: url, encoding: NSUTF8StringEncoding) as String
            self.lines = urlContentsString.lines
        } catch {
            self.lines = []
        }
    }
    
    init(filePath: String) {
        do {
            let urlContentsString = try NSString(contentsOfFile: filePath, encoding: NSUTF8StringEncoding) as String
            self.lines = urlContentsString.lines
        } catch {
            self.lines = []
        }
    }
    
    init(stringContents: String) {
        self.lines = stringContents.lines
    }
    
    func next() -> String? {
        defer { self.lineOffset += 1 }
        return lineOffset < self.lines.count ? self.lines[lineOffset] : nil
    }
    
//    func foreachRow(_ rowParcer: ((String, Int) -> Void)) {
//        for (offset, line) in lines[self.lineOffset...].enumerated() {
//            self.lineOffset += 1
//            rowParcer(line, offset)
//        }
//    }
//    
//    func foreachRow(_ rowParcer: ((String, Int) throws -> Void)) rethrows {
//        for (offset, line) in lines[self.lineOffset...].enumerated() {
//            self.lineOffset += 1
//            try rowParcer(line, offset)
//        }
//    }
//    
//    func foreachRow(skipLines: Int, _ rowParcer: ((String, Int) -> Void)) {
//        self.lineOffset += skipLines
//        for (offset, line) in lines[self.lineOffset...].enumerated() {
//            self.lineOffset += 1
//            rowParcer(line, offset)
//        }
//    }

}
