//
//  TextStreamer.swift
//  
//
//  Created by Cory Kornowicz on 11/16/23.
//

import Foundation

struct TextStreamer: TextOutputStream {
    
    var filepath: String
    
    lazy var fileHandle: FileHandle? =  {
        if !FileManager.default.fileExists(atPath: filepath) {
            FileManager.default.createFile(atPath: filepath, contents: nil)
        }
        let fileHandle = FileHandle(forWritingAtPath: filepath)
        return fileHandle
    }()
    
    mutating func write(_ string: String) {
        fileHandle?.seekToEndOfFile()
        fileHandle?.write(string.data(using:.utf8)!)
    }
}
