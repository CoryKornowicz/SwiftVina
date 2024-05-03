//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/16/23.
//

import Foundation

extension InputStream {

    func readLine(bufferSize: Int = 4096, lineCallback: (String) -> Void) {
        var buffer = [UInt8](repeating: 0, count: bufferSize)
        var currentLineData = Data()

        open()

        while hasBytesAvailable {
            let readCount = read(&buffer, maxLength: bufferSize)
            if readCount < 0 {
                // Handle error
                if let error = streamError {
                    print("Stream read error: \(error.localizedDescription)")
                }
                break
            }

            guard readCount > 0 else {
                continue
            }

            currentLineData.append(buffer, count: readCount)

            while let lineRange = currentLineData.range(of: Data("\n".utf8)) {
                let lineData = currentLineData.subdata(in: 0..<lineRange.lowerBound)
                if let line = String(data: lineData, encoding: .utf8) {
                    lineCallback(line)
                }
                currentLineData.removeSubrange(0..<lineRange.upperBound)
            }
        }

        // Process any remaining data as the last line (if not empty)
        if !currentLineData.isEmpty, let line = String(data: currentLineData, encoding: .utf8) {
            lineCallback(line)
        }

        close()
    }
    
    func readLine(bufferSize: Int = 4096, lineCallback: (String) throws -> Void) rethrows {
        var buffer = [UInt8](repeating: 0, count: bufferSize)
        var currentLineData = Data()

        open()

        while hasBytesAvailable {
            let readCount = read(&buffer, maxLength: bufferSize)
            if readCount < 0 {
                // Handle error
                if let error = streamError {
                    print("Stream read error: \(error.localizedDescription)")
                }
                break
            }

            guard readCount > 0 else {
                continue
            }

            currentLineData.append(buffer, count: readCount)

            while let lineRange = currentLineData.range(of: Data("\n".utf8)) {
                let lineData = currentLineData.subdata(in: 0..<lineRange.lowerBound)
                if let line = String(data: lineData, encoding: .utf8) {
                    try lineCallback(line)
                }
                currentLineData.removeSubrange(0..<lineRange.upperBound)
            }
        }

        // Process any remaining data as the last line (if not empty)
        if !currentLineData.isEmpty, let line = String(data: currentLineData, encoding: .utf8) {
            try lineCallback(line)
        }

        close()
    }
    
}
