//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/16/23.
//

import Foundation

protocol Incrementable {
    func increment()
}

final class ProgressBar {
    private var total: UInt
    private var current: UInt
    
    init(total: UInt) {
        self.total = total
        self.current = 0
    }
    
    func update() -> UInt {
        self.current += 1
        printProgressBar()
        return self.current
    }
    
    private func printProgressBar() {
        let percentage = (current * 100) / total
        let maxSymbols: UInt = 50 // Number of symbols for 100% progress
        let currentSymbols: Int = Int((percentage * maxSymbols) / 100)
        let remainingSymbols: Int = Int(maxSymbols) - currentSymbols
        
        let progressBar = String(repeating: "*", count: currentSymbols) +
                          String(repeating: "-", count: remainingSymbols)
        let progressLine = "0% |\(progressBar)| 100%"
        
        print("\u{1B}[2K\r\(progressLine) \(percentage)%", terminator: "")
        fflush(stdout)
    }
}

final class ParallelProgress: Incrementable {
    
    private let p: ProgressBar
    private let callback: ((fl) -> Void)?
    private let count: UInt
    
    private let messageQueue: DispatchQueue = DispatchQueue(label: "parallelPrintStreamQueue", qos: .utility)

    init(callback: ((fl) -> Void)?, count: UInt) {
        self.callback = callback
        self.count = count
        self.p = ProgressBar(total: count)
    }
    
    func increment() {
        self.messageQueue.async {
            // print asynchronously
            let printCount = self.p.update()
            if let callback = self.callback {
                callback(fl( printCount / self.count))
            }
        }
    }
}



