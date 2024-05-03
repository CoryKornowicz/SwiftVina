//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

class DataLoader {
    
    let filename: String
    let subdir = "Data"
    
    init(filename: String) {
        self.filename = filename
        guard let filePath = Bundle.path(forResource: self.filename, ofType: "dat", inDirectory: subdir) else { return }
        print(filePath)
    }
    
}
