//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

protocol RowIterable<T> {
    associatedtype T
    func next() -> T?
}

extension URL {
    func foreachRow(_ rowParcer:((String, Int) -> Void) )
    {
        //Here we should use path not the absoluteString (wich contains file://)
        let path = self.path
        let m = "r"
        guard let cfilePath = (path as NSString).utf8String else {return}
        
        //Open file with specific mode (just use "r")
        guard let file = fopen(cfilePath, m)
        else {
            print("fopen can't open file: \"\(path)\", mode: \"\(m)\"")
            return
        }
        
        //Row capacity for getline()
        var cap = 0
        
        var row_index = 0
        
        //Row container for getline()
        var cline:UnsafeMutablePointer<CChar>? = nil
        
        //Free memory and close file at the end
        defer{free(cline); fclose(file)}
                    
        while getline(&cline, &cap, file) > 0
        {
            if let crow = cline,
               // the output line may contain '\n' that's why we filtered it
                let s = String(utf8String: crow)?.trimmingCharacters(in: .newlines)
            {
                rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
    
    func foreachRow(_ rowParcer: ((String, Int) throws -> Void)) rethrows {
        //Here we should use path not the absoluteString (wich contains file://)
        let path = self.path
        let m = "r"
        guard let cfilePath = (path as NSString).utf8String else {return}
        
        //Open file with specific mode (just use "r")
        guard let file = fopen(cfilePath, m)
        else {
            print("fopen can't open file: \"\(path)\", mode: \"\(m)\"")
            return
        }
        
        //Row capacity for getline()
        var cap = 0
        
        var row_index = 0
        
        //Row container for getline()
        var cline:UnsafeMutablePointer<CChar>? = nil
        
        //Free memory and close file at the end
        defer{free(cline); fclose(file)}
                    
        while getline(&cline, &cap, file) > 0
        {
            if let crow = cline,
               // the output line may contain '\n' that's why we filtered it
                let s = String(utf8String: crow)?.trimmingCharacters(in: .newlines)
            {
                try rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
    
    func foreachRow(skipLines: Int, _ rowParcer:((String, Int) -> Void) )
    {
        //Here we should use path not the absoluteString (wich contains file://)
        let path = self.path
        let m = "r"
        guard let cfilePath = (path as NSString).utf8String else {return}
        
        //Open file with specific mode (just use "r")
        guard let file = fopen(cfilePath, m)
        else {
            print("fopen can't open file: \"\(path)\", mode: \"\(m)\"")
            return
        }
        
        //Row capacity for getline()
        var cap = 0
        
        var row_index = 0
        
        //Row container for getline()
        var cline:UnsafeMutablePointer<CChar>? = nil
        
        //Free memory and close file at the end
        defer{free(cline); fclose(file)}
        
        // Skip the specified number of lines
        for _ in 0..<skipLines {
            if getline(&cline, &cap, file) <= 0 {
                break // End of file or error encountered
            }
        }
                    
        while getline(&cline, &cap, file) > 0
        {
            if let crow = cline,
               // the output line may contain '\n' that's why we filtered it
//               let s = String(utf8String: crow)?.filter({($0.asciiValue ?? 0) >= 32})
                let s = String(utf8String: crow)?.trimmingCharacters(in: .newlines)
            {
                rowParcer(s, row_index)
            }
            
            row_index += 1
        }
    }
}
