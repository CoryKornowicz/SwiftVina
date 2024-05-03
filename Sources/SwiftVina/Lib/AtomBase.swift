//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/12/23.
//

import Foundation

class AtomBase: AtomType {
    var charge: fl
    
    override init() {
        self.charge = 0.0
        super.init()
    }
    
    override func copy(with zone: NSZone? = nil) -> Any {
        let newAtomBase = AtomBase()
        newAtomBase.charge = self.charge
        newAtomBase.xs = self.xs
        newAtomBase.el = self.el
        newAtomBase.ad = self.ad
        newAtomBase.sy = self.sy
        return newAtomBase
    }
    
}
