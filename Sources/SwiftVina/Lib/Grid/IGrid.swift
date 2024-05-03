//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

protocol igrid { // grids interface (that cache, etc. conform to)
    func eval(m: Model, v: fl) -> fl       // needs m.coords // clean up
    func eval_intra(m: Model, v: fl) -> fl // only flexres-grids
    func eval_deriv(m: Model, v: fl) -> fl // needs m.coords, sets m.minus_forces // clean up
}
