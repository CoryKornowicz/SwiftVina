//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

final class QuasiNewtonAux { 
    var m: Model
    var p: precalculate_byatom
    var ig: igrid
    var v: vec
    
    init(m: inout Model, p: precalculate_byatom, ig: igrid, v: vec) {
        self.m = m
        self.p = p
        self.ig = ig
        self.v = v
    }
    
    func operate(_ c: inout Conf, _ g: inout Change) -> fl {
        // Before evaluating conf, we have to update model
        m.set(&c)
        let tmp = m.eval_deriv(p, ig, v, &g)
        return tmp
    }
}

final class QuasiNewton {
    var max_steps: UInt
    var average_required_improvement: fl
    
    init(max_steps: UInt = 1000, average_required_improvement: fl = 0.0) {
        self.max_steps = max_steps
        self.average_required_improvement = average_required_improvement
    }
    
    func operate(_ m: inout Model, _ p: precalculate_byatom, _ ig: igrid, _ out: inout OutputType, _ g: inout Change, _ v: vec, _ evalcount: inout Int) {
        var aux = QuasiNewtonAux(m: &m, p: p, ig: ig, v: v)
        let res = bfgs(&aux, x: &out.c, g: &g, max_steps: max_steps, average_required_improvement: average_required_improvement, over: 10, evalcount: &evalcount)
        // Update model a last time after optimization
        m.set(&out.c)
        out.e = res
    }
}
