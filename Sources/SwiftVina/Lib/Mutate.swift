//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/16/23.
//

import Foundation

func count_mutable_entites(_ c: Conf) -> sz {
    var counter: Int = 0
    for i in 0..<c.ligands.count {
        counter += 2 + c.ligands[i].torsions.count
    }
    for i in 0..<c.flex.count {
        counter += c.flex[i].torsions.count
    }
    return sz(counter)
}


func mutate_conf(_ c: inout Conf, m: Model, amplitude: fl) {
    let mutatble_entities_num = count_mutable_entites(c)
    if mutatble_entities_num == 0 { return }
    var which: sz = random_sz(0, mutatble_entities_num - 1)
    assert(which < mutatble_entities_num)
    
    for i in 0..<c.ligands.count {
        if(which == 0) {
            c.ligands[i].rigid.position += amplitude * random_inside_sphere()
            return
        }
        which -= 1
        if(which == 0) {
            let gr: fl = m.gyration_radius(sz(i))
            if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
                let rotation: vec = amplitude / gr * random_inside_sphere()
                quaternion_increment(q: &c.ligands[i].rigid.orientation, rotation: rotation)
            }
            return
        }
        which -= 1
        if(which < c.ligands[i].torsions.count) {
            c.ligands[i].torsions[which] = random_fl(-pi, pi)
            return
        }
        which -= sz(c.ligands[i].torsions.count)
    }
    
}
