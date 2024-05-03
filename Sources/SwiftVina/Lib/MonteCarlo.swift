//
//  MonteCarlo.swift
//  
//
//  Created by Cory Kornowicz on 11/13/23.
//

import Foundation

@inlinable
func metropolis_accept(_ old_fl: fl, _ new_f: fl, temperature: fl) -> Bool {
    if new_f < old_fl { return true }
    let acceptance_probabiliy: fl = exp((old_fl - new_f) / temperature)
    return random_fl(0, 1) < acceptance_probabiliy
}

final class MonteCarlo {
    
    var max_evals: UInt = 0
    var global_steps: UInt = 2500
    var temperature: fl = 1.2
    var hunt_cap: vec = vec(10, 1.5, 10)
    var min_rmsd: fl = 0.5
    var num_saved_mins: UInt = 50
    var mutation_amplitude: fl = 2
    var local_steps: UInt = 0
    
    var debugLevel: VerbosityLevel
    
    init(debugLevel: VerbosityLevel = .critical) {
        self.debugLevel = debugLevel
    }
    // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  global_steps = 50*lig_atoms = 2500
    
    func operate<T: Incrementable>(m: inout Model, p: precalculate_byatom, ig: igrid, corner1: vec, corner2: vec, increment_me: T?) -> OutputType {
        var tmp: output_container = []
        self.operate(m: &m, out: &tmp, p: p, ig: ig, corner1: corner1, corner2: corner2, increment_me: increment_me)
        assert(!tmp.isEmpty)
        return tmp.first!
    }
    
    // out is sorted
    func operate<T: Incrementable>(m: inout Model, out: inout output_container, p: precalculate_byatom, ig: igrid, corner1: vec, corner2: vec, increment_me: T?) {
        var evalcount: Int = 0
        let authentic_v: vec = vec(1000,1000,1000) // FIXME? this is here to avoid max_fl/max_fl
        let s: ConfSize = m.get_size()
        var g: Change = Change(s)
        var tmp: OutputType = OutputType(Conf(s), 0)
        tmp.c.randomize(corner1, corner2)
        
        var best_e: fl = fl.greatestFiniteMagnitude
        assert(local_steps > 0)
        let quasi_newton_par: QuasiNewton = QuasiNewton(max_steps: local_steps)
        
        for step in 0..<global_steps {
            if let increment = increment_me { increment.increment() } // increment
            
            if((max_evals > 0) && (evalcount > max_evals)) { break }
            
            var candidate: OutputType = tmp
            mutate_conf(&candidate.c, m: m, amplitude: mutation_amplitude)
            
            quasi_newton_par.operate(&m, p, ig, &candidate, &g, hunt_cap, &evalcount)
            
            if(step == 0 || metropolis_accept(tmp.e, candidate.e, temperature: temperature)) {
                tmp = candidate
                m.set(&tmp.c) // FIXME? useless?
                // FIXME only for very promising ones
                if(tmp.e < best_e || out.count < num_saved_mins) {
                    quasi_newton_par.operate(&m, p, ig, &tmp, &g, authentic_v, &evalcount)
                    m.set(&tmp.c) // FIXME? useless?
                    tmp.coords = m.get_heavy_atom_movable_coords()
                    add_to_output_container(&out, tmp, min_rmsd, sz(num_saved_mins)) // 20 - max size
                    if(tmp.e < best_e) {
                        best_e = tmp.e
                        if self.debugLevel <= .log {
                            print("New Best E: \(best_e) Step: \(step)")
                        }
                    }
                }
            }
        }
        assert(!out.isEmpty)
        assert(out.first!.e <= out.last!.e)
    }
}

typealias parallel_mc_task_container = ContiguousArray<ParallelMCTask>

final class ParallelMCTask {
    var m: Model
    var out: output_container
    
    init(m: Model) {
        self.m = m.copy() as! Model
        self.out = output_container()
    }
}

final class ParallelMCAux {
    let mc: MonteCarlo
    let p: precalculate_byatom
    let ig: igrid
    let corner1: vec
    let corner2: vec
    let pg: ParallelProgress?
    
    init(mc: inout MonteCarlo, p: precalculate_byatom, ig: igrid, corner1: vec, corner2: vec, pg: ParallelProgress?) {
        self.mc = mc
        self.p = p
        self.ig = ig
        self.corner1 = corner1
        self.corner2 = corner2
        self.pg = pg
    }
    
    func operate(_ t: ParallelMCTask) {
        mc.operate(m: &t.m, out: &t.out, p: p, ig: ig, corner1: corner1, corner2: corner2, increment_me: pg)
    }
    
}

fileprivate func merge_output_containers(in_: output_container, out: inout output_container, min_rmsd: fl, max_size: sz) {
    for i in 0..<in_.count {
        add_to_output_container(&out, in_[i], min_rmsd, max_size)
    }
}

fileprivate func merge_output_containers(many: parallel_mc_task_container, out: inout output_container, min_rmsd: fl = 2.0, max_size: sz) {
    // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
    for i in 0..<many.count {
        merge_output_containers(in_: many[i].out, out: &out, min_rmsd: min_rmsd, max_size: max_size)
    }
    out.sort()
}

final class ParallelMC {
    var mc: MonteCarlo = MonteCarlo()
    var num_tasks: UInt = 8
    var num_threads: UInt = 1
    var display_progress: Bool = true
        
    func operate(m: Model, out: inout output_container, p: precalculate_byatom, ig: igrid,
                 corner1: vec, corner2: vec, progress_callback: ((fl) -> Void)?) {
        
        let pp: ParallelProgress = ParallelProgress(callback: progress_callback, count: num_tasks * mc.global_steps)
        let parallel_mc_aux_instance = ParallelMCAux(mc: &mc, p: p, ig: ig, corner1: corner1, corner2: corner2, pg: display_progress ? pp : nil)
        var task_container: parallel_mc_task_container = parallel_mc_task_container()
        
        for _ in 0..<num_tasks {
            task_container.append(ParallelMCTask(m: m))
        }
        
        DispatchQueue.concurrentPerform(iterations: Int(num_tasks)) { i in
            parallel_mc_aux_instance.operate(task_container[i])
        }

        merge_output_containers(many: task_container, out: &out, min_rmsd: mc.min_rmsd, max_size: sz(mc.num_saved_mins))
    }
    
}
