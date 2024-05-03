//
//  DockingTests.swift
//  
//
//  Created by Cory Kornowicz on 11/18/23.
//
import SwiftVina
import XCTest

final class DockingTests: XCTestCase {
    
    func testBasicDocking() throws {
        
        //Default File Paths for Vina Scoring Function Docking
        
        let receptorPath: String = "/Users/corykornowicz/Downloads/1iep_receptor.pdbqt"
        let ligandPath: String   = "/Users/corykornowicz/Downloads/1iep_ligand.pdbqt"
        let outputPath: String   = "/Users/corykornowicz/Downloads/1iep_ligand_vina_out.pdbqt"
//        let outputMapsPath: String = "/Users/corykornowicz/Downloads/1iep_receptor"
        
        let vinaEngine: Vina = try Vina(sf_name: .SF_VINA, verbosity: .debug, no_refine: false, progress_callback: nil)
        
        try vinaEngine.set_receptor(rigid_filepath: receptorPath)
        print("Receptor Initialized.")

        try vinaEngine.set_ligand_from_filepath(ligand_filepath: ligandPath)
        print("Ligand Initialized.")
        // Will compute maps only for Vina atom types in the ligand(s)
        // In the case users ask for score and local only with the autobox arg, we compute the optimal box size for it/them.
//        if ((score_only || local_only) && autobox) {
//            std::vector<double> dim = v.grid_dimensions_from_ligand(buffer_size);
//            v.compute_vina_maps(dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], grid_spacing, force_even_voxels);
//        } else {
//            v.compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z, grid_spacing, force_even_voxels);
//        }

        let center_x: fl = 15.190
        let center_y: fl = 53.903
        let center_z: fl = 16.917
        let size_x: fl = 20.0
        let size_y: fl = 20.0
        let size_z: fl = 20.0
        
        let gridSpacing: fl = 0.375
        let forceEvenVoxels: Bool = false
        
        try vinaEngine.computeVinaMaps(center_x: center_x, center_y: center_y, center_z: center_z,
                                       size_x: size_x, size_y: size_y, size_z: size_z,
                                       granularity: gridSpacing, force_even_voxels: forceEvenVoxels)
        
//        if (vm.count("write_maps"))
//            try vinaEngine.write_maps(map_prefix: outputMapsPath)
        
        print("Vina Maps Initialized.")
        
        try vinaEngine.global_search(exhaustiveness: 5, n_poses: 9, min_rmsd: 1.0, max_evals: 0)
        
        print("Attempting to write results")
        try vinaEngine.write_poses(output_name: outputPath)
        
        // Score
        let energies = try vinaEngine.score()
        vinaEngine.show_score(energies: energies)
        
        print("Test Complete.")
    }

//    func testPerformanceExample() throws {
//        // This is an example of a performance test case.
//        self.measure {
//            // Put the code you want to measure the time of here.
//        }
//    }

}
