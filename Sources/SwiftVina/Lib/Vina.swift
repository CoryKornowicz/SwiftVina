//
//  File.swift
//  
//
//  Created by Cory Kornowicz on 11/11/23.
//

import Foundation
import simd

public enum VerbosityLevel: Int, Comparable {
    
    case debug = 0
    case log = 1
    case critical = 2
    case silent = 3
    
    public static func < (lhs: VerbosityLevel, rhs: VerbosityLevel) -> Bool {
        lhs.rawValue < rhs.rawValue
    }
    
}

public final class Vina {
    
    // model and poses
    var m_receptor: Model?
    var m_model: Model?
    var m_poses: output_container = []
    var m_receptor_initialized: Bool = false
    var m_ligand_initialized: Bool = false
    //scoring function
    var m_sf_choice: SFChoice
    var m_weights: flv = []
    var m_scoring_function: ScoringFunction?
    var m_precalculated_byatom: precalculate_byatom?
    var m_precalculated_sf: precalculate?
    // maps
    var m_grid: Cache?
    // var m_ad4grid: ad4cache
    var m_non_cache: Noncache?
    var m_map_initialized: Bool = false
    // global search
    //others
    var m_verbosity: VerbosityLevel = .critical
    var m_no_refine: Bool = false
    var m_progress_callback: ((fl) -> Void)? = nil
    
    public init(sf_name: SFChoice = .SF_VINA, verbosity: VerbosityLevel = .critical, no_refine: Bool = false, progress_callback: ((fl) -> Void)? = nil) throws {
        self.m_verbosity = verbosity
        self.m_receptor_initialized = false
        self.m_ligand_initialized = false
        self.m_map_initialized = false
        self.m_no_refine = no_refine
        self.m_progress_callback = progress_callback
            
        switch sf_name {
        case .SF_VINA:
            self.m_sf_choice = .SF_VINA
            self.setVinaWeights()
        case .SF_AD42:
            self.m_sf_choice = .SF_VINARDO
            self.setVinardoWeights()
        case .SF_VINARDO:
            self.m_sf_choice = .SF_AD42
            self.setAd4Weights()
        }
    }
    
    func setVinaWeights(weight_gauss1: fl = -0.035579, weight_gauss2: fl = -0.005156,
                          weight_repulsion: fl = 0.840245, weight_hydrophobic: fl = -0.035069,
                          weight_hydrogen: fl = -0.587439, weight_glue: fl = 50,
                          weight_rot: fl = 0.05846) {
        
        var weights: flv = []

        if (m_sf_choice == .SF_VINA) {
            weights.append(weight_gauss1)
            weights.append(weight_gauss2)
            weights.append(weight_repulsion)
            weights.append(weight_hydrophobic)
            weights.append(weight_hydrogen)
            weights.append(weight_glue)
            weights.append(5 * weight_rot / 0.1 - 1)
            
            // Store in Vina object
            m_weights = weights
            
            // Since we set (different) weights, we automatically initialize the forcefield
            setForceField()
        }
    }
    
    func setVinardoWeights(weight_gauss1: fl = -0.045,
                           weight_repulsion: fl = 0.8, weight_hydrophobic: fl = -0.035,
                           weight_hydrogen: fl = -0.600, weight_glue: fl = 50,
                           weight_rot: fl = 0.05846) {
        

        var weights: flv = []

        if (m_sf_choice == .SF_VINARDO) {
            weights.append(weight_gauss1)
            weights.append(weight_repulsion)
            weights.append(weight_hydrophobic)
            weights.append(weight_hydrogen)
            weights.append(weight_glue)
            weights.append(5 * weight_rot / 0.1 - 1)
            
            // Store in Vina object
            m_weights = weights
            
            // Since we set (different) weights, we automatically initialize the forcefield
            setForceField()
        }
        
    }
    
    
    func setAd4Weights(weight_ad4_vdw: fl = 0.1662, weight_ad4_hb: fl = 0.1209,
                       weight_ad4_elec: fl = 0.1406, weight_ad4_dsolv: fl = 0.1322,
                       weight_glue: fl = 50, weight_ad4_rot: fl = 0.2983) {
        
        var weights: flv = []

        if (m_sf_choice == .SF_AD42) {
            weights.append(weight_ad4_vdw)
            weights.append(weight_ad4_hb)
            weights.append(weight_ad4_elec)
            weights.append(weight_ad4_dsolv)
            weights.append(weight_glue)
            weights.append(weight_ad4_rot)
            
            // Store in Vina object
            m_weights = weights
            
            // Since we set (different) weights, we automatically initialize the forcefield
            setForceField()
        }
        
    }

    private func setForceField() {
        let scoringFunction = ScoringFunction(scoring_function_choice: self.m_sf_choice, weights: self.m_weights)
        self.m_scoring_function = scoringFunction
    }
        
    private func vina_remarks(pose: OutputType, lb: fl, ub: fl) -> String {
        var remark: String = ""

        remark += "REMARK VINA RESULT: "
        remark += String(format: "%9.3f", pose.e)
        remark += "  "
        remark += String(format: "%9.3f", lb)
        remark += "  "
        remark += String(format: "%9.3f", ub)
        remark += "\n"

        remark += "REMARK INTER + INTRA:    "
        remark += String(format: "%12.3f", pose.total)
        remark += "\n"
        remark += "REMARK INTER:            "
        remark += String(format: "%12.3f", pose.inter)
        remark += "\n"
        remark += "REMARK INTRA:            "
        remark += String(format: "%12.3f", pose.intra)
        remark += "\n"
        if (m_sf_choice == .SF_AD42) {
            remark += "REMARK CONF_INDEPENDENT: "
            remark += String(format: "%12.3f", pose.conf_independent)
            remark += "\n"
        }
        remark += "REMARK UNBOUND:          "
        remark += String(format: "%12.3f", pose.unbound)
        remark += "\n"

        return remark
    }
        
    func cite() { 
        var cite_message = ""
        cite_message += "#################################################################\n"
        cite_message += "# If you used AutoDock Vina in your work, please cite:          #\n"
        cite_message += "#                                                               #\n"
        cite_message += "# J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli  #\n"
        cite_message += "# AutoDock Vina 1.2.0: New Docking Methods, Expanded Force      #\n"
        cite_message += "# Field, and Python Bindings, J. Chem. Inf. Model. (2021)       #\n"
        cite_message += "# DOI 10.1021/acs.jcim.1c00203                                  #\n"
        cite_message += "#                                                               #\n"
        cite_message += "# O. Trott, A. J. Olson,                                        #\n"
        cite_message += "# AutoDock Vina: improving the speed and accuracy of docking    #\n"
        cite_message += "# with a new scoring function, efficient optimization and       #\n"
        cite_message += "# multithreading, J. Comp. Chem. (2010)                         #\n"
        cite_message += "# DOI 10.1002/jcc.21334                                         #\n"
        cite_message += "#                                                               #\n"
        cite_message += "# Please see https://github.com/ccsb-scripps/AutoDock-Vina for  #\n"
        cite_message += "# more information.                                             #\n"
        cite_message += "#################################################################\n"

        print(cite_message)
    }

    public func set_receptor(rigid_filepath: String, flex_filepath: String? = nil) throws {
        // Read the receptor PDBQT file
        /* CONDITIONS:
            - 1. AD4/Vina rigid  NO, flex  NO: FAIL
            - 2. AD4      rigid YES, flex YES: FAIL
            - 3. AD4      rigid YES, flex  NO: FAIL
            - 4. AD4      rigid  NO, flex YES: SUCCESS (need to read maps later)
            - 5. Vina     rigid YES, flex YES: SUCCESS
            - 6. Vina     rigid YES, flex  NO: SUCCESS
            - 7. Vina     rigid  NO, flex YES: SUCCESS (need to read maps later)
        */
        if rigid_filepath.isEmpty && (flex_filepath != nil && !flex_filepath!.isEmpty) && m_sf_choice == .SF_VINA {
            // CONDITION 1
            throw VinaError("ERROR: No (rigid) receptor or flexible residues were specified. (vina.cpp)\n")
        } else if m_sf_choice == .SF_AD42 && !rigid_filepath.isEmpty {
            // CONDITIONS 2, 3
            throw VinaError("ERROR: Only flexible residues allowed with the AD4 scoring function. No (rigid) receptor.\n")
        }

        // CONDITIONS 4, 5, 6, 7 (rigid_name and flex_name are empty strings per default)
        m_receptor = try parse_receptor_pdbqt(rigid: rigid_filepath, flex: flex_filepath, 
                                              atype: m_scoring_function!.get_atom_typing())
        m_model = m_receptor
        m_receptor_initialized = true
        // If we are reading another receptor we should not consider the ligand and the map as initialized anymore
        m_ligand_initialized = false
        m_map_initialized = false
    }
    
    public func set_ligand_from_filepath(ligand_filepath: String) throws {
        // Read ligand PDBQT string and add it to the model
        if ligand_filepath.isEmpty {
            throw VinaError("ERROR: Cannot read ligand file. Ligand string is empty.\n")
        }

        guard let atom_typing = m_scoring_function?.get_atom_typing() else {
            throw VinaError("ERROR: Cannot get atom typing.\n")
        }

        if (!m_receptor_initialized) {
            // This situation will happen if we don't need a receptor and we are using affinity maps
            let m = Model(atom_typing)
            m_model = m
            m_receptor = m.copy() as? Model
        } else {
            // Replace current model with receptor and reinitialize poses
            m_model = m_receptor!.copy() as? Model
        }

        // ... and add ligand to the model
        
        guard m_model != nil else { throw VinaError("ERROR: Model is nil.\n")}
        
        m_model!.append(try parse_ligand_pdbqt_from_file(filePath: ligand_filepath, atype: atom_typing))

        // Because we precalculate ligand atoms interactions
        let precalculated_byatom: precalculate_byatom = precalculate_byatom(sf: m_scoring_function!, model: m_model!)

        // Check that all atom types are in the grid (if initialized)
        if (m_map_initialized) {
            let atom_types: szv = m_model!.get_movable_atom_types(atom_typing)

            if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
                guard let cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            } else {
                guard let ad4_cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!ad4_cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            }
        }

        // Store in Vina object
        m_poses = output_container()
        m_precalculated_byatom = precalculated_byatom
        m_ligand_initialized = true
    }
    
    func set_ligand_from_string(ligand_string: String) throws {
        // Read ligand PDBQT string and add it to the model
        if ligand_string.isEmpty {
            throw VinaError("ERROR: Cannot read ligand file. Ligand string is empty.\n")
        }

        guard let atom_typing = m_scoring_function?.get_atom_typing() else {
            throw VinaError("ERROR: Cannot get atom typing.\n")
        }

        if (!m_receptor_initialized) {
            // This situation will happen if we don't need a receptor and we are using affinity maps
            let m = Model(atom_typing)
            m_model = m
            m_receptor = m.copy() as? Model
        } else {
            // Replace current model with receptor and reinitialize poses
            m_model = m_receptor!.copy() as? Model
        }

        // ... and add ligand to the model
        
        guard m_model != nil else { throw VinaError("ERROR: Model is nil.\n")}
        
        m_model!.append(try parse_ligand_pdbqt_from_string(stringContents: ligand_string, atype: atom_typing))

        // Because we precalculate ligand atoms interactions
        let precalculated_byatom: precalculate_byatom = precalculate_byatom(sf: m_scoring_function!, model: m_model!)

        // Check that all atom types are in the grid (if initialized)
        if (m_map_initialized) {
            let atom_types: szv = m_model!.get_movable_atom_types(atom_typing)

            if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
                guard let cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            } else {
                guard let ad4_cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!ad4_cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            }
        }

        // Store in Vina object
        m_poses = output_container()
        m_precalculated_byatom = precalculated_byatom;
        m_ligand_initialized = true;
    }
    
    func set_ligand_from_string(ligand_strings: [String]) throws {
        if ligand_strings.isEmpty {
            throw VinaError("ERROR: Cannot read ligand files. Ligand files are empty.\n")
        }

        guard let atom_typing = m_scoring_function?.get_atom_typing() else {
            throw VinaError("ERROR: Cannot get atom typing.\n")
        }

        if (!m_receptor_initialized) {
            // This situation will happen if we don't need a receptor and we are using affinity maps
            let m = Model(atom_typing)
            m_model = m
            m_receptor = m.copy() as? Model
        } else {
            // Replace current model with receptor and reinitialize poses
            m_model = m_receptor!.copy() as? Model
        }
        
        for ligand_string in ligand_strings {
            m_model!.append(try parse_ligand_pdbqt_from_string(stringContents: ligand_string, atype: atom_typing))
        }
        
        // Because we precalculate ligand atoms interactions
        let precalculated_byatom: precalculate_byatom = precalculate_byatom(sf: m_scoring_function!, model: m_model!)

        // Check that all atom types are in the grid (if initialized)
        if (m_map_initialized) {
            let atom_types: szv = m_model!.get_movable_atom_types(atom_typing)

            if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
                guard let cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            } else {
                guard let ad4_cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!ad4_cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            }
        }

        // Store in Vina object
        m_poses = output_container()
        m_precalculated_byatom = precalculated_byatom;
        m_ligand_initialized = true;
    }

    
    func set_ligand_from_file(ligand_file_path: String) throws {
        let urlContentsString = try NSString(contentsOfFile: ligand_file_path, encoding: NSUTF8StringEncoding) as String
        try set_ligand_from_string(ligand_string: urlContentsString)
    }

    func set_ligand_from_file(ligand_file_paths: [String]) throws {
        var parsableStrings: [String] = [String]()
        for ligand_file_path in ligand_file_paths {
            let urlContentsString = try NSString(contentsOfFile: ligand_file_path, encoding: NSUTF8StringEncoding) as String
            parsableStrings.append(urlContentsString)
        }
        try set_ligand_from_string(ligand_strings: parsableStrings)
    }
    
    func grid_dimensions_from_ligand(buffer_size: fl) -> [fl] {
        var box_dimensions: [fl] = [fl](repeating: 0, count: 6)
        var box_center: [fl] = [fl](repeating: 0, count: 3)
        var max_distance: [fl] = [fl](repeating: 0, count: 3)

        // The center of the ligand will be the center of the box
        box_center = m_model!.center()

        // Get the furthest atom coordinates from the center in each dimensions
        for i in 0..<m_model!.num_movable_atoms {
            let atom_coords = m_model!.get_coords(i)

            for j in 0..<3 {
                let distance = abs(box_center[j] - atom_coords[j])

                if (max_distance[j] < distance) {
                    max_distance[j] = distance
                }
            }
        }

        // Get the final dimensions of the box
        box_dimensions[0] = box_center[0]
        box_dimensions[1] = box_center[1]
        box_dimensions[2] = box_center[2]
        box_dimensions[3] = ceil((max_distance[0] + buffer_size) * 2)
        box_dimensions[4] = ceil((max_distance[1] + buffer_size) * 2)
        box_dimensions[5] = ceil((max_distance[2] + buffer_size) * 2)

        return box_dimensions
    }
    
    public func computeVinaMaps(center_x: fl, center_y: fl, center_z: fl,
                         size_x: fl, size_y: fl, size_z: fl,
                         granularity: fl = 0.5, force_even_voxels: Bool = false) throws {
        // Setup the search box
        // Check first that the receptor was added
        if m_sf_choice == .SF_AD42 {
            throw VinaError("Cannot compute Vina affinity maps using the AD4 scoring function.")
        } else if !m_receptor_initialized {
            throw VinaError("Cannot compute Vina or Vinardo affinity maps. The (rigid) receptor was not initialized.")
        } else if (size_x <= 0 || size_y <= 0 || size_z <= 0) {
            throw VinaError("Grid box dimensions must be greater than 0 Angstrom")
        } else if (size_x * size_y * size_z > 27e3) {
            if m_verbosity <= .critical {
                print("WARNING: Search space volume is greater than 27000 Angstrom^3 (See FAQ)")
            }
        }
        
        var gd: grid_dims = grid_dims(repeating: grid_dim(), count: 3)
        let span: vec = vec(size_x, size_y, size_z)
        let center: vec = vec(center_x, center_y, center_z)
        let slope: fl = 1e6 // FIXME: too large? used to be 100
        var atom_types: szv
        let atom_typing: AtomType.T = m_scoring_function!.get_atom_typing()
        
        /* Atom types initialization
         If a ligand was defined before, we only use those present in the ligand
         otherwise we use all the atom types present in the forcefield
         */
        
        if (m_ligand_initialized) {
            atom_types = m_model!.get_movable_atom_types(atom_typing)
        } else {
            atom_types = m_scoring_function!.get_atom_types()
        }
        
        // Grid dimensions
        
        for i in 0..<3 { // grid_dims will always have 3 dimensions
            gd[i].n_voxels = sz(ceil(span[i] / granularity))
            
            // If odd n_voxels increment by 1
            if (force_even_voxels && (gd[i].n_voxels % 2 == 1)) {
                // because sample points (npts) == n_voxels + 1
                gd[i].n_voxels += 1
            }
            
            let real_span = granularity * fl(gd[i].n_voxels)
            gd[i].begin = center[i] - real_span / 2
            gd[i].end = gd[i].begin + real_span
        }
        
        // Initialize the scoring function
        let precalculated_sf: precalculate = precalculate(sf: m_scoring_function!)
        // Store it now in Vina object because of non_cache
        m_precalculated_sf = precalculated_sf
        
        if (m_sf_choice == .SF_VINA) {
            if m_verbosity <= .debug {
                print("Computing Vina grid")
            }
        } else {
            if m_verbosity <= .debug {
                print("Computing Vinardo grid")
            }
        }
        
        // Compute the Vina grids
        
        m_grid = Cache(gd, slope)
        m_grid!.populate(m: m_model!, p: precalculated_sf, atom_types_needed: atom_types)
        
        if m_verbosity <= .debug {
            print("Grid initialized")
        }
        
        // create non_cache for scoring with explicit receptor atoms (instead of grids)
        if (!m_no_refine) {
            m_non_cache = Noncache(m: m_model!, gd_: gd, p_: m_precalculated_sf!, slope_: slope)
        }
        
        // Store in Vina object
        m_map_initialized = true
    }

    func load_maps(maps: String) throws { 
        let slope: fl = 1e6 // FIXME: too large? used to be 100

        if m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO { 
            print("Reading Vina maps")
            let grid: Cache = Cache(slope)
            try grid.read(maps)
            print("Done reading Vina maps")
            m_grid = grid
        } else { 
            fatalError("Not Implemented Yet")
        }

        // Check that all the affinity map are present for ligands/flex residues (if initialized already)
        if m_ligand_initialized {
            let atom_typing: AtomType.T = m_scoring_function!.get_atom_typing()
            let atom_types: szv = m_model!.get_movable_atom_types(atom_typing)

            if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
                guard let cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            } else {
                guard let ad4_cache = m_grid else {
                    throw VinaError("Cache is not initialized")
                }
                if(!ad4_cache.are_atom_types_grid_initialized(atom_types)) {
                    throw VinaError("Atom Type Grids are not initialized")
                }
            }
        }

        // Store in Vina object
        m_map_initialized = true;
    }

    public func randomize(max_steps: Int = 10000) throws {
        // Randomize ligand/flex residues conformation
        // Check the box was defined
        if !m_ligand_initialized {
            throw VinaError("Cannot do ligand randomization. Ligand(s) was(ere) not initialized.")
        } else if !m_map_initialized {
            throw VinaError("Cannot do ligand randomization. Affinity maps were not initialized.")
        }

        var c: Conf 
        var penalty: fl = 0
        var best_clash_penalty: fl = 0

        // It's okay to take the initial conf since we will randomize it
        let init_conf: Conf = m_model!.get_initial_conf()
        var best_conf: Conf = init_conf

        if m_verbosity <= .log {
            print("Randomize conformation...")
        }

        for i in 0..<max_steps { 
            c = init_conf
            c.randomize(m_grid!.corner1(), m_grid!.corner2())
            m_model!.set(&c)
            penalty = m_model!.clash_penalty()
            if m_verbosity <= .debug {
                print("Clash Penalty: \(penalty)")
            }
            if (i == 0 || penalty < best_clash_penalty) {
                best_conf = c
                best_clash_penalty = penalty
            }
        }
        
        if m_verbosity <= .log {
            print("Done.")
        }

        m_model!.set(&best_conf)

        if m_verbosity <= .log {
            print("Clash penalty: \(best_clash_penalty)")
        }
    }

    func remove_redundant(in_: output_container, min_rmsd: fl) -> output_container {
        var tmp: output_container = output_container() 
        for i in 0..<in_.count { 
            add_to_output_container(&tmp, in_[i], min_rmsd, sz(in_.count))
        }
        return tmp
    }
    
    public func score() throws -> flv { 
        // Score the current conf in the model
        // Check if ff and ligand were initialized
        // Check if the ligand is not outside the box

        if !m_ligand_initialized {
            throw VinaError("Cannot score the pose. Ligand(s) was(ere) not initialized.")
        } else if !m_map_initialized {
            throw VinaError("Cannot score the pose. Affinity maps were not initialized.")
        } else if !m_grid!.is_in_grid(m_model!) {
            throw VinaError("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.")
        }

        var intramolecular_energy: fl = 0
        let authentic_v: vec = vec(1000, 1000, 1000)

        if m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO {
            intramolecular_energy = m_model!.eval_intramolecular(m_precalculated_byatom!, ig: m_grid!, v: authentic_v)
        }

        let energies = score(intramolecular_energy: intramolecular_energy)
        return energies
    }

    func score(intramolecular_energy: fl) -> flv { 
        // Score the current conf in the model

        var total: fl = 0
        var inter: fl = 0
        var intra: fl = 0
        var all_grids: fl = 0 // ligand & flex
        var lig_grids: fl = 0
        var flex_grids: fl = 0
        var lig_intra: fl = 0
        var conf_independent: fl = 0
        var inter_pairs: fl = 0
        var intra_pairs: fl = 0
        let authentic_v: vec = vec(1000, 1000, 1000)
        var energies: flv = []
        
        if m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO { 

            // Inter
            if m_no_refine || !m_receptor_initialized { 
                all_grids = m_grid!.eval(m: m_model!, v: authentic_v[1]) // [1] ligand & flex -- grid
            } else { 
                all_grids = m_non_cache!.eval(m: m_model!, v: authentic_v[1]) // [1] ligand & flex -- grid
            }
            inter_pairs = m_model!.eval_inter(m_precalculated_byatom!, authentic_v) // [1] ligand -- flex
            // Intra
            if m_no_refine || !m_receptor_initialized {
                flex_grids = m_grid!.eval_intra(m: m_model!, v: authentic_v[1]) // [1] flex -- grid
            } else { 
                flex_grids = m_non_cache!.eval_intra(m: m_model!, v: authentic_v[1]) // [1] flex -- grid
            }
            intra_pairs = m_model!.evalo(m_precalculated_byatom!, authentic_v) // [1] flex_i -- flex_i and flex_i -- flex_j
            lig_grids = all_grids - flex_grids
            inter = lig_grids + inter_pairs
            lig_intra = m_model!.evali(m_precalculated_byatom!, v: authentic_v) // [2] ligand_i -- ligand_i
            intra = flex_grids + intra_pairs + lig_intra
            // Total
            total = m_scoring_function!.conf_independent(m_model!, inter + intra - intramolecular_energy) // we pass intermolecular energy from the best pose
            // Torsion, we want to know how much torsion penalty was added to the total energy
            conf_independent = total - (inter + intra - intramolecular_energy)
        } else { 
            fatalError("Not Implemented Yet")
        //     // Inter
        //     lig_grids = m_ad4grid.eval(m_model, authentic_v[1]) // [1] ligand -- grid
        //     inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // [1] ligand -- flex
        //     inter = lig_grids + inter_pairs;
        //     // Intra
        //     flex_grids = m_ad4grid.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
        //     intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1] flex_i -- flex_i and flex_i -- flex_j
        //     lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2] ligand_i -- ligand_i
        //     intra = flex_grids + intra_pairs + lig_intra;
        //     // Torsion
        //     conf_independent = m_scoring_function.conf_independent(m_model, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
        //     // Total
        //     total = inter + conf_independent; // (+ intra - intra)
        }

        energies.append(total)
        energies.append(lig_grids)
        energies.append(inter_pairs)
        energies.append(flex_grids)
        energies.append(intra_pairs)
        energies.append(lig_intra)
        energies.append(conf_independent)    

        if (m_sf_choice == .SF_VINA  || m_sf_choice == .SF_VINARDO) {
            energies.append(intramolecular_energy)
        } else {
            energies.append(intra)
        }

        return energies
    }

    public func optimize(max_steps: UInt = 0) throws -> flv {
        // Local optimization of the ligand conf
        // Check if ff, box and ligand were initialized
        // Check if the ligand is not outside the box

        if !m_ligand_initialized {
            throw VinaError("Cannot do the optimization. Ligand(s) was(ere) not initialized.")
        } else if !m_map_initialized {
            throw VinaError("Cannot do the optimization. Affinity maps were not initialized.")
        } else if !m_grid!.is_in_grid(m_model!) {
            throw VinaError("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.")
        }

        var e: fl = 0
        var c: Conf

        if !m_poses.isEmpty { 
            // if m_poses is not empty, it means that we did a docking before
            // But it is really that useful to minimize after docking?
            e = m_poses[0].e
            c = m_poses[0].c
        } else { 
            c = m_model!.get_initial_conf()
        }

        var out: OutputType = OutputType(c, e)

        // std::vector<double> energies = optimize(out, max_steps);
        let energies = try optimize(&out, max_steps: max_steps)

        return energies
    }

    private func optimize(_ out: inout OutputType, max_steps: UInt = 0) throws -> ContiguousArray<fl> {
        // Local optimization of the ligand conf
        var g: Change = Change(m_model!.get_size())
        let quasi_newton_par: QuasiNewton = QuasiNewton(max_steps: max_steps)
        let authentic_v: vec = vec(1000, 1000, 1000)
        var energies_before_opt: ContiguousArray<fl> = []
        var energies_after_opt: ContiguousArray<fl> = []
        var evalcount: Int = 0

        // Define the number minimization steps based on the number moving atoms

        if (max_steps == 0) {
            let new_max_steps = ((25 + m_model!.num_movable_atoms) / 3)
            if (m_verbosity <= .debug) {
                print("Number of local optimization steps: \(new_max_steps)")
            }
            quasi_newton_par.max_steps = UInt(new_max_steps)
        } else { 
            quasi_newton_par.max_steps = max_steps
        }

        if m_verbosity <= .debug { 
            print("Before local optimization: ")
            energies_before_opt = try score()
            show_score(energies: energies_before_opt)
        }

        // doing("Performing local search", m_verbosity, 0);
        if m_verbosity <= .log { 
            print("Performing local search...")
        }
        
        // Try 5 five times to optimize locally the conformation
        for _ in 0..<5 { 
            if m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO { 
                quasi_newton_par.operate(&m_model!, m_precalculated_byatom!, m_grid!, &out, &g, authentic_v, &evalcount)
                // Break if we succeed to bring (back) the ligand within the grid
                if (m_grid!.is_in_grid(m_model!)) {
                    break
                }
            } else { 
                fatalError("Not Implemented")
            }
        }

        if m_verbosity <= .log { 
            print("Done with optimization.")
        }

        energies_after_opt = try score()

        return energies_after_opt
    }

    public func global_search(exhaustiveness: Int = 8, n_poses: Int = 20, min_rmsd: fl = 1.0, max_evals: Int = 0) throws {
        // Vina search (Monte-carlo and local optimization)
        // Check if ff, box and ligand were initialized

        if (!m_ligand_initialized) {
            throw VinaError("Cannot do the global search. Ligand(s) was(ere) not initialized.")
        } else if (!m_map_initialized) {
            throw VinaError("Cannot do the global search. Affinity maps were not initialized.")
        } else if (exhaustiveness < 1) {
            throw VinaError("Exhaustiveness must be 1 or greater")
        }

        if (exhaustiveness < 8) {
            if m_verbosity <= .critical {
                print("WARNING: At low exhaustiveness, parallization will be under utilized.")
            }
        }

        var intramolecular_energy: fl = 0
        let authentic_v: vec = vec(1000, 1000, 1000)
        var poses: output_container = output_container()

        // Setup Monte-Carlo search
        let parallelmc: ParallelMC = ParallelMC()
        let heuristic: sz = m_model!.num_movable_atoms + 10 * m_model!.get_size().num_degrees_of_freedom()
        parallelmc.mc.global_steps = UInt(70 * 3 * (50 + heuristic) / 2) // 2 * 70 -> 8 * 20 // FIXME
        parallelmc.mc.local_steps = UInt((25 + m_model!.num_movable_atoms) / 3)
        // MARK: Reduced from 2500 for faster testing
//        parallelmc.mc.global_steps = 25
//        parallelmc.mc.local_steps = 5
        
        parallelmc.mc.max_evals = UInt(max_evals)
        parallelmc.mc.min_rmsd = min_rmsd
        parallelmc.mc.num_saved_mins = UInt(n_poses)
        parallelmc.mc.hunt_cap = vec(10, 10, 10)
        parallelmc.num_tasks = UInt(exhaustiveness)
        parallelmc.num_threads = 8 // Maybe make this changable for iOS devices
        parallelmc.display_progress = (m_verbosity < .silent)
        parallelmc.mc.debugLevel = self.m_verbosity

        // Docking search
        print("Performing docking...")
        
        if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
            parallelmc.operate(m: m_model!, out: &poses, p: m_precalculated_byatom!, ig: m_grid!, corner1: m_grid!.corner1(), corner2: m_grid!.corner2(), progress_callback: m_progress_callback)
        } else {
            fatalError("Not Implemented Yet")
        }

        print("Done.")

        // Docking post-processing and rescoring
        poses = remove_redundant(in_: poses, min_rmsd: min_rmsd)

        if (!poses.isEmpty) {
            // For the Vina scoring function, we take the intramolecular energy from the best pose
            // the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
            if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
                // Refine poses if no_refine is false and got receptor
                if (!m_no_refine && m_receptor_initialized) {
                    var g: Change = Change(m_model!.get_size())
                    let quasi_newton_par: QuasiNewton = QuasiNewton(max_steps: UInt((25 + m_model!.num_movable_atoms) / 3))
                    var evalcount: Int = 0
                    let slope: fl = 1e6
                    m_non_cache!.slope = slope
                    for i in 0..<poses.count { 
                        for p in 0..<5 { 
                            m_non_cache!.slope = 100 * pow(10.0, 2.0 * fl(p))
                            quasi_newton_par.operate(&m_model!, m_precalculated_byatom!, m_non_cache!, &poses[i], &g, authentic_v, &evalcount)
                            if (m_non_cache!.within(m: m_model!)) {
                                if self.m_verbosity <= . log {
                                    print("Pose \(i) within non_cache!")
                                }
                                break
                            }
                        }
                        poses[i].coords = m_model!.get_heavy_atom_movable_coords()
                        m_non_cache!.slope = slope
                        // rescoring in case a ligand or flex sidechain atom is outside box
                        // ensuring poses will be sorted with the same slope (a.k.a. out of
                        // box penalty) that will be used to calculate final energies.
                        m_model!.set(&poses[i].c)
                        let all_grids = m_non_cache!.eval(m: m_model!, v: authentic_v[1])
                        let inter_pairs = m_model!.eval_inter(m_precalculated_byatom!, authentic_v) // ligand -- flex
                        let intra_pairs = m_model!.evalo(m_precalculated_byatom!, authentic_v)      // flex_i -- flex_i and flex_i -- flex_j
                        let lig_intra = m_model!.evali(m_precalculated_byatom!, v: authentic_v)     // ligand_i -- ligand_i
                        poses[i].e = all_grids + inter_pairs + intra_pairs + lig_intra
                    }
                }

                poses.sort() // order often changes after non_cache refinement
                m_model!.set(&poses[0].c)
                if (m_no_refine || !m_receptor_initialized) {
                    intramolecular_energy = m_model!.eval_intramolecular(m_precalculated_byatom!, ig: m_grid!, v: authentic_v)
                } else {
                    intramolecular_energy = m_model!.eval_intramolecular(m_precalculated_byatom!, ig: m_non_cache!, v: authentic_v)
                }       
            }

            for i in 0..<poses.count { 
                if m_verbosity <= .debug {
                    print("ENERGY FROM SEARCH: \(poses[i].e)")
                }
                m_model!.set(&poses[i].c)
                // For AD42 intramolecular_energy is equal to 0
                let energies = score(intramolecular_energy: intramolecular_energy)
                // Store energy components in current pose
                poses[i].e = energies[0] // specific to each scoring function
                poses[i].inter = energies[1] + energies[2]
                poses[i].intra = energies[3] + energies[4] + energies[5]
                poses[i].total = poses[i].inter + poses[i].intra // cost function for optimization
                poses[i].conf_independent = energies[6] // "torsion"
                poses[i].unbound = energies[7] // specific to each scoring function
                if m_verbosity <= .debug {
                    print("FINAL ENERGY: ")
                    show_score(energies: energies)
                }
            }

            // In AD4, the unbound energy is intra for each pose, so order may have changed
            // The order does not change in Vina because unbound is intra of the 1st pose
            if (m_sf_choice == .SF_AD42) {
                poses.sort()
            }
            
            // Now compute RMSD from the best model
            // Necessary to do it in two pass for AD4 scoring function
            m_model!.set(&poses[0].c)
            let best_model = m_model!.copy() as! Model

            if m_verbosity <= .log {
                print("mode |   affinity | dist from best mode")
                print("     | (kcal/mol) | rmsd l.b.| rmsd u.b.")
                print("-----+------------+----------+----------")
            }

            for i in 0..<poses.count { 
                m_model!.set(&poses[i].c)
                // Get RMSD between current pose and best_model
                poses[i].lb = m_model!.rmsd_lower_bound(best_model)
                poses[i].ub = m_model!.rmsd_upper_bound(best_model)

                if m_verbosity <= .log {
                    print("\(i + 1)    \(poses[i].e)  \(poses[i].lb)  \(poses[i].ub)")
                }
            }

            // Clean up by putting back the best pose in model
            m_model!.set(&poses[0].c)
        } else {
            throw VinaError("WARNING: Zero poses in output container after global search. This should not be happening and is likely a bug.\nWARNING: If possible, please file a bug report with your input files and random seed on GitHub.")
        }

        // Store results in Vina object
        m_poses = poses
    }

    func get_poses(_ how_many: Int = 9, energy_range: fl = 3.0) throws -> String { 

        var n: Int = 0
        var best_energy: fl = 0
        var out: String = ""
        var remarks: String = ""

        if how_many < 0 {
            throw VinaError("number of poses asked must be greater than zero.")
        }

        if energy_range < 0 {
            throw VinaError("energy range must be greater than zero.")
        }

        if (!m_poses.isEmpty) {
            // Get energy from the best conf
            best_energy = m_poses[0].e

            for i in 0..<m_poses.count { 
                /* Stop if:
                    - We wrote the number of conf asked
                    - If there is no conf to write
                    - The energy of the current conf is superior than best_energy + energy_range
                */
                if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range) { 
                    break
                }

                // Push the current pose to model
                m_model!.set(&m_poses[i].c)

                // Write conf
                remarks = vina_remarks(pose: m_poses[i], lb: m_poses[i].lb, ub: m_poses[i].ub)
                out += m_model!.write_model(model_number: sz(n + 1), remark: remarks)

                n += 1
            }

            // Push back the best conf in model
            m_model!.set(&m_poses[0].c)
        } else {
            throw VinaError("WARNING: Could not find any poses. No poses were written.")
        }

        return out
    }

    func get_poses_coordinates(_ how_many: Int = 9, energy_range: fl = 3.0) throws -> [[fl]] {
        var n: Int = 0
        var best_energy: fl = 0
        var coordinates: [[fl]] = []

        if how_many < 0 {
            throw VinaError("number of poses asked must be greater than zero.")
        }

        if energy_range < 0 {
            throw VinaError("energy range must be greater than zero.")
        }

        if (!m_poses.isEmpty) {
            // Get energy from the best conf
            best_energy = m_poses[0].e

            for i in 0..<m_poses.count { 
                /* Stop if:
                    - We wrote the number of conf asked
                    - If there is no conf to write
                    - The energy of the current conf is superior than best_energy + energy_range
                */
                if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range) { 
                    break
                }
                // Push the current pose to model
                m_model!.set(&m_poses[i].c)
                coordinates.append(m_model!.get_ligand_coords())
                n += 1
            }

            // Push back the best conf in model
            m_model!.set(&m_poses[0].c)
        } else {
            throw VinaError("WARNING: Could not find any pose coordinaates.")
        }

        return coordinates
    }
    
    func get_pose_energies(_ how_many: Int = 9, energy_range: fl = 3.0) throws -> [[fl]] {
        var n: Int = 0
        var best_energy: fl = 0
        var energies: [[fl]] = []
        
        if how_many < 0 {
            throw VinaError("number of poses asked must be greater than zero.")
        }
        
        if energy_range < 0 {
            throw VinaError("energy range must be greater than zero.")
        }
        
        if (!m_poses.isEmpty) {
            // Get energy from the best conf
            best_energy = m_poses[0].e
            
            for i in 0..<m_poses.count {
                /* Stop if:
                 - We wrote the number of conf asked
                 - If there is no conf to write
                 - The energy of the current conf is superior than best_energy + energy_range
                 */
                if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range) {
                    break // FIXME: check energy_range sanity?
                }
                
                // Push the current pose to model
                energies.append([m_poses[i].e,
                                 m_poses[i].inter, m_poses[i].intra,
                                 m_poses[i].conf_independent, m_poses[i].unbound])
                
                n += 1
            }
        } else {
            throw VinaError("WARNING: Could not find any pose energies.")
        }
        
        return energies
    }

    public func write_poses(output_name: String, _ how_many: Int = 9, _ energy_range: fl = 3.0) throws {
        if (!m_poses.isEmpty) {
            // Open output file with TextStreamer
            // maybe make into file path and use filename as input
            var f: TextStreamer = TextStreamer(filepath: output_name)
            let out = try get_poses(how_many, energy_range: energy_range)
            print(out, to: &f)
        } else {
            throw VinaError("WARNING: Could not find any poses. No poses were written.")
        }
    }
    
    public func write_pose(output_name: String, remark: String = "") { 
        var format_remark: String = ""

        // Add REMARK keyword to be PDB valid
        if(!remark.isEmpty){
            format_remark += "REMARK " + remark + " \n"
        }

        var f: TextStreamer = TextStreamer(filepath: output_name)
        m_model!.write_structure(&f, format_remark)
    }

    func write_maps(map_prefix: String = "receptor", gpf_filename: String? = nil,
                    fld_filename: String? = nil, receptor_filename: String? = nil) throws { 
        if !m_map_initialized {
            throw VinaError("Cannot write affinity maps. Affinity maps were not initialized.")
        }
        var atom_types: szv = []
        let atom_typing: AtomType.T = m_scoring_function!.get_atom_typing()

        if (m_ligand_initialized) {
            atom_types = m_model!.get_movable_atom_types(atom_typing)
        } else {
            atom_types = m_scoring_function!.get_atom_types()
        }

        if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
            print("Writing Vina maps")
            try m_grid!.write(out_prefix: map_prefix, atom_types: atom_types, gpf_filename: gpf_filename, fld_filename: fld_filename, receptor_filename: receptor_filename)
            print("Done writing Vina maps")
        } else {
            fatalError("Not Implemented Yet")
            // // Add electrostatics and desolvation maps
            // atom_types.append(AD_TYPE_SIZE)
            // atom_types.append(AD_TYPE_SIZE + 1)
            // print("Writing AD4.2 maps")
            // try m_ad4grid!.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename)
            // print("Done writing AD4.2 maps")
        }
        
    }

    public func show_score(energies: ContiguousArray<fl>) {
        print("Estimated Free Energy of Binding   : \(energies[0]) (kcal/mol) [=(1)+(2)+(3)-(4)]")
        print("(1) Final Intermolecular Energy    : \(energies[1] + energies[2]) (kcal/mol)")
        print("    Ligand - Receptor              : \(energies[1]) (kcal/mol)")
        print("    Ligand - Flex side chains      : \(energies[2]) (kcal/mol)")
        print("(2) Final Total Internal Energy    : \(energies[3] + energies[4] + energies[5]) (kcal/mol)")
        print("    Ligand                         : \(energies[5]) (kcal/mol)")
        print("    Flex   - Receptor              : \(energies[3]) (kcal/mol)")
        print("    Flex   - Flex side chains      : \(energies[4]) (kcal/mol)")
        print("(3) Torsional Free Energy          : \(energies[6]) (kcal/mol)")
        if (m_sf_choice == .SF_VINA || m_sf_choice == .SF_VINARDO) {
            print("(4) Unbound System's Energy        : \(energies[7]) (kcal/mol)")
        } else {
            print("(4) Unbound System's Energy [=(2)] : \(energies[7]) (kcal/mol)")
        }
    }
    
    public func debugModel() {
        m_model?.print_stuff()
    }
    
}
