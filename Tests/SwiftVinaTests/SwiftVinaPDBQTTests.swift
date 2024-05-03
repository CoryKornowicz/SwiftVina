import XCTest
@testable import SwiftVina

final class SwiftVinaPDBQTTests: XCTestCase {
    
    // Test Loading Ligand PDBQT File
    func testPDBQTParseLigandFromFile() throws {
        //file path
        let filePath: String = "/Users/corykornowicz/Downloads/1iep_ligand.pdbqt"
        let ligand: Model = try parse_ligand_pdbqt_from_file(filePath: filePath, atype: .XS)
        XCTAssert(!ligand.coords.isEmpty)
    }
    
    // Test Loading Ligand PDBQT File
    func testPDBQTParseLigandFromString() throws {
        //file path
        let filePath: String = "/Users/corykornowicz/Downloads/1iep_ligand.pdbqt"
        let fileString = try NSString(contentsOfFile: filePath, encoding: NSUTF8StringEncoding) as String
        let ligand: Model = try parse_ligand_pdbqt_from_string(stringContents: fileString, atype: .XS)
        XCTAssert(!ligand.coords.isEmpty)
    }
    
    // Test Loading Ligand PDBQT File
    func testPDBQTParseReceptorFromFile() throws {
        //file path
        let filePath: String = "/Users/corykornowicz/Downloads/1iep_receptor.pdbqt"
        let receptor: Model = try parse_receptor_pdbqt(rigid: filePath, flex: nil)
        XCTAssert(!receptor.grid_atoms.map(\.coords).isEmpty)
    }
    
}
