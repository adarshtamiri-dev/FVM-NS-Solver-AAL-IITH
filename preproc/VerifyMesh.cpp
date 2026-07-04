#include "PreProcessor.h"

#include <iostream>  // For std::cout
#include <iomanip>   // For std::setw, std::left, and std::scientific
#include <vector>    // For std::vector (handling cell types and counts)
#include <string>    // For std::string and std::to_string
#include <cctype>    // For std::toupper (used in name parsing)

void PreProcessor::VerifyMesh(const Mesh& mesh) {
    int num_nodes = static_cast<int>(mesh.node_coordinate_x.size()); // Total node array size
    int num_cells = static_cast<int>(mesh.cell_type.size()); // Total cell topology size
    int num_faces = static_cast<int>(mesh.face_area.size()); // Total processed face count

    std::cout << "\n\n    Detail                                Value\n"; // Consistent master column headers setup
    std::cout << "________________________________________________________\n\n"; // Visual boundary line separator

    // --- [1] Basic Grid Topology ---
    std::cout << "[1] Basic Grid Topology\n";
    std::cout << "      " << std::left << std::setw(36) << "Number of Nodes" << num_nodes << "\n"; // Align data with the primary value column
    std::cout << "      " << std::left << std::setw(36) << "Number of Cells" << num_cells << "\n"; // Align data with the primary value column
    std::cout << "      " << std::left << std::setw(36) << "Number of Faces" << num_faces << "\n\n"; // Align data with the primary value column

    // --- [2] Cell Type Breakdown ---
    std::vector<std::string> cell_type_names = { "Tetrahedral (TETRA)", "Pyramid (PYRA)", "Prism (PENTA)", "Hexahedral (HEXA)", "Polyhedron (POLY)" };
    std::vector<int> cell_type_counts(5, 0); // Histogram array initialization for element shapes
    for (int type : mesh.cell_type) { 
        if (type >= 0 && type < 5) cell_type_counts[type]++; // Increment matching bucket counter
    }

    std::cout << "[2] Cell Type Breakdown\n";
    for (int i = 0; i < 5; ++i) {
        if (cell_type_counts[i] > 0) {
            std::cout << "      " << std::left << std::setw(36) << cell_type_names[i] << cell_type_counts[i] << "\n"; // Print element counts under the same alignment value column
        }
    }
    std::cout << "\n";

    // --- [3] Geometric Volume Assessment ---
    double min_vol = 1e30, max_vol = -1e30; // Establish baseline comparison scale limits
    bool volumes_strictly_positive = true; // Inverted flag tracking negative or degenerate cells
    for (double vol : mesh.cell_volume) {
        if (vol < min_vol) min_vol = vol; // Lower limit bound capture update
        if (vol > max_vol) max_vol = vol; // Upper limit bound capture update
        if (vol <= 0.0) volumes_strictly_positive = false; // Mark failure flag if an invalid cell volume occurs
    }

    std::cout << "[3] Geometric Volume Assessment\n";
    std::cout << "      " << std::left << std::setw(36) << "Minimum Cell Volume" << std::scientific << min_vol << "\n"; // Print geometric properties aligned nicely
    std::cout << "      " << std::left << std::setw(36) << "Maximum Cell Volume" << std::scientific << max_vol << "\n"; // Print geometric properties aligned nicely
    std::cout << "      " << std::left << std::setw(36) << "Non-Negative Volume Test" << (volumes_strictly_positive ? "PASSED" : "FAILED") << "\n\n"; // Verification test evaluation status alignment

    // --- [4] Boundary Patch Inventory ---
    std::cout << "[4] Defined Boundary Patches\n";
    std::cout << "      " << std::left << std::setw(36) << "Boundary Type" << "Face Count\n"; // Standardized column labels removing the raw CGNS field entirely

    int total_reconstructed_boundary_faces = 0; // Accumulator verifying global boundary components
    for (size_t b = 0; b < mesh.boundary_name.size(); ++b) {
        int count = mesh.boundary_face_count[b]; // Grab face count within current boundary patch index
        total_reconstructed_boundary_faces += count; // Accumulate global bounds face counts
        
        std::string cgns_type_str = "BCTypeNull/UserDefined"; // Fallback string attribute designation
        int internal_id = -1; // Fallback category integer flag identification indicator
        std::string internal_name_str = "UNKNOWN"; // Fallback identifier text label string

        std::string upper_name = mesh.boundary_name[b]; // Create mutable copy string for uppercase processing
        for (auto &c : upper_name) c = toupper(c); // Force conversion characters string transformation loops

        if (upper_name.find("INLET") != std::string::npos) {
            cgns_type_str = "BCInflow";
            internal_id = 1;
            internal_name_str = "INFLOW";
        } 
        else if (upper_name.find("OUTLET") != std::string::npos) {
            cgns_type_str = "BCOutflow";
            internal_id = 2;
            internal_name_str = "OUTFLOW";
        } 
        else if (upper_name.find("WALL") != std::string::npos) {
            cgns_type_str = "BCWall";
            internal_id = 0;
            internal_name_str = "WALL";
        } 
        else if (upper_name.find("FARFIELD") != std::string::npos) {
            cgns_type_str = "BCFarfield";
            internal_id = 4;
            internal_name_str = "FARFIELD";
        } 
        else if (upper_name.find("SYM") != std::string::npos) {
            cgns_type_str = "BCSymmetry";
            internal_id = 3;
            internal_name_str = "SYMMETRY";
        }

        int raw_file_enum = mesh.boundary_type[b]; // Load original native input file identification codes
        if (raw_file_enum == 7) {
            cgns_type_str = "BCFarfield";
            internal_id = 4;
            internal_name_str = "FARFIELD";
        } 
        else if (raw_file_enum == 9 || raw_file_enum == 10 || raw_file_enum == 11) {
            cgns_type_str = "BCInflow";
            internal_id = 1;
            internal_name_str = "INFLOW";
        } 
        else if (raw_file_enum == 13 || raw_file_enum == 14 || raw_file_enum == 15) {
            cgns_type_str = "BCOutflow";
            internal_id = 2;
            internal_name_str = "OUTFLOW";
        } 
        else if (raw_file_enum == 16 || raw_file_enum == 17) {
            cgns_type_str = "BCSymmetryPlane";
            internal_id = 3;
            internal_name_str = "SYMMETRY";
        } 
        else if (raw_file_enum == 19 || raw_file_enum == 20 || raw_file_enum == 21) {
            cgns_type_str = "BCWall";
            internal_id = 0;
            internal_name_str = "WALL";
        }

        if (internal_id == -1) {
            internal_name_str = "UNKNOWN, " + std::to_string(raw_file_enum); // Append missing tracking indicator code digits string
        }

        std::cout << "      " << std::left << std::setw(36) << internal_name_str << count << "\n"; // Clean, structured face counts print statement execution
    }
    std::cout << "________________________________________________________\n"; // Final table closure boundary design highlight line
}