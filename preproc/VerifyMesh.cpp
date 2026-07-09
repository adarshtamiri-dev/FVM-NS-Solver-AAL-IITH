#include "PreProcessor.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <limits>

void PreProcessor::VerifyMesh(const Mesh& mesh) {
    // --- 1. Basic Grid Topology Tally ---
    int64_t num_nodes = static_cast<int64_t>(mesh.node_coordinate_x.size());
    int64_t num_cells = static_cast<int64_t>(mesh.cell_type.size());
    int64_t num_faces = static_cast<int64_t>(mesh.face_owner_cell_index.size());

    // --- 2. Cell Type Breakdown ---
    int64_t num_tetra = 0;
    int64_t num_pyra  = 0;
    int64_t num_penta = 0;
    int64_t num_hexa  = 0;

    for (int type : mesh.cell_type) {
        if (type == 0)      num_tetra++;
        else if (type == 1) num_pyra++;
        else if (type == 2) num_penta++;
        else if (type == 3) num_hexa++;
    }

    // --- 3. Geometric Volume Assessment ---
    double min_vol = std::numeric_limits<double>::max();
    double max_vol = -std::numeric_limits<double>::max();
    bool passed_volume_test = true;

    for (double vol : mesh.cell_volume) {
        if (vol < min_vol) min_vol = vol;
        if (vol > max_vol) max_vol = vol;
        if (vol <= 0.0) {
            passed_volume_test = false;
        }
    }
    // Clean fallback if cell size is empty
    if(mesh.cell_volume.empty()) { min_vol = 0.0; max_vol = 0.0; passed_volume_test = false; }

    // --- 4. Render Validation Report ---
    std::cout << std::left << std::setw(35) << "Detail" << "Value\n";
    std::cout << std::string(60, '-') << "\n\n";

    std::cout << "[1] Basic Grid Topology\n";
    std::cout << "    " << std::setw(31) << "Number of Nodes" << num_nodes << "\n";
    std::cout << "    " << std::setw(31) << "Number of Cells" << num_cells << "\n";
    std::cout << "    " << std::setw(31) << "Number of Faces" << num_faces << "\n\n";

    std::cout << "[2] Cell Type Breakdown\n";
    if (num_tetra > 0) std::cout << "    " << std::setw(31) << "Tetrahedral (TETRA)" << num_tetra << "\n";
    if (num_pyra > 0)  std::cout << "    " << std::setw(31) << "Pyramidal (PYRA)"    << num_pyra  << "\n";
    if (num_penta > 0) std::cout << "    " << std::setw(31) << "Prismatic (PENTA)"   << num_penta << "\n";
    if (num_hexa > 0)  std::cout << "    " << std::setw(31) << "Hexahedral (HEXA)"    << num_hexa  << "\n";
    std::cout << "\n";

    std::cout << "[3] Geometric Volume Assessment\n";
    std::cout << "    " << std::setw(31) << "Minimum Cell Volume" << std::scientific << min_vol << "\n";
    std::cout << "    " << std::setw(31) << "Maximum Cell Volume" << std::scientific << max_vol << "\n";
    std::cout << std::defaultfloat; // Reset float formatting
    std::cout << "    " << std::setw(31) << "Non-Negative Volume Test" << (passed_volume_test ? "PASSED" : "FAILED") << "\n\n";

    std::cout << "[4] Defined Boundary Patches\n";
    std::cout << "    " << std::left << std::setw(27) << "Boundary Type" << "Face Count\n";
    
    // CGNS BCType enum mapper helper function string names
    auto get_bc_type_name = [](int type_id) -> std::string {
        switch (type_id) {
            case 2:  return "BCOutflow"; 
            case 3:  return "BCInflow";
            case 7:  return "BCFarfield";
            case 8:  return "BCWall";
            default: return "PATCH_" + std::to_string(type_id);
        }
    };

    for (size_t b = 0; b < mesh.boundary_name.size(); ++b) {
        // Render name clean or fallback to enum translation logic matching your style output
        std::string label = mesh.boundary_name[b];
        // Convert label to uppercase to cleanly match your visual criteria layout format
        std::transform(label.begin(), label.end(), label.begin(), ::toupper);

        std::cout << "    " << std::left << std::setw(27) << label << mesh.boundary_face_count[b] << "\n";
    }
    std::cout << std::string(60, '-') << "\n";
}
