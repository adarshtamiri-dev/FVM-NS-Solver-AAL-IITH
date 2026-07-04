#pragma once

#include <iostream>
#include <string>
#include <vector>


struct Mesh  // Fully Flattened SoA Mesh Container
{
    std::vector<double>             node_coordinate_x;
    std::vector<double>             node_coordinate_y;
    std::vector<double>             node_coordinate_z;

    // 2. FACES
    std::vector<int>                face_owner_cell_index;
    std::vector<int>                face_neighbour_cell_index;
    std::vector<double>             face_area;                 
    
    std::vector<double>             face_normal_x;
    std::vector<double>             face_normal_y;
    std::vector<double>             face_normal_z;

    std::vector<double>             face_center_x;
    std::vector<double>             face_center_y;
    std::vector<double>             face_center_z;

    std::vector<int>                face_nodes;           
    std::vector<int>                face_nodes_offsets;   

    // 3. CELLS
    std::vector<int>                cell_type;            
    std::vector<double>             cell_volume;
    
    std::vector<double>             cell_center_x;
    std::vector<double>             cell_center_y;
    std::vector<double>             cell_center_z;

    std::vector<int>                cell_nodes;
    std::vector<int>                cell_nodes_offsets;

    std::vector<int>                cell_faces;
    std::vector<int>                cell_faces_offsets;   

    // 4. BOUNDARIES
    std::vector<std::string>        boundary_name;
    std::vector<int>                boundary_type;       
    std::vector<int>                boundary_face_start_index;
    std::vector<int>                boundary_face_count;
};

// --- RAW INT REPRESENTATIONS FOR CELL TYPES ---
// 0 = TETRA_4, 1 = PYRA_5, 2 = PENTA_6, 3 = HEXA_8, 4 = POLYHEDRON
// --- RAW INT REPRESENTATIONS FOR BOUNDARY TYPES ---
// 0 = WALL, 1 = INFLOW, 2 = OUTFLOW, 3 = SYMMETRY, 4 = FARFIELD