#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

struct Mesh  // Fully Flattened 64-bit Safe Structure of Arrays (SoA) Mesh Container
{

    // 1. NODES
    std::vector<double>             node_coordinate_x;
    std::vector<double>             node_coordinate_y;
    std::vector<double>             node_coordinate_z;


    // 2. FACES
    std::vector<int64_t>            face_owner_cell_index;      // Left cell parent tracking index
    std::vector<int64_t>            face_neighbour_cell_index;  // Right cell parent tracking index (-1 for boundary)
    std::vector<double>             face_area;                  // Calculated via Newell's Method
    
    std::vector<double>             face_normal_x;              // Normalized normal vector component (X)
    std::vector<double>             face_normal_y;              // Normalized normal vector component (Y)
    std::vector<double>             face_normal_z;              // Normalized normal vector component (Z)

    std::vector<double>             face_center_x;              // Centroid spatial coordinate (X)
    std::vector<double>             face_center_y;              // Centroid spatial coordinate (Y)
    std::vector<double>             face_center_z;              // Centroid spatial coordinate (Z)

    std::vector<int64_t>            face_nodes;                 // Flat index sequence containing unique nodes per face
    std::vector<int64_t>            face_nodes_offsets;         // CSR offset limits array for face_nodes mapping


    // 3. CELLS
    std::vector<int>                cell_type;                  // 0 = TETRA_4, 1 = PYRA_5, 2 = PENTA_6, 3 = HEXA_8
    std::vector<double>             cell_volume;                // Evaluated precisely via sub-pyramid summation
    
    std::vector<double>             cell_center_x;              // Perfect 3D centroid position coordinate (X)
    std::vector<double>             cell_center_y;              // Perfect 3D centroid position coordinate (Y)
    std::vector<double>             cell_center_z;              // Perfect 3D centroid position coordinate (Z)

    std::vector<int64_t>            cell_nodes;                 // Flat vector list grouping vertex numbers per cell
    std::vector<int64_t>            cell_nodes_offsets;         // CSR offset mapping boundaries for cell_nodes vector

    std::vector<int64_t>            cell_faces;                 // Flat vector list grouping unique face IDs enclosing cell
    std::vector<int64_t>            cell_faces_offsets;         // CSR offset mapping boundaries for cell_faces vector


    // 4. EXPLICIT 2D BOUNDARY ELEMENTS
    std::vector<int>                boco_element_type;          // 0 = TRI_3, 1 = QUAD_4
    std::vector<int64_t>            boco_element_nodes;         // Flat layout of nodes defining the raw 2D elements
    std::vector<int64_t>            boco_element_offsets;       // CSR tracking pointers array for boco_element_nodes
    int64_t                         boco_elements_global_start = 0; // Tracks the baseline file element counter sequence


    // 5. BOUNDARY CONDITIONS
    std::vector<std::string>        boundary_name;              // Names retrieved from file (e.g., "Inlet", "Wall")
    std::vector<int>                boundary_type;              // Converted CGNS BCType_t enum identifiers
    std::vector<int64_t>            boundary_face_start_index;  // Real target front start ID inside solver's face arrays
    std::vector<int64_t>            boundary_face_count;        // The absolute tally count of elements defining this BC patch
    
    // Explicit lists tracking raw mesh descriptors prior to face generator mapping
    std::vector<std::vector<int64_t>> boundary_explicit_elements; // Retains data from files utilizing PointList mapping
    std::vector<int>                  boundary_grid_location;    // Stores native source CGNS GridLocation_t states
};
