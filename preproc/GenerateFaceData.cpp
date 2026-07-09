#include "PreProcessor.h"
#include "cgnslib.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <map>

struct FaceHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t hash = 0; // Tracks the running hash calculation
        for (int id : v) {
            hash ^= std::hash<int>()(id) + 0x9e3779b9 + (hash << 6) + (hash >> 2); // Combine node IDs using bit-mixing
        }
        return hash; // Yields a uniform hash for the unique node combination
    }
};

struct FaceHash64 {
    std::size_t operator()(const std::vector<int64_t>& v) const {
        std::size_t hash = 0;
        for (int64_t id : v) {
            hash ^= std::hash<int64_t>()(id) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

struct FaceMatchData {
    int global_face_index = -1; // Tracking ID assigned to the face globally
    int owner_cell = -1; // First element that indexes and claims this face
    int neighbor_cell = -1; // Opposite element; stays -1 for external boundaries
    std::vector<int64_t> original_nodes; // Keeps vertex sequence intact for geometric orientation
};

void PreProcessor::GenerateFaceData(Mesh& mesh) {
    // Ensure we are working with modern 64-bit safe types matching our upgraded Mesh struct
    int64_t num_cells = static_cast<int64_t>(mesh.cell_type.size()); 
    std::unordered_map<std::vector<int64_t>, FaceMatchData, FaceHash64> face_registry; 
    int64_t face_counter = 0; 

    auto get_local_faces = [](int cell_type, const std::vector<int64_t>& cell_nodes_list) -> std::vector<std::vector<int64_t>> {
        if (cell_type == 0) { // TETRA_4
            return {
                {cell_nodes_list[0], cell_nodes_list[2], cell_nodes_list[1]},
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[3]},
                {cell_nodes_list[1], cell_nodes_list[2], cell_nodes_list[3]},
                {cell_nodes_list[2], cell_nodes_list[0], cell_nodes_list[3]}
            };
        } else if (cell_type == 1) { // PYRA_5
            return {
                {cell_nodes_list[0], cell_nodes_list[3], cell_nodes_list[2], cell_nodes_list[1]},
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[4]},
                {cell_nodes_list[1], cell_nodes_list[2], cell_nodes_list[4]},
                {cell_nodes_list[2], cell_nodes_list[3], cell_nodes_list[4]},
                {cell_nodes_list[3], cell_nodes_list[0], cell_nodes_list[4]}
            };
        } else if (cell_type == 2) { // PENTA_6
            return {
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[2]},
                {cell_nodes_list[3], cell_nodes_list[5], cell_nodes_list[4]},
                {cell_nodes_list[0], cell_nodes_list[2], cell_nodes_list[5], cell_nodes_list[3]},
                {cell_nodes_list[1], cell_nodes_list[0], cell_nodes_list[3], cell_nodes_list[4]},
                {cell_nodes_list[2], cell_nodes_list[1], cell_nodes_list[4], cell_nodes_list[5]}
            };
        } else { // HEXA_8
            return {
                {cell_nodes_list[0], cell_nodes_list[3], cell_nodes_list[2], cell_nodes_list[1]},
                {cell_nodes_list[4], cell_nodes_list[5], cell_nodes_list[6], cell_nodes_list[7]},
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[5], cell_nodes_list[4]},
                {cell_nodes_list[1], cell_nodes_list[2], cell_nodes_list[6], cell_nodes_list[5]},
                {cell_nodes_list[2], cell_nodes_list[3], cell_nodes_list[7], cell_nodes_list[6]},
                {cell_nodes_list[3], cell_nodes_list[0], cell_nodes_list[4], cell_nodes_list[7]}
            };
        }
    };

    // 2. PRE-PROCESS THE READER'S 2D SURFACE ELEMENTS FOR COMPATIBILITY MATCHING
    // Maps sorted face nodes -> global file 2D element index
    std::unordered_map<std::vector<int64_t>, int64_t, FaceHash64> boco_elements_lookup;
    int64_t total_boco_elements = static_cast<int64_t>(mesh.boco_element_type.size());
    
    for (int64_t e = 0; e < total_boco_elements; ++e) {
        int64_t start = mesh.boco_element_offsets[e];
        int64_t end = mesh.boco_element_offsets[e + 1];
        std::vector<int64_t> boco_nodes(mesh.boco_element_nodes.begin() + start, mesh.boco_element_nodes.begin() + end);
        std::sort(boco_nodes.begin(), boco_nodes.end());
        
        // Map back to its original 0-based file tracking location identifier
        boco_elements_lookup[boco_nodes] = mesh.boco_elements_global_start + e;
    }

    // 3. Populate the Face Registry from 3D Cells
    for (int64_t c = 0; c < num_cells; ++c) {
        int64_t start_idx = mesh.cell_nodes_offsets[c]; 
        int64_t end_idx = mesh.cell_nodes_offsets[c + 1]; 
        std::vector<int64_t> cell_nodes_list(mesh.cell_nodes.begin() + start_idx, mesh.cell_nodes.begin() + end_idx); 
        auto local_faces = get_local_faces(mesh.cell_type[c], cell_nodes_list); 

        for (auto& face_nodes : local_faces) {
            std::vector<int64_t> sorted_nodes = face_nodes; 
            std::sort(sorted_nodes.begin(), sorted_nodes.end()); 
            auto it = face_registry.find(sorted_nodes); 
            if (it == face_registry.end()) {
                FaceMatchData data; 
                data.global_face_index = face_counter++; 
                data.owner_cell = c; 
                data.neighbor_cell = -1; 
                data.original_nodes = face_nodes; 
                face_registry[sorted_nodes] = data; 
            } else {
                it->second.neighbor_cell = c; 
            }
        }
    }

    // 4. SORT AND GROUP BOUNDARY FACES ACCORDING TO THEIR ASSIGNED NAMED PATCHES
    std::vector<FaceMatchData> internal_faces;
    // Map: BC Patch Index -> Vector of FaceMatchData belonging to it
    std::map<size_t, std::vector<FaceMatchData>> boundary_patches_registry;
    std::vector<FaceMatchData> unmapped_boundary_faces; // Fallback container

    size_t num_bcos = mesh.boundary_name.size();

    for (auto& [sorted_key, data] : face_registry) {
        if (data.neighbor_cell != -1) {
            internal_faces.push_back(data);
        } else {
            // It's an exterior boundary face. Find its 2D element reference ID from the file node layout
            auto lookup_it = boco_elements_lookup.find(sorted_key);
            bool matched_to_patch = false;

            if (lookup_it != boco_elements_lookup.end()) {
                int64_t file_element_id = lookup_it->second;

                // Query which boundary condition block owns this file element identifier
                for (size_t b = 0; b < num_bcos; ++b) {
                    if (mesh.boundary_face_start_index[b] != -1) {
                        // Case A: Contiguous PointRange allocation check
                        int64_t start = mesh.boundary_face_start_index[b];
                        int64_t count = mesh.boundary_face_count[b];
                        if (file_element_id >= start && file_element_id < (start + count)) {
                            boundary_patches_registry[b].push_back(data);
                            matched_to_patch = true;
                            break;
                        }
                    } else {
                        // Case B: Scattered PointList lookups check
                        const auto& explicit_list = mesh.boundary_explicit_elements[b];
                        if (std::find(explicit_list.begin(), explicit_list.end(), file_element_id) != explicit_list.end()) {
                            boundary_patches_registry[b].push_back(data);
                            matched_to_patch = true;
                            break;
                        }
                    }
                }
            }

            if (!matched_to_patch) {
                unmapped_boundary_faces.push_back(data); // Kept safe for free surfaces or symmetry voids
            }
        }
    }

    // 5. ASSEMBLE CONTIGUOUS ORDERED FACES STREAM (Internal front, grouped boundaries back)
    std::vector<FaceMatchData> ordered_faces = std::move(internal_faces);
    int64_t num_internal = static_cast<int64_t>(ordered_faces.size());

    // Re-index ranges for the mesh metadata arrays cleanly
    for (size_t b = 0; b < num_bcos; ++b) {
        if (mesh.boundary_face_start_index[b] != -1) {
            // Re-map the continuous block start pointing directly inside our solver's face arrays
            mesh.boundary_face_start_index[b] = static_cast<int64_t>(ordered_faces.size());
        }
        
        const auto& patch_faces = boundary_patches_registry[b];
        mesh.boundary_face_count[b] = static_cast<int64_t>(patch_faces.size());
        
        // Insert them sequentially into the flat ordered stream array
        ordered_faces.insert(ordered_faces.end(), patch_faces.begin(), patch_faces.end());
    }
    
    // Push remaining unmapped boundary faces to the absolute tail end
    ordered_faces.insert(ordered_faces.end(), unmapped_boundary_faces.begin(), unmapped_boundary_faces.end());

    int64_t total_faces = static_cast<int64_t>(ordered_faces.size());

    // 6. Populate Flat Mesh Fields, Calculate Newell Normals & Area Metas
    mesh.face_owner_cell_index.resize(total_faces); 
    mesh.face_neighbour_cell_index.resize(total_faces); 
    mesh.face_area.resize(total_faces); 
    mesh.face_normal_x.resize(total_faces); mesh.face_normal_y.resize(total_faces); mesh.face_normal_z.resize(total_faces); 
    mesh.face_center_x.resize(total_faces); mesh.face_center_y.resize(total_faces); mesh.face_center_z.resize(total_faces); 
    mesh.face_nodes_offsets.resize(total_faces + 1, 0); 
    mesh.face_nodes.reserve(total_faces * 4); 

    std::vector<std::vector<int64_t>> cell_to_faces_temp(num_cells); 
    int64_t current_node_offset = 0; 

    for (int64_t f = 0; f < total_faces; ++f) {
        const auto& f_data = ordered_faces[f]; 
        mesh.face_owner_cell_index[f] = f_data.owner_cell; 
        mesh.face_neighbour_cell_index[f] = f_data.neighbor_cell; 
        cell_to_faces_temp[f_data.owner_cell].push_back(f); 
        if (f_data.neighbor_cell != -1) {
            cell_to_faces_temp[f_data.neighbor_cell].push_back(f); 
        }

        for (int64_t node : f_data.original_nodes) mesh.face_nodes.push_back(node); 
        current_node_offset += static_cast<int64_t>(f_data.original_nodes.size()); 
        mesh.face_nodes_offsets[f + 1] = current_node_offset; 

        int64_t num_face_nodes = static_cast<int64_t>(f_data.original_nodes.size()); 
        double cx = 0.0, cy = 0.0, cz = 0.0; 
        for (int64_t node : f_data.original_nodes) {
            cx += mesh.node_coordinate_x[node]; cy += mesh.node_coordinate_y[node]; cz += mesh.node_coordinate_z[node]; 
        }
        mesh.face_center_x[f] = cx / num_face_nodes; mesh.face_center_y[f] = cy / num_face_nodes; mesh.face_center_z[f] = cz / num_face_nodes; 

        double nx = 0.0, ny = 0.0, nz = 0.0; 
        for (int64_t i = 0; i < num_face_nodes; ++i) {
            int64_t n0 = f_data.original_nodes[i], n1 = f_data.original_nodes[(i + 1) % num_face_nodes]; 
            nx += (mesh.node_coordinate_y[n0] - mesh.node_coordinate_y[n1]) * (mesh.node_coordinate_z[n0] + mesh.node_coordinate_z[n1]); 
            ny += (mesh.node_coordinate_z[n0] - mesh.node_coordinate_z[n1]) * (mesh.node_coordinate_x[n0] + mesh.node_coordinate_x[n1]); 
            nz += (mesh.node_coordinate_x[n0] - mesh.node_coordinate_x[n1]) * (mesh.node_coordinate_y[n0] + mesh.node_coordinate_y[n1]); 
        }
        double mag = std::sqrt(nx*nx + ny*ny + nz*nz); 
        if (mag > 1e-14) {
            mesh.face_area[f] = 0.5 * mag; mesh.face_normal_x[f] = nx/mag; mesh.face_normal_y[f] = ny/mag; mesh.face_normal_z[f] = nz/mag; 
        }
    }

    mesh.cell_faces_offsets.resize(num_cells + 1, 0); 
    int64_t current_cf_offset = 0; 
    for (int64_t c = 0; c < num_cells; ++c) {
        for (int64_t face_id : cell_to_faces_temp[c]) {
            mesh.cell_faces.push_back(face_id); 
        }
        current_cf_offset += static_cast<int64_t>(cell_to_faces_temp[c].size()); 
        mesh.cell_faces_offsets[c + 1] = current_cf_offset; 
    }
}
