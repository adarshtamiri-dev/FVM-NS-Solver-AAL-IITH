#include "PreProcessor.h"
#include "cgnslib.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <numeric>

struct FaceHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t hash = 0; // Tracks the running hash calculation
        for (int id : v) {
            hash ^= std::hash<int>()(id) + 0x9e3779b9 + (hash << 6) + (hash >> 2); // Combine node IDs using bit-mixing
        }
        return hash; // Yields a uniform hash for the unique node combination
    }
};

struct FaceMatchData {
    int global_face_index = -1; // Tracking ID assigned to the face globally
    int owner_cell = -1; // First element that indexes and claims this face
    int neighbor_cell = -1; // Opposite element; stays -1 for external boundaries
    std::vector<int> original_nodes; // Keeps vertex sequence intact for geometric orientation
};

void PreProcessor::GenerateFaceData(Mesh& mesh) {
    int num_cells = static_cast<int>(mesh.cell_type.size()); // Grab total element count from mesh storage
    std::unordered_map<std::vector<int>, FaceMatchData, FaceHash> face_registry; // Unique map using sorted keys to find pairs
    int face_counter = 0; // Incremental tracker to hand out unique face IDs

    auto get_local_faces = [](int cell_type, const std::vector<int>& cell_nodes_list) -> std::vector<std::vector<int>> {
        if (cell_type == 0) { // TETRA_4 element mapping definition
            return {
                {cell_nodes_list[0], cell_nodes_list[2], cell_nodes_list[1]},
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[3]},
                {cell_nodes_list[1], cell_nodes_list[2], cell_nodes_list[3]},
                {cell_nodes_list[2], cell_nodes_list[0], cell_nodes_list[3]}
            };
        } else if (cell_type == 1) { // PYRA_5 topology definition
            return {
                {cell_nodes_list[0], cell_nodes_list[3], cell_nodes_list[2], cell_nodes_list[1]},
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[4]},
                {cell_nodes_list[1], cell_nodes_list[2], cell_nodes_list[4]},
                {cell_nodes_list[2], cell_nodes_list[3], cell_nodes_list[4]},
                {cell_nodes_list[3], cell_nodes_list[0], cell_nodes_list[4]}
            };
        } else if (cell_type == 2) { // PENTA_6 (Prisms) mapping layout
            return {
                {cell_nodes_list[0], cell_nodes_list[1], cell_nodes_list[2]},
                {cell_nodes_list[3], cell_nodes_list[5], cell_nodes_list[4]},
                {cell_nodes_list[0], cell_nodes_list[2], cell_nodes_list[5], cell_nodes_list[3]},
                {cell_nodes_list[1], cell_nodes_list[0], cell_nodes_list[3], cell_nodes_list[4]},
                {cell_nodes_list[2], cell_nodes_list[1], cell_nodes_list[4], cell_nodes_list[5]}
            };
        } else { // HEXA_8 brick element arrangement layout
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

    for (int c = 0; c < num_cells; ++c) {
        int start_idx = mesh.cell_nodes_offsets[c]; // CSR start location marker for cell data
        int end_idx = mesh.cell_nodes_offsets[c + 1]; // CSR end location marker for cell data
        std::vector<int> cell_nodes_list(mesh.cell_nodes.begin() + start_idx, mesh.cell_nodes.begin() + end_idx); // Extract local connectivity array chunk
        auto local_faces = get_local_faces(mesh.cell_type[c], cell_nodes_list); // Decompose current element configuration into faces

        for (auto& face_nodes : local_faces) {
            std::vector<int> sorted_nodes = face_nodes; // Work with a duplicate for sorting tasks
            std::sort(sorted_nodes.begin(), sorted_nodes.end()); // Enforce invariant sorting order to form unique lookup key
            auto it = face_registry.find(sorted_nodes); // Query lookup index table for existing shared side match
            if (it == face_registry.end()) {
                FaceMatchData data; // First match instantiation branch
                data.global_face_index = face_counter++; // Generate unique index sequentially
                data.owner_cell = c; // Record element as base layout owner
                data.neighbor_cell = -1; // Unbound marker fallback position
                data.original_nodes = face_nodes; // Keep raw loop direction for geometric normals math
                face_registry[sorted_nodes] = data; // Cache newly built metadata record in grid index map
            } else {
                it->second.neighbor_cell = c; // Face already logged, link this current cell as neighbor
            }
        }
    }

    std::vector<FaceMatchData> internal_faces, boundary_faces; // Temporary tracking arrays for separation step
    for (const auto& [key, data] : face_registry) {
        if (data.neighbor_cell == -1) boundary_faces.push_back(data); // Unmatched neighbor tracks go into boundary collection
        else internal_faces.push_back(data); // Internal connection match confirmed
    }

    int num_internal = static_cast<int>(internal_faces.size()); // Final internal interface totals count
    int num_boundary = static_cast<int>(boundary_faces.size()); // Final external surface boundary totals count
    int total_faces = num_internal + num_boundary; // Combined total structural grid lines count

    mesh.face_owner_cell_index.resize(total_faces); // Sizing contiguous storage array for face owners
    mesh.face_neighbour_cell_index.resize(total_faces); // Sizing contiguous storage array for face neighbors
    mesh.face_area.resize(total_faces); // Size geometric face metric array
    mesh.face_normal_x.resize(total_faces); mesh.face_normal_y.resize(total_faces); mesh.face_normal_z.resize(total_faces); // Normalize vector component channels sizing
    mesh.face_center_x.resize(total_faces); mesh.face_center_y.resize(total_faces); mesh.face_center_z.resize(total_faces); // Grid positions array allocation
    mesh.face_nodes_offsets.resize(total_faces + 1, 0); // Setup CSR offset structures for tracking face node bounds
    mesh.face_nodes.reserve(total_faces * 4); // Pre-allocate assuming a quad-dominant geometry distribution pattern

    std::vector<FaceMatchData> ordered_faces = std::move(internal_faces); // Reorder stream so internal faces occupy the front of the array
    ordered_faces.insert(ordered_faces.end(), boundary_faces.begin(), boundary_faces.end()); // Chain exterior layout elements directly down to the backend

    std::vector<std::vector<int>> cell_to_faces_temp(num_cells); // Inverse transient tracker linking cells back to their structural faces
    int current_node_offset = 0; // Flat tracking accumulator tracking array indexes position

    for (int f = 0; f < total_faces; ++f) {
        const auto& f_data = ordered_faces[f]; // Access clean structured mapping reference
        mesh.face_owner_cell_index[f] = f_data.owner_cell; // Save master parent element index relationship
        mesh.face_neighbour_cell_index[f] = f_data.neighbor_cell; // Save opposite neighbor element index mapping relationship
        cell_to_faces_temp[f_data.owner_cell].push_back(f); // Log face index down onto the primary owner lookup tracker
        if (f_data.neighbor_cell != -1) {
            cell_to_faces_temp[f_data.neighbor_cell].push_back(f); // Log link mapping relationship inside the neighbor cell roster
        }

        for (int node : f_data.original_nodes) mesh.face_nodes.push_back(node); // Stack structural IDs sequentially into flat array
        current_node_offset += static_cast<int>(f_data.original_nodes.size()); // Tally total coordinates handled for track updates
        mesh.face_nodes_offsets[f + 1] = current_node_offset; // Record offset break boundaries index limit

        int num_face_nodes = static_cast<int>(f_data.original_nodes.size()); // Cache local vertex count sizing profile
        double cx = 0.0, cy = 0.0, cz = 0.0; // Clear coordinate sum registers for center location calculations
        for (int node : f_data.original_nodes) {
            cx += mesh.node_coordinate_x[node]; cy += mesh.node_coordinate_y[node]; cz += mesh.node_coordinate_z[node]; // Accumulate vertex tracking coordinates
        }
        mesh.face_center_x[f] = cx / num_face_nodes; mesh.face_center_y[f] = cy / num_face_nodes; mesh.face_center_z[f] = cz / num_face_nodes; // Center location math calculations

        double nx = 0.0, ny = 0.0, nz = 0.0; // Clear geometric accumulator registers for Newell's area calculation formula
        for (int i = 0; i < num_face_nodes; ++i) {
            int n0 = f_data.original_nodes[i], n1 = f_data.original_nodes[(i + 1) % num_face_nodes]; // Wrap indices cleanly via modular operations
            nx += (mesh.node_coordinate_y[n0] - mesh.node_coordinate_y[n1]) * (mesh.node_coordinate_z[n0] + mesh.node_coordinate_z[n1]); // Vector step math mapping loop X
            ny += (mesh.node_coordinate_z[n0] - mesh.node_coordinate_z[n1]) * (mesh.node_coordinate_x[n0] + mesh.node_coordinate_x[n1]); // Vector step math mapping loop Y
            nz += (mesh.node_coordinate_x[n0] - mesh.node_coordinate_x[n1]) * (mesh.node_coordinate_y[n0] + mesh.node_coordinate_y[n1]); // Vector step math mapping loop Z
        }
        double mag = std::sqrt(nx*nx + ny*ny + nz*nz); // Calculate absolute directional vector length magnitude
        if (mag > 1e-14) {
            mesh.face_area[f] = 0.5 * mag; mesh.face_normal_x[f] = nx/mag; mesh.face_normal_y[f] = ny/mag; mesh.face_normal_z[f] = nz/mag; // Normalise data elements profiles
        }
    }

    mesh.cell_faces_offsets.resize(num_cells + 1, 0); // Initialize cell-to-face CSR layout maps
    int current_cf_offset = 0; // Cumulative tracker tracking indexing bounds
    for (int c = 0; c < num_cells; ++c) {
        for (int face_id : cell_to_faces_temp[c]) {
            mesh.cell_faces.push_back(face_id); // Flatten nested list references down into linear structures space layout
        }
        current_cf_offset += static_cast<int>(cell_to_faces_temp[c].size()); // Advance tracking values count indices limits
        mesh.cell_faces_offsets[c + 1] = current_cf_offset; // Commit boundaries tracking details markers
    }

    int current_search_start = num_internal; // Boundary structures index starts directly where inner channels end
    for (size_t b = 0; b < mesh.boundary_name.size(); ++b) {
        mesh.boundary_face_start_index[b] = current_search_start; // Map absolute location base markers positions
        current_search_start += mesh.boundary_face_count[b]; // Push cursor limits out cleanly for the next physical patch blocks tracking
    }
}