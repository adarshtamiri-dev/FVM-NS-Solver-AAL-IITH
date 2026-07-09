#include "cgnslib.h"
#include "PreProcessor.h"

bool PreProcessor::CGNSReader(const std::string& filename, Mesh& mesh) {
    int file_index = 0; 

    // Open file in read mode
    if (cg_open(filename.c_str(), CG_MODE_READ, &file_index) != CG_OK) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return false; 
    }

    int base_id = 1; 
    int zone_id = 1; // Enforced single-zone tracking for now
    char zone_name[33]; 
    cgsize_t zone_sizes[3] = {0, 0, 0}; 

    if (cg_zone_read(file_index, base_id, zone_id, zone_name, zone_sizes) != CG_OK) {
        std::cerr << "Error reading zone information.\n";
        cg_close(file_index);
        return false; 
    }

    cgsize_t num_nodes = zone_sizes[0]; 
    cgsize_t num_cells = zone_sizes[1]; // Total 3D elements reported by the zone

    // 1. READ COORDINATES
    mesh.node_coordinate_x.resize(num_nodes); 
    mesh.node_coordinate_y.resize(num_nodes);
    mesh.node_coordinate_z.resize(num_nodes); 

    cgsize_t start_index = 1; 
    cgsize_t end_index = num_nodes; 
    cg_coord_read(file_index, base_id, zone_id, "CoordinateX", CGNS_ENUMV(RealDouble), &start_index, &end_index, mesh.node_coordinate_x.data()); 
    cg_coord_read(file_index, base_id, zone_id, "CoordinateY", CGNS_ENUMV(RealDouble), &start_index, &end_index, mesh.node_coordinate_y.data()); 
    cg_coord_read(file_index, base_id, zone_id, "CoordinateZ", CGNS_ENUMV(RealDouble), &start_index, &end_index, mesh.node_coordinate_z.data()); 

    int num_sections = 0; 
    cg_nsections(file_index, base_id, zone_id, &num_sections); 

    // Setup cell tracking structures
    mesh.cell_type.resize(num_cells); 
    mesh.cell_nodes_offsets.resize(num_cells + 1, 0); 
    mesh.cell_nodes_offsets[0] = 0;

    // Setup 2D boundary element tracking structures
    mesh.boco_element_offsets.push_back(0);

    // Dynamic sizing variables
    int64_t cell_counter = 0;
    int64_t cell_connectivity_offset = 0;
    int64_t boco_element_counter = 0;
    int64_t boco_connectivity_offset = 0;

    bool first_2d_element_seen = false;

    // 2. READ SECTIONS (3D Volumes & 2D Surface Elements)
    for (int sect = 1; sect <= num_sections; ++sect) {
        char section_name[33]; 
        CGNS_ENUMT(ElementType_t) type; 
        cgsize_t ele_start = 0, ele_end = 0; 
        int nbndry = 0, parent_flag = 0; 

        if (cg_section_read(file_index, base_id, zone_id, sect, section_name, &type, &ele_start, &ele_end, &nbndry, &parent_flag) != CG_OK) {
            continue; 
        }

        cgsize_t num_elements = ele_end - ele_start + 1;

        // Determine type classification
        int nodes_per_element = 0; 
        int raw_element_type = -1; 
        bool is_3d_volume = false;
        bool is_2d_surface = false;

        // 3D Volume Cells
        if (type == CGNS_ENUMV(TETRA_4))       { nodes_per_element = 4; raw_element_type = 0; is_3d_volume = true; } 
        else if (type == CGNS_ENUMV(PYRA_5))   { nodes_per_element = 5; raw_element_type = 1; is_3d_volume = true; } 
        else if (type == CGNS_ENUMV(PENTA_6))  { nodes_per_element = 6; raw_element_type = 2; is_3d_volume = true; } 
        else if (type == CGNS_ENUMV(HEXA_8))   { nodes_per_element = 8; raw_element_type = 3; is_3d_volume = true; } 
        // 2D Surface Boundary Elements
        else if (type == CGNS_ENUMV(TRI_3))    { nodes_per_element = 3; raw_element_type = 0; is_2d_surface = true; }
        else if (type == CGNS_ENUMV(QUAD_4))   { nodes_per_element = 4; raw_element_type = 1; is_2d_surface = true; }
        else {
            continue; // Ignore line elements or unsupported shapes safely
        }

        std::vector<cgsize_t> element_connectivity(num_elements * nodes_per_element); 
        cg_elements_read(file_index, base_id, zone_id, sect, element_connectivity.data(), nullptr); 

        if (is_3d_volume) {
            for (cgsize_t e = 0; e < num_elements; ++e) {
                if (cell_counter >= num_cells) break; 

                mesh.cell_type[cell_counter] = raw_element_type; 
                for (int n = 0; n < nodes_per_element; ++n) {
                    mesh.cell_nodes.push_back(static_cast<int64_t>(element_connectivity[e * nodes_per_element + n] - 1)); 
                }

                cell_connectivity_offset += nodes_per_element; 
                mesh.cell_nodes_offsets[cell_counter + 1] = cell_connectivity_offset; 
                cell_counter++; 
            }
        } 
        else if (is_2d_surface) {
            if (!first_2d_element_seen) {
                mesh.boco_elements_global_start = static_cast<int64_t>(ele_start - 1);
                first_2d_element_seen = true;
            }

            for (cgsize_t e = 0; e < num_elements; ++e) {
                mesh.boco_element_type.push_back(raw_element_type);
                for (int n = 0; n < nodes_per_element; ++n) {
                    mesh.boco_element_nodes.push_back(static_cast<int64_t>(element_connectivity[e * nodes_per_element + n] - 1));
                }
                boco_connectivity_offset += nodes_per_element;
                mesh.boco_element_offsets.push_back(boco_connectivity_offset);
                boco_element_counter++;
            }
        }
    }

    // 3. READ ALL BOUNDARY CONDITIONS
    int num_bocos = 0; 
    if (cg_nbocos(file_index, base_id, zone_id, &num_bocos) == CG_OK) {
        mesh.boundary_name.resize(num_bocos); 
        mesh.boundary_type.resize(num_bocos); 
        mesh.boundary_face_start_index.resize(num_bocos); 
        mesh.boundary_face_count.resize(num_bocos); 
        mesh.boundary_explicit_elements.resize(num_bocos);
        mesh.boundary_grid_location.resize(num_bocos);

        std::vector<cgsize_t> pnts; 

        for (int bc = 1; bc <= num_bocos; ++bc) {
            char bc_name[33]; 
            CGNS_ENUMT(BCType_t) bc_type;
            CGNS_ENUMT(PointSetType_t) pt_set_type;
            cgsize_t npts = 0;
            int normal_index = 0; cgsize_t normal_list_size = 0;
            CGNS_ENUMT(DataType_t) normal_data_type; int ndataset = 0;

            if (cg_boco_info(file_index, base_id, zone_id, bc, bc_name, &bc_type, &pt_set_type, &npts, 
                             &normal_index, &normal_list_size, &normal_data_type, &ndataset) != CG_OK) {
                continue;
            }
            
            mesh.boundary_name[bc - 1] = bc_name;
            mesh.boundary_type[bc - 1] = static_cast<int>(bc_type);
            mesh.boundary_face_count[bc - 1] = static_cast<int64_t>(npts);

            pnts.resize(npts); 
            cg_boco_read(file_index, base_id, zone_id, bc, pnts.data(), nullptr);

            CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(FaceCenter); 
            cg_boco_gridlocation_read(file_index, base_id, zone_id, bc, &location);
            mesh.boundary_grid_location[bc - 1] = static_cast<int>(location);

            if (pt_set_type == CGNS_ENUMV(PointRange)) {
                mesh.boundary_face_start_index[bc - 1] = static_cast<int64_t>(pnts[0] - 1); 
                // FIX: Calculate the true count of elements in the range sequence
                mesh.boundary_face_count[bc - 1] = static_cast<int64_t>(pnts[1] - pnts[0] + 1); 
            } 
            else if (pt_set_type == CGNS_ENUMV(PointList)) {
                mesh.boundary_face_start_index[bc - 1] = -1; // Flags layout as list-based
                mesh.boundary_face_count[bc - 1] = static_cast<int64_t>(npts); // npts is correct for lists
                
                mesh.boundary_explicit_elements[bc - 1].reserve(npts);
                for (cgsize_t i = 0; i < npts; ++i) {
                    mesh.boundary_explicit_elements[bc - 1].push_back(static_cast<int64_t>(pnts[i] - 1));
                }
            }
        }
    }

    cg_close(file_index); 
    return true; 
}
