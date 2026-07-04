#include "cgnslib.h"
#include "cgns/include/cgnslib.h"
#include "PreProcessor.h"

bool PreProcessor::CGNSReader(const std::string& filename, Mesh& mesh) {
    int file_index = 0; // File descriptor handle used by the tracking library

    if (cg_open(filename.c_str(), CG_MODE_READ, &file_index) != CG_OK) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return false; // Interrupt execution if the target file cannot be reached
    }

    int base_id = 1; // Default assumption for standard single-base database layouts
    int zone_id = 1; // Default assumption for standalone single-zone grid topologies
    char zone_name[33]; // Character buffer to hold the descriptive zone name string
    cgsize_t zone_sizes[3] = {0, 0, 0}; // Sizing array to capture node, cell, and boundary counts

    if (cg_zone_read(file_index, base_id, zone_id, zone_name, zone_sizes) != CG_OK) {
        std::cerr << "Error reading zone information.\n";
        cg_close(file_index);
        return false; // Close out file tracking descriptors to avoid system leaks before leaving
    }

    cgsize_t num_nodes = zone_sizes[0]; // Total vertex data count recorded inside the file block
    cgsize_t num_cells = zone_sizes[1]; // Total volumetric element tracking count reported by the zone

    mesh.node_coordinate_x.resize(num_nodes); // Shape internal storage layout arrays to match vertex limits
    mesh.node_coordinate_y.resize(num_nodes); // Adjust coordinate space footprint sizing arrays
    mesh.node_coordinate_z.resize(num_nodes); // Match structural limits for vertical coordinate dimensions

    cgsize_t start_index = 1; // Set reading pointer cursor to base-1 index start rule
    cgsize_t end_index = num_nodes; // Set tracking pointer upper limit boundary to terminal index marker
    cg_coord_read(file_index, base_id, zone_id, "CoordinateX", CGNS_ENUMV(RealDouble), &start_index, &end_index, mesh.node_coordinate_x.data()); // Load continuous raw structural x data
    cg_coord_read(file_index, base_id, zone_id, "CoordinateY", CGNS_ENUMV(RealDouble), &start_index, &end_index, mesh.node_coordinate_y.data()); // Load continuous raw structural y data
    cg_coord_read(file_index, base_id, zone_id, "CoordinateZ", CGNS_ENUMV(RealDouble), &start_index, &end_index, mesh.node_coordinate_z.data()); // Load continuous raw structural z data

    int num_sections = 0; // Counter allocation used to inventory topological chunk blocks
    cg_nsections(file_index, base_id, zone_id, &num_sections); // Extract total discrete elements section blocks count

    mesh.cell_type.resize(num_cells); // Size the structural tracking array for layout identifier flags
    mesh.cell_nodes_offsets.resize(num_cells + 1); // Expand the compression offset index tracking boundary vector
    mesh.cell_nodes_offsets[0] = 0; // Initialize root starting tracking marker positions directly to zero
    mesh.cell_nodes.reserve(num_cells * 8); // Pre-allocate storage matrix maps to comfortably mitigate reallocations overhead

    cgsize_t cell_counter = 0; // Running structural index tracks tracking total cells parsed
    int current_offset = 0; // Accumulator tracker keeping track of global connectivity storage lookups

    for (int sect = 1; sect <= num_sections; ++sect) {
        char section_name[33]; // String character buffer used to capture local block name strings
        CGNS_ENUMT(ElementType_t) type; // Enumeration tracker field to intercept element shape keys
        cgsize_t ele_start = 0, ele_end = 0; // Tracking variables indicating localized continuous entry indices bounds
        int nbndry = 0, parent_flag = 0; // Auxiliary tracking indicators mapped out by the core library formatting

        if (cg_section_read(file_index, base_id, zone_id, sect, section_name, &type, &ele_start, &ele_end, &nbndry, &parent_flag) != CG_OK) {
            continue; // Skip past any corrupt data headers encountered within the section sequence
        }

        int nodes_per_element = 0; // Sizing footprint map to scale cell parsing loops
        int raw_cell_type = -1; // Local structural code indicator mapping system flags

        if (type == CGNS_ENUMV(TETRA_4)) {
            nodes_per_element = 4;
            raw_cell_type = 0; // Map tetrahedral configurations down to system code zero
        } else if (type == CGNS_ENUMV(PYRA_5)) {
            nodes_per_element = 5;
            raw_cell_type = 1; // Map pyramidal configurations down to system code one
        } else if (type == CGNS_ENUMV(PENTA_6)) {
            nodes_per_element = 6;
            raw_cell_type = 2; // Map prismatic configurations down to system code two
        } else if (type == CGNS_ENUMV(HEXA_8)) {
            nodes_per_element = 8;
            raw_cell_type = 3; // Map hexahedral configurations down to system code three
        } else {
            continue; // Safely bypass surface elements or non-supported polyhedral boundaries sections
        }

        cgsize_t num_elements = ele_end - ele_start + 1; // Compute absolute volume elements contained in current segment block
        std::vector<cgsize_t> element_connectivity(num_elements * nodes_per_element); // Allocate transient collection array space
        cg_elements_read(file_index, base_id, zone_id, sect, element_connectivity.data(), nullptr); // Pull raw connectivity maps flat into memory space

        for (cgsize_t e = 0; e < num_elements; ++e) {
            if (cell_counter >= num_cells) {
                break; // Halted if processing index maps breach predefined array sizes boundaries
            }

            mesh.cell_type[cell_counter] = raw_cell_type; // Log the layout mapping key straight to data arrays
            for (int n = 0; n < nodes_per_element; ++n) {
                mesh.cell_nodes.push_back(static_cast<int>(element_connectivity[e * nodes_per_element + n] - 1)); // Downshift 1-based indexing structure elements safely to 0-based layout arrays
            }

            current_offset += nodes_per_element; // Shift total item tracking counts cursor forward across space arrays
            mesh.cell_nodes_offsets[cell_counter + 1] = current_offset; // Commit index break locations markers to tracking structure maps
            cell_counter++; // Process next tracking index position safely
        }
    }

    int num_bocos = 0; // Variable allocation to detect boundary configuration records blocks
    if (cg_nbocos(file_index, base_id, zone_id, &num_bocos) == CG_OK) {
        mesh.boundary_name.resize(num_bocos); // Sync structure arrays sizing for boundary block string layouts
        mesh.boundary_type.resize(num_bocos); // Sync sizing configurations maps trackers for identification identifiers
        mesh.boundary_face_start_index.resize(num_bocos); // Scale structural coordinate markers lookup indexes arrays
        mesh.boundary_face_count.resize(num_bocos); // Scale allocation arrays detailing total boundaries group limits

        for (int bc = 1; bc <= num_bocos; ++bc) {
            char bc_name[33]; // Character storage array tracking patch descriptor identifiers names
            CGNS_ENUMT(BCType_t) bc_type; // Enumeration type flag mapping group behaviors properties
            CGNS_ENUMT(PointSetType_t) pt_set_type; // Format descriptor tracking structural indexing maps layouts
            cgsize_t npts = 0; // Total nodes or face entities bound to this target patch identifier group
            int normal_index = 0; // Structural tracking cursor used for explicit coordinate properties
            cgsize_t normal_list_size = 0; // Track variable evaluating raw normal vector streams size bounds
            CGNS_ENUMT(DataType_t) normal_data_type; // Structural enumeration tracing precise scalar floating precisions markers
            int ndataset = 0; // Count tracking internal custom property structures attached here

            cg_boco_info(file_index, base_id, zone_id, bc, bc_name, &bc_type, &pt_set_type, &npts, &normal_index, &normal_list_size, &normal_data_type, &ndataset); // Query core metadata details for the patch structure records
            mesh.boundary_name[bc - 1] = std::string(bc_name); // Convert and assign text attributes straight to mesh storage space properties
            mesh.boundary_type[bc - 1] = static_cast<int>(bc_type); // Log condition properties directly via tracking indicator scalar variables

            std::vector<cgsize_t> pnts(npts); // Temporary vector allocated to intercept entity tracking ID numbers arrays
            cg_boco_read(file_index, base_id, zone_id, bc, pnts.data(), nullptr); // Pull referencing identity arrays from file storage blocks directly to memory space
            if (pt_set_type == CGNS_ENUMV(PointRange)) {
                mesh.boundary_face_start_index[bc - 1] = static_cast<int>(pnts[0]); // Save initial start position index details directly
                mesh.boundary_face_count[bc - 1] = static_cast<int>(pnts[1] - pnts[0] + 1); // Calculate span parameters for continuous linear sequences loops
            } else if (pt_set_type == CGNS_ENUMV(PointList)) {
                mesh.boundary_face_start_index[bc - 1] = -1; // Flag layout to indicate scattered non-linear listing layouts array maps
                mesh.boundary_face_count[bc - 1] = static_cast<int>(npts); // Track total count bounds directly from elements parsed fields
            }
        }
    }

    cg_close(file_index); // Destroy system file locks and stream objects safely to free lingering memory allocations tracking
    return true; // Execution sequence completed successfully
}