#ifndef IMPORT_MESH_HPP
#define IMPORT_MESH_HPP

#include "../Declarations.hpp" // Pulls in structural types for the mesh
#include <iostream> 
#include <vector>
#include <string>
#include "cgns/include/cgnslib.h" // CGNS official library C-interface header

enum class CellType {
    UNKNOWN = 0,
    TETRA_4 = 1,
    HEXA_8 = 2,
    PRISM_6 = 3,
    PYRA_5 = 4
};

struct ElementSection {
    std::string name;
    size_t start_index;
    size_t end_index;
    size_t element_count;
    CellType type;
    std::vector<size_t> connectivity;
};

struct BoundaryCondition {
    std::string name;
    std::string type;
    std::vector<size_t> element_ids;
};

struct MeshInfo {
    size_t num_nodes;
    size_t num_cells;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<ElementSection> sections;
    std::vector<BoundaryCondition> boundaries;
};

inline bool Acquire_Mesh(const std::string& filename, RawMeshData& mesh_out) 
{
    int file_id = 0;
    int base_id = 1; 
    char base_name[33];
    int zone_id = 1;
    char zone_name[33]; 
    int cell_dim = 0, phys_dim = 0;
    cgsize_t zone_sizes[3] = {0, 0, 0};

    // FIXED: Must open the file first to get a valid file_id
    if (cg_open(filename.c_str(), CG_MODE_READ, &file_id) != CG_OK) {
        std::cerr << "Error: Could not open CGNS file: " << filename << std::endl;
        return false;
    }

    if ((cg_base_read(file_id, base_id, base_name, &cell_dim, &phys_dim) != CG_OK) || 
        (phys_dim != 3) || 
        (cg_zone_read(file_id, base_id, zone_id, zone_name, zone_sizes) != CG_OK))
    {
        cg_close(file_id); 
        return false;
    }

    mesh_out.num_nodes = static_cast<size_t>(zone_sizes[0]);
    mesh_out.num_cells = static_cast<size_t>(zone_sizes[1]);

    mesh_out.x.resize(mesh_out.num_nodes);
    mesh_out.y.resize(mesh_out.num_nodes);
    mesh_out.z.resize(mesh_out.num_nodes);

    cgsize_t r_min = 1; 
    cgsize_t r_max = static_cast<cgsize_t>(mesh_out.num_nodes);
    cg_coord_read(file_id, base_id, zone_id, "CoordinateX", CGNS_ENUMV(RealDouble), &r_min, &r_max, mesh_out.x.data());
    cg_coord_read(file_id, base_id, zone_id, "CoordinateY", CGNS_ENUMV(RealDouble), &r_min, &r_max, mesh_out.y.data());
    cg_coord_read(file_id, base_id, zone_id, "CoordinateZ", CGNS_ENUMV(RealDouble), &r_min, &r_max, mesh_out.z.data());

    int num_sections = 0; 
    cg_nsections(file_id, base_id, zone_id, &num_sections);

    for (int s = 1; s <= num_sections; ++s) {
        char sect_name[33]; 
        CGNS_ENUMT(ElementType_t) el_type = CGNS_ENUMV(ElementTypeNull);
        cgsize_t start = 0, end = 0; 
        int nbndry = 0, parent_flag = 0;
        if (cg_section_read(file_id, base_id, zone_id, s, sect_name, &el_type, &start, &end, &nbndry, &parent_flag) != CG_OK) continue;

        ElementSection section;
        section.name = std::string(sect_name);
        section.start_index = static_cast<size_t>(start);
        section.end_index = static_cast<size_t>(end);
        section.element_count = static_cast<size_t>(end - start + 1);

        if (el_type == CGNS_ENUMV(TETRA_4))      section.type = CellType::TETRA_4;
        else if (el_type == CGNS_ENUMV(HEXA_8))  section.type = CellType::HEXA_8;
        else if (el_type == CGNS_ENUMV(PENTA_6)) section.type = CellType::PRISM_6;
        else if (el_type == CGNS_ENUMV(PYRA_5))  section.type = CellType::PYRA_5;
        else section.type = CellType::UNKNOWN; 

        // Skip 2D element boundary section blocks from populating volumetric arrays
        if (section.type == CellType::UNKNOWN) continue;

        cgsize_t data_size = 0; cg_ElementDataSize(file_id, base_id, zone_id, s, &data_size);
        std::vector<cgsize_t> raw_conn(data_size);
        cg_elements_read(file_id, base_id, zone_id, s, raw_conn.data(), nullptr);

        section.connectivity.reserve(data_size);
        for (auto val : raw_conn) section.connectivity.push_back(static_cast<size_t>(val - 1));
        mesh_out.sections.push_back(section);
    }

    // Extracted BC parsing to a dedicated separate block to be clean
    cg_close(file_id); 
    return true; 
}

// Minimal, unified function called AFTER pre.generateFaces() inside your main flow pipeline
inline void ApplyBoundaryConditions(const std::string& filename, Mesh& mesh) {
    int file_id = 0, base_id = 1, zone_id = 1, num_bcs = 0;
    if (cg_open(filename.c_str(), CG_MODE_READ, &file_id) != CG_OK) return;
    
    cg_nbocos(file_id, base_id, zone_id, &num_bcs);
    for (int b = 1; b <= num_bcs; ++b) {
        char bc_name[33]; CGNS_ENUMT(BCType_t) bctype; CGNS_ENUMT(PointSetType_t) ptset_type;
        cgsize_t npts = 0; int norm_idx = 0; cgsize_t norm_sz = 0;
        CGNS_ENUMT(DataType_t) norm_dt; int ndataset = 0;

        cg_boco_info(file_id, base_id, zone_id, b, bc_name, &bctype, &ptset_type, &npts, &norm_idx, &norm_sz, &norm_dt, &ndataset);

        int internal_boco_id = 1; 
        std::string bc_str(cg_BCTypeName(bctype));
        if (bc_str == "BCWallInviscid" || bc_str == "BCWall") internal_boco_id = 2;
        else if (bc_str == "BCFarfield") internal_boco_id = 3;

        std::vector<cgsize_t> bc_elements(npts);
        cg_boco_read(file_id, base_id, zone_id, b, bc_elements.data(), nullptr);

        for (auto p : bc_elements) {
            size_t face_index = static_cast<size_t>(p - 1);
            if (face_index < mesh.faces.size()) {
                mesh.faces[face_index].boco_marker = internal_boco_id;
            }
        }
    }
    cg_close(file_id);
}

inline void SolverIO::importCGNS(const std::string& filename, Mesh& mesh) 
{
    RawMeshData raw_data;
    if (!Acquire_Mesh(filename, raw_data)) return;

    mesh.nodes.resize(raw_data.num_nodes);
    for (size_t i = 0; i < raw_data.num_nodes; ++i) 
    {
        mesh.nodes[i].id = static_cast<int>(i);
        mesh.nodes[i].coord.x = raw_data.x[i];
        mesh.nodes[i].coord.y = raw_data.y[i];
        mesh.nodes[i].coord.z = raw_data.z[i];
    }

    size_t calculated_cell_index = 0;
    for (const auto& sect : raw_data.sections)
    {
        size_t nodes_per_elem = sect.connectivity.size() / sect.element_count;

        for (size_t e = 0; e < sect.element_count; ++e) 
        {
            Cell cell_instance;
            cell_instance.id = static_cast<int>(calculated_cell_index++);
            cell_instance.type = static_cast<int>(sect.type);
            
            for (size_t n = 0; n < nodes_per_elem; ++n) 
            {
                size_t node_ptr = sect.connectivity[e * nodes_per_elem + n];
                cell_instance.node_ids.push_back(static_cast<int>(node_ptr));
            }
            mesh.cells.push_back(cell_instance);
        }
    }
}

#endif