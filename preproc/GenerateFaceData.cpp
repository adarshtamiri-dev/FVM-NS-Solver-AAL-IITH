#include "PreProcessor.h"
#include <vector>
#include <cmath>
#include <iostream>

void PreProcessor::ComputeProperties(Mesh& mesh) {
    int64_t num_cells = static_cast<int64_t>(mesh.cell_type.size());

    // 1. Allocate space using modern 64-bit index limits matching our updated SoA Mesh
    mesh.cell_volume.assign(num_cells, 0.0);
    mesh.cell_center_x.assign(num_cells, 0.0);
    mesh.cell_center_y.assign(num_cells, 0.0);
    mesh.cell_center_z.assign(num_cells, 0.0);

    // 2. Loop through cells using sub-pyramid decomposition (Highly robust for all unstructured shapes)
    for (int64_t c = 0; c < num_cells; ++c) {
        int64_t face_start = mesh.cell_faces_offsets[c];
        int64_t face_end = mesh.cell_faces_offsets[c + 1];

        // Step A: Calculate a temporary local anchor point to use as a virtual apex.
        // The arithmetic mean of face centers works perfectly and keeps numbers stable.
        double ax = 0.0, ay = 0.0, az = 0.0;
        int64_t num_cell_faces = face_end - face_start;
        
        for (int64_t i = face_start; i < face_end; ++i) {
            int64_t f = mesh.cell_faces[i];
            ax += mesh.face_center_x[f];
            ay += mesh.face_center_y[f];
            az += mesh.face_center_z[f];
        }
        if (num_cell_faces > 0) {
            ax /= num_cell_faces; ay /= num_cell_faces; az /= num_cell_faces;
        }

        double total_volume = 0.0;
        double cx_weighted = 0.0;
        double cy_weighted = 0.0;
        double cz_weighted = 0.0;

        // Step B: Integrate sub-pyramids formed by each face and the anchor apex
        for (int64_t i = face_start; i < face_end; ++i) {
            int64_t f = mesh.cell_faces[i];

            double af = mesh.face_area[f];
            double fx = mesh.face_center_x[f];
            double fy = mesh.face_center_y[f];
            double fz = mesh.face_center_z[f];
            
            double nx = mesh.face_normal_x[f];
            double ny = mesh.face_normal_y[f];
            double nz = mesh.face_normal_z[f];

            // Orient normal outward relative to the cell container boundary
            double sign = 1.0;
            if (mesh.face_neighbour_cell_index[f] == c) {
                sign = -1.0;
            }

            double oriented_nx = nx * sign;
            double oriented_ny = ny * sign;
            double oriented_nz = nz * sign;

            // Height of pyramid = projection of vector (face_center - anchor) onto oriented face normal
            double h = ((fx - ax) * oriented_nx) + 
                       ((fy - ay) * oriented_ny) + 
                       ((fz - az) * oriented_nz);

            // Volume of a pyramid = (Base Area * Height) / 3
            double pyramid_vol = (af * h) / 3.0;

            // Centroid of a pyramid is located exactly 3/4 of the way from the apex to the face base center
            double pyr_cx = ax + 0.75 * (fx - ax);
            double pyr_cy = ay + 0.75 * (fy - ay);
            double pyr_cz = az + 0.75 * (fz - az);

            // Accumulate volume-weighted properties
            total_volume += pyramid_vol;
            cx_weighted  += pyr_cx * pyramid_vol;
            cy_weighted  += pyr_cy * pyramid_vol;
            cz_weighted  += pyr_cz * pyramid_vol;
        }

        // Finalize cell data metrics
        if (total_volume > 1e-15) {
            mesh.cell_volume[c]   = total_volume;
            mesh.cell_center_x[c] = cx_weighted / total_volume;
            mesh.cell_center_y[c] = cy_weighted / total_volume;
            mesh.cell_center_z[c] = cz_weighted / total_volume;
        } else {
            std::cerr << "Warning: Cell " << c << " has an invalid or near-zero volume (" << total_volume << ").\n";
            mesh.cell_volume[c]   = 0.0;
            mesh.cell_center_x[c] = ax; // Fallback to anchor position
            mesh.cell_center_y[c] = ay;
            mesh.cell_center_z[c] = az;
        }
    }
}
