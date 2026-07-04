#include "PreProcessor.h"
#include <vector>
#include <cmath>
#include <iostream>

void PreProcessor::ComputeProperties(Mesh& mesh) {
    int num_cells = static_cast<int>(mesh.cell_type.size());

    // 1. Allocate space for cell properties
    mesh.cell_volume.assign(num_cells, 0.0);
    mesh.cell_center_x.assign(num_cells, 0.0);
    mesh.cell_center_y.assign(num_cells, 0.0);
    mesh.cell_center_z.assign(num_cells, 0.0);

    // 2. Loop through every cell to compute volume and centers using the Divergence Theorem
    for (int c = 0; c < num_cells; ++c) {
        int face_start = mesh.cell_faces_offsets[c];
        int face_end = mesh.cell_faces_offsets[c + 1];

        double vol_accumulator = 0.0;
        double cx_accumulator = 0.0;
        double cy_accumulator = 0.0;
        double cz_accumulator = 0.0;

        for (int i = face_start; i < face_end; ++i) {
            int f = mesh.cell_faces[i];

            // Extract face geometry
            double af = mesh.face_area[f];
            double fx = mesh.face_center_x[f];
            double fy = mesh.face_center_y[f];
            double fz = mesh.face_center_z[f];
            
            double nx = mesh.face_normal_x[f];
            double ny = mesh.face_normal_y[f];
            double nz = mesh.face_normal_z[f];

            // Check orientation: The normal must point OUTWARD from the current cell.
            // By convention, face_normals point from owner to neighbor.
            // If this cell is the neighbor, we flip the normal direction sign.
            double sign = 1.0;
            if (mesh.face_neighbour_cell_index[f] == c) {
                sign = -1.0;
            }

            double oriented_nx = nx * sign;
            double oriented_ny = ny * sign;
            double oriented_nz = nz * sign;

            // Dot product: x_face_center · n_oriented
            double dot_product = (fx * oriented_nx) + (fy * oriented_ny) + (fz * oriented_nz);

            // Accumulate Volumetric Component
            // V = 1/3 * sum( (xf · n) * Area )
            double face_vol_contribution = dot_product * af;
            vol_accumulator += face_vol_contribution;

            // Accumulate Centroid Components
            // x_c = 1/(4V) * sum( xf * (xf · n) * Area )
            cx_accumulator += fx * face_vol_contribution;
            cy_accumulator += fy * face_vol_contribution;
            cz_accumulator += fz * face_vol_contribution;
        }

        // Finalize cell calculations
        double total_volume = vol_accumulator / 3.0;

        if (total_volume > 1e-15) {
            mesh.cell_volume[c] = total_volume;
            mesh.cell_center_x[c] = cx_accumulator / (4.0 * total_volume);
            mesh.cell_center_y[c] = cy_accumulator / (4.0 * total_volume);
            mesh.cell_center_z[c] = cz_accumulator / (4.0 * total_volume);
        } else {
            std::cerr << "Warning: Cell " << c << " has an invalid or zero volume (" << total_volume << ").\n";
            mesh.cell_volume[c] = 0.0;
            mesh.cell_center_x[c] = 0.0;
            mesh.cell_center_y[c] = 0.0;
            mesh.cell_center_z[c] = 0.0;
        }
    }

}