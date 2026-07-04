#include "../Declarations.hpp"


// Compute physical Euler flux vector projected onto a face normal direction
ConservedState compute_physical_flux(const ConservedState& q, const Vector3D& normal) {
    double u = q.rho_u / q.rho; // Extract velocity X component
    double v = q.rho_v / q.rho; // Extract velocity Y component
    double w = q.rho_w / q.rho; // Extract velocity Z component
    
    double velocity_squared = u * u + v * v + w * w; // Magnitude of velocity vector squared
    double pressure = (GAMMA - 1.0) * (q.energy - 0.5 * q.rho * velocity_squared); // Ideal gas law equation
    
    double v_normal = u * normal.x + v * normal.y + w * normal.z; // Project velocity onto face normal
    
    ConservedState flux;
    flux.rho    = q.rho * v_normal; // Continuity flux component
    flux.rho_u  = q.rho_u * v_normal + pressure * normal.x; // X-momentum flux component
    flux.rho_v  = q.rho_v * v_normal + pressure * normal.y; // Y-momentum flux component
    flux.rho_w  = q.rho_w * v_normal + pressure * normal.z; // Z-momentum flux component
    flux.energy = (q.energy + pressure) * v_normal; // Energy flux component
    return flux;
}

// Compute local acoustic speed of sound
double compute_speed_of_sound(const ConservedState& q) {
    double u = q.rho_u / q.rho;
    double v = q.rho_v / q.rho;
    double w = q.rho_w / q.rho;
    double velocity_squared = u * u + v * v + w * w;
    double pressure = (GAMMA - 1.0) * (q.energy - 0.5 * q.rho * velocity_squared);
    return std::sqrt(GAMMA * pressure / q.rho); // Standard thermodynamic speed definition
}

// Handle specialized boundary state based on physical patch type markers
ConservedState ghost_state_generator(const ConservedState& internal, BCType bc) {
    if (bc == BCType::WALL) {
        ConservedState ghost = internal;
        ghost.rho_u = -internal.rho_u; // Reflect X momentum to enforce slip/no-penetration condition
        ghost.rho_v = -internal.rho_v; // Reflect Y momentum to enforce slip/no-penetration condition
        ghost.rho_w = -internal.rho_w; // Reflect Z momentum to enforce slip/no-penetration condition
        return ghost;
    }
    return internal; // Simple extrapolation fallback for transmissive boundaries
}



void SolveExplicitEuler(const Mesh& mesh, FlowFields& fields, double dt) 
{
    
    // 1. Reset all cell residuals before accumulating fluxes
    for (size_t c = 0; c < mesh.cells.size(); ++c) {
        fields.residual[c] = ConservedState{0.0, 0.0, 0.0, 0.0, 0.0}; // Clear past data arrays
    }
    
    // 2. Compute and distribute fluxes across all Internal Faces
    for (size_t f = 0; f < mesh.faces.size(); ++f) {
        const Face& face = mesh.faces[f];
        if (face.is_boundary()) continue; // Skip boundary indices during interior collection loop
        
        int id_l = face.owner_cell_id; // Left cell spatial state mapping index
        int id_r = face.neighbour_cell_id; // Right cell spatial state mapping index
        
        const ConservedState& q_l = fields.state[id_l];
        const ConservedState& q_r = fields.state[id_r];
        
        // Evaluate physical spatial vectors from left and right viewpoints
        ConservedState flux_l = compute_physical_flux(q_l, face.normal);
        ConservedState flux_r = compute_physical_flux(q_r, face.normal);
        
        // Calculate maximum spectral radius for Lax-Friedrichs wave tracking
        double v_n_l = (q_l.rho_u * face.normal.x + q_l.rho_v * face.normal.y + q_l.rho_w * face.normal.z) / q_l.rho;
        double v_n_r = (q_r.rho_u * face.normal.x + q_r.rho_v * face.normal.y + q_r.rho_w * face.normal.z) / q_r.rho;
        double wave_speed_l = std::abs(v_n_l) + compute_speed_of_sound(q_l);
        double wave_speed_r = std::abs(v_n_r) + compute_speed_of_sound(q_r);
        double max_lambda = std::max(wave_speed_l, wave_speed_r); // Pick dominating wave speed limit
        
        // Assemble final Rusanov numerical dissipative framework array values
        ConservedState numerical_flux;
        double area_scale = face.area; // Face scale multiplier tracking metric
        numerical_flux.rho    = 0.5 * (flux_l.rho    + flux_r.rho)    - 0.5 * max_lambda * (q_r.rho    - q_l.rho);
        numerical_flux.rho_u  = 0.5 * (flux_l.rho_u  + flux_r.rho_u)  - 0.5 * max_lambda * (q_r.rho_u  - q_l.rho_u);
        numerical_flux.rho_v  = 0.5 * (flux_l.rho_v  + flux_r.rho_v)  - 0.5 * max_lambda * (q_r.rho_v  - q_l.rho_v);
        numerical_flux.rho_w  = 0.5 * (flux_l.rho_w  + flux_r.rho_w)  - 0.5 * max_lambda * (q_r.rho_w  - q_l.rho_w);
        numerical_flux.energy = 0.5 * (flux_l.energy + flux_r.energy) - 0.5 * max_lambda * (q_r.energy - q_l.energy);
        
        // Distribute calculated surface metric loads symmetrically to shared volumes
        fields.residual[id_l].rho    += numerical_flux.rho    * area_scale; // Sinks fluid from owner
        fields.residual[id_l].rho_u  += numerical_flux.rho_u  * area_scale;
        fields.residual[id_l].rho_v  += numerical_flux.rho_v  * area_scale;
        fields.residual[id_l].rho_w  += numerical_flux.rho_w  * area_scale;
        fields.residual[id_l].energy += numerical_flux.energy * area_scale;
        
        fields.residual[id_r].rho    -= numerical_flux.rho    * area_scale; // Sources fluid into neighbor
        fields.residual[id_r].rho_u  -= numerical_flux.rho_u  * area_scale;
        fields.residual[id_r].rho_v  -= numerical_flux.rho_v  * area_scale;
        fields.residual[id_r].rho_w  -= numerical_flux.rho_w  * area_scale;
        fields.residual[id_r].energy -= numerical_flux.energy * area_scale;
    }
    
    // 3. Process localized boundary loops grouped sequentially by patches
    for (const auto& patch : mesh.boundaries) {
        for (int i = 0; i < patch.face_count; ++i) {
            int face_id = patch.start_face_id + i; // Extract unique face element tracking address
            const Face& b_face = mesh.faces[face_id];
            
            int id_l = b_face.owner_cell_id; // Identify bounding inner physical core region
            const ConservedState& q_l = fields.state[id_l];
            ConservedState q_r = ghost_state_generator(q_l, patch.type); // Create dummy state
            
            ConservedState flux_l = compute_physical_flux(q_l, b_face.normal);
            ConservedState flux_r = compute_physical_flux(q_r, b_face.normal);
            
            double v_n_l = (q_l.rho_u * b_face.normal.x + q_l.rho_v * b_face.normal.y + q_l.rho_w * b_face.normal.z) / q_l.rho;
            double max_lambda = std::abs(v_n_l) + compute_speed_of_sound(q_l); // Simple acoustic scaling speed bound
            
            ConservedState b_flux;
            double area_scale = b_face.area;
            b_flux.rho    = 0.5 * (flux_l.rho    + flux_r.rho)    - 0.5 * max_lambda * (q_r.rho    - q_l.rho);
            b_flux.rho_u  = 0.5 * (flux_l.rho_u  + flux_r.rho_u)  - 0.5 * max_lambda * (q_r.rho_u  - q_l.rho_u);
            b_flux.rho_v  = 0.5 * (flux_l.rho_v  + flux_r.rho_v)  - 0.5 * max_lambda * (q_r.rho_v  - q_l.rho_v);
            b_flux.rho_w  = 0.5 * (flux_l.rho_w  + flux_r.rho_w)  - 0.5 * max_lambda * (q_r.rho_w  - q_l.rho_w);
            b_flux.energy = 0.5 * (flux_l.energy + flux_r.energy) - 0.5 * max_lambda * (q_r.energy - q_l.energy);
            
            fields.residual[id_l].rho    += b_flux.rho    * area_scale; // Accumulate boundary load onto owner
            fields.residual[id_l].rho_u  += b_flux.rho_u  * area_scale;
            fields.residual[id_l].rho_v  += b_flux.rho_v  * area_scale;
            fields.residual[id_l].rho_w  += b_flux.rho_w  * area_scale;
            fields.residual[id_l].energy += b_flux.energy * area_scale;
        }
    }
    
    // 4. Update the solution fields forward in time using Explicit Euler integration
    for (size_t c = 0; c < mesh.cells.size(); ++c) 
    {
        double vol = mesh.cells[c].volume; // Extract static geometric sizing index
        double time_factor = dt / vol; // Compute scaling multiplier scalar index
        
        fields.state[c].rho    -= fields.residual[c].rho    * time_factor; // Update global density field
        fields.state[c].rho_u  -= fields.residual[c].rho_u  * time_factor; // Update global X momentum field
        fields.state[c].rho_v  -= fields.residual[c].rho_v  * time_factor; // Update global Y momentum field
        fields.state[c].rho_w  -= fields.residual[c].rho_w  * time_factor; // Update global Z momentum field
        fields.state[c].energy -= fields.residual[c].energy * time_factor; // Update global total energy field
    }
}
