#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include "../Declarations.hpp" // Pulls in structural types for the solver
#include "../io/cgns/include/cgnslib.h" // CGNS official library C-interface header
#include <iostream> // For tracking initialization logs
#include <cmath> // For processing vector velocity projections
#include <algorithm> // Required for node sorting and matching routines

// Sizes field arrays and populates them with a physical starting state
FlowFields initialize_flow_fields(const Mesh& mesh, double initial_rho, Vector3D initial_vel, double initial_pressure) {
    FlowFields fields;
    
    int number_of_cells = static_cast<int>(mesh.cells.size()); // Count total cells from imported mesh
    
    // 1. Allocate exact memory spaces to prevent runtime resizing
    fields.state.resize(number_of_cells); // Size the state array
    fields.residual.resize(number_of_cells); // Size the parallel residual array
    
    // 2. Compute initial conserved variables from primitive physical inputs
    double vel_sqr = initial_vel.x * initial_vel.x + initial_vel.y * initial_vel.y + initial_vel.z * initial_vel.z; // Velocity squared
    double initial_energy = (initial_pressure / (GAMMA - 1.0)) + (0.5 * initial_rho * vel_sqr); // Total energy calculation
    
    // 3. Uniformly fill every cell field with the starting state
    for (int c = 0; c < number_of_cells; ++c) 
    {
        fields.state[c].rho    = initial_rho; // Set density
        fields.state[c].rho_u  = initial_rho * initial_vel.x; // Set X momentum
        fields.state[c].rho_v  = initial_rho * initial_vel.y; // Set Y momentum
        fields.state[c].rho_w  = initial_rho * initial_vel.z; // Set Z momentum
        fields.state[c].energy = initial_energy; // Set total energy
        
        // Residuals default to zero via structural constructor rules
    }
    
    return fields; // Return ready-to-solve data block
}
