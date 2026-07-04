#pragma once
#include <vector>

struct FlowField 
{
    // --- 1. CONSERVATIVE FIELDS (Stored at Cell Centers) ---
    std::vector<double> rho;               // Density (Mass)
    std::vector<double> rhou;              // X-Direction Momentum 
    std::vector<double> rhov;              // Y-Direction Momentum
    std::vector<double> rhow;              // Z-Direction Momentum
    std::vector<double> rhoE;              // Total Energy

    // --- 2. RESIDUAL ACCUMULATORS (Net Fluxes per Cell) ---
    std::vector<double> rho_residual;     
    std::vector<double> rhou_residual;     // Net change tracking for X-Momentum
    std::vector<double> rhov_residual;     // Net change tracking for Y-Momentum
    std::vector<double> rhow_residual;     // Net change tracking for Z-Momentum
    std::vector<double> rhoE_residual;     // Net change tracking for Energy

    // Helper function to size everything to match your mesh cell count
    void allocate(size_t num_cells) 
    {
        rho.resize(num_cells, 0.0);
        rhou.resize(num_cells, 0.0);
        rhov.resize(num_cells, 0.0);
        rhow.resize(num_cells, 0.0);
        rhoE.resize(num_cells, 0.0);
        
        rho_residual.resize(num_cells, 0.0);
        rhou_residual.resize(num_cells, 0.0);
        rhov_residual.resize(num_cells, 0.0);
        rhow_residual.resize(num_cells, 0.0);
        rhoE_residual.resize(num_cells, 0.0);
    }
};