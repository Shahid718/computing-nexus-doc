//----------------------------------------------------------------------------------------------------
//  Cahn-Hilliard Phase-Field Solver
//  High-Performance Finite Difference Implementation with OpenMP Parallelization
//
//  Solves: ∂c/∂t = M∇²(∂f/∂c - κ∇²c)  |  Domain: Nx × Ny  |  BC: Periodic
//
//  (c) 2026 Shahid Maqbool
//  Built with: C++17, OpenMP 4.5+, GCC/Clang
//  Usage: OMP_NUM_THREADS=N ./cahn_hilliard
//----------------------------------------------------------------------------------------------------

#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <omp.h>

// Use modern C++ features without pulling everything into global namespace
using namespace std::chrono;

//--- simulation cell parameters
constexpr int Nx = 2000;
constexpr int Ny = 2000; //128;
constexpr double dx = 1.0;
constexpr double dy = 1.0;

//--- time integration parameters
constexpr int nsteps = 1000; //2000;
constexpr int nprint = 1000;
constexpr double dt = 0.01;

//--- material specific parameters
constexpr double con_0 = 0.4;
constexpr double mobility = 1.0;
constexpr double grad_coef = 0.5;

//--- microstructure parameters
constexpr double noise = 0.02;
constexpr double A = 1.0;

int ip,im,jp,jm;

int main() 
{
    // ✅ ALLOCATE ARRAYS 
    Matrix con(Nx, std::vector<double>(Ny, 0.0));
    Matrix con_new(Nx, std::vector<double>(Ny, 0.0));
    Matrix dfdcon(Nx, std::vector<double>(Ny, 0.0));
    Matrix lap_con(Nx, std::vector<double>(Ny, 0.0));
    Matrix dummy_con(Nx, std::vector<double>(Ny, 0.0));
    Matrix lap_dummy(Nx, std::vector<double>(Ny, 0.0));
    
    // Start timing
    auto start_time = high_resolution_clock::now();

    //--- initial microstructure
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double r = dist(rng);
            con[i][j] = con_0 + noise * (0.5 - r);
        }
    }

    //--- start microstructure evolution
    for (int tsteps = 1; tsteps <= nsteps; tsteps++)
    {
        // First parallel loop: Compute dfdcon, lap_con, and dummy_con
        #pragma omp parallel for collapse(2) \
        shared(con, dfdcon, lap_con, dummy_con) \
        private(ip, im, jp, jm)
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                //--- free energy derivative
                dfdcon[i][j] = A * (2.0 * con[i][j] * (1.0 - con[i][j]) * (1.0 - con[i][j]) - 
                     2.0 * con[i][j] * con[i][j] * (1.0 - con[i][j]));

                //--- laplace evaluation for con
                int jp = j + 1;
                int jm = j - 1;
                int ip = i + 1;
                int im = i - 1;

                if (im == -1) im = Nx - 1;
                if (ip == Nx) ip = 0;
                if (jm == -1) jm = Ny - 1;
                if (jp == Ny) jp = 0;

                lap_con[i][j] = (con[ip][j] + con[im][j] + con[i][jm] +
                                 con[i][jp] - 4.0 * con[i][j]) / (dx * dx);

                dummy_con[i][j] = dfdcon[i][j] - grad_coef * lap_con[i][j];
            }
        }
        
        // Second parallel loop: Compute lap_dummy and update con_new
        #pragma omp parallel for collapse(2) \
        shared(con, con_new, dummy_con, lap_dummy, dt, mobility) \
        private(ip, im, jp, jm)        
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {                
                //--- laplace evaluation for dummy_con
                int jp = j + 1;
                int jm = j - 1;
                int ip = i + 1;
                int im = i - 1;

                if (im == -1) im = Nx - 1;
                if (ip == Nx) ip = 0;
                if (jm == -1) jm = Ny - 1;
                if (jp == Ny) jp = 0;
                
                // Compute Laplacian of dummy_con
                lap_dummy[i][j] = (dummy_con[ip][j] + dummy_con[im][j] + 
                                   dummy_con[i][jm] + dummy_con[i][jp] - 
                                   4.0 * dummy_con[i][j]) / (dx * dx);
                
                // Time integration
                con_new[i][j] = con[i][j] + dt * mobility * lap_dummy[i][j];
                con_new[i][j] = std::clamp(con_new[i][j], 0.00001, 0.99999);
            }
        }
        
        // Update con
        std::swap(con, con_new);

        //--- print steps
        if (tsteps % nprint == 0) {
            std::cout << "Done steps = " << tsteps << std::endl;
        }
    }
    
    //--- End timing
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
       
    //--- Print timing results
    std::cout << "\n---------------------------------" << std::endl;
    std::cout << "  Time       = " << duration.count() / 1000.0 << " seconds." << std::endl;
    std::cout << "---------------------------------" << std::endl;
    
    //--- Write results to file using our helper function
    save_matrix_to_file(con, "ch.dat", 6);
    
    return EXIT_SUCCESS;
}