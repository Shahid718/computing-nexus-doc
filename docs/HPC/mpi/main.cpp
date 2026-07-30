//-------------------------------------------------------------------------------
//   MPI Finite Difference Phase Field Code of Cahn-Hilliard Equation (Eigen Optimized)
//   Compile: mpic++ main.cpp -O3 -march=native -std=c++20 -I /usr/include/eigen3/ -o main
//   run    : mpirun -np 4 ./main
//   			SYNCRHONIZED
//
//   Compile:  mpiicpx -I"C:\libs\eigen-5.0.0" main.cpp
//   run    :  mpiexec -n 4 a.exe
//
//-------------------------------------------------------------------------------

#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <mpi.h>
#include <Eigen/Dense>
#include <fstream>
#include <format>
#include <ranges>
#include <span>
#include <iomanip>


using namespace std::chrono;

//--- Scaled simulation cell parameters
constexpr int Nx = 5000;          // Reduce for testing
constexpr int Ny = 5000;
constexpr float dx = 1.0f;
constexpr float dy = 1.0f;

constexpr int nsteps = 1000;
constexpr int nprint = 1000;
constexpr float dt = 0.01f;  

constexpr float con_0 = 0.4f;
constexpr float mobility = 1.0f;
constexpr float grad_coef = 0.5f;
constexpr float noise = 0.02f;
constexpr float A = 1.0f;

// Clamping constants
constexpr float CLAMP_LOW = 1e-5f;
constexpr float CLAMP_HIGH = 0.99999f;

// Use RowMajor Array to easily pass rows across MPI
using EigenMatrix = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Helper function for periodic y-boundary conditions
inline void apply_periodic_y(EigenMatrix& field, int local_Nx) {
    // Apply periodic boundary conditions in y-direction
    for (int i = 1; i <= local_Nx; ++i) {
        field(i, 0) = field(i, Ny - 2);
        field(i, Ny - 1) = field(i, 1);
    }
}

// Helper function for MPI halo exchange in x-direction
inline void exchange_x_halos(EigenMatrix& field, int local_Nx, int left, int right, 
                             MPI_Comm comm, MPI_Status& status) {
    // Send row 1 to left neighbor, receive into row local_Nx + 1 from right
    MPI_Sendrecv(field.row(1).data(), Ny, MPI_FLOAT, left, 0,
                 field.row(local_Nx + 1).data(), Ny, MPI_FLOAT, right, 0,
                 comm, &status);

    // Send row local_Nx to right neighbor, receive into row 0 from left
    MPI_Sendrecv(field.row(local_Nx).data(), Ny, MPI_FLOAT, right, 1,
                 field.row(0).data(), Ny, MPI_FLOAT, left, 1,
                 comm, &status);
}

// Optimized stencil computation
inline void compute_laplacian_and_chemical_potential(
    const EigenMatrix& con, 
    EigenMatrix& dfdcon,
    EigenMatrix& lap_con,
    EigenMatrix& dummy_con,
    int local_Nx,
    float dx, float dy, float A, float grad_coef
) {
    const float inv_dxdy = 1.0f / (dx * dy);
    
    // Precompute neighbor indices for y-direction (periodic)
    static std::vector<int> jp(Ny), jm(Ny);
    static bool initialized = false;
    if (!initialized) {
        for (int j = 0; j < Ny; ++j) {
            jp[j] = (j + 1) % Ny;
            jm[j] = (j - 1 + Ny) % Ny;
        }
        initialized = true;
    }
    
    for (int i = 1; i <= local_Nx; ++i) {
        const auto& con_i = con.row(i);
        const auto& con_ip1 = con.row(i + 1);
        const auto& con_im1 = con.row(i - 1);
        auto dfdcon_i = dfdcon.row(i);
        auto lap_con_i = lap_con.row(i);
        auto dummy_con_i = dummy_con.row(i);
        
        for (int j = 0; j < Ny; ++j) {
            const float c = con_i(j);
            
            // Correct chemical potential derivative
            dfdcon_i(j) = A * (2.0f * c * (1.0f - c) * (1.0f - 2.0f * c));
            
            // Laplacian with periodic y-boundary
            lap_con_i(j) = (con_ip1(j) + con_im1(j) + 
                           con_i(jp[j]) + con_i(jm[j]) - 
                           4.0f * c) * inv_dxdy;
            
            dummy_con_i(j) = dfdcon_i(j) - grad_coef * lap_con_i(j);
        }
    }
}

// Optimized second Laplacian and update
inline void update_field(
    const EigenMatrix& con,
    const EigenMatrix& dummy_con,
    EigenMatrix& con_new,
    int local_Nx,
    float dt, float mobility, float inv_dxdy
) {
    static std::vector<int> jp(Ny), jm(Ny);
    static bool initialized = false;
    if (!initialized) {
        for (int j = 0; j < Ny; ++j) {
            jp[j] = (j + 1) % Ny;
            jm[j] = (j - 1 + Ny) % Ny;
        }
        initialized = true;
    }
    
    const float factor = dt * mobility * inv_dxdy;
    
    for (int i = 1; i <= local_Nx; ++i) {
        const auto& con_i = con.row(i);
        const auto& dummy_i = dummy_con.row(i);
        const auto& dummy_ip1 = dummy_con.row(i + 1);
        const auto& dummy_im1 = dummy_con.row(i - 1);
        auto con_new_i = con_new.row(i);
        
        for (int j = 0; j < Ny; ++j) {
            const float lap_dummy = (dummy_ip1(j) + dummy_im1(j) + 
                                    dummy_i(jp[j]) + dummy_i(jm[j]) - 
                                    4.0f * dummy_i(j));
            
            float c_new = con_i(j) + factor * lap_dummy;
            con_new_i(j) = std::max(CLAMP_LOW, std::min(c_new, CLAMP_HIGH));
        }
    }
}

int main(int argc, char** argv) 
{
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        std::cout << "Total processors running: " << size << std::endl;
        std::cout << "Grid size: " << Nx << " x " << Ny << std::endl;
        std::cout << "Using float precision" << std::endl;
    }
    
    // Domain decomposition along x-direction
    int base_Nx = Nx / size;
    int remainder = Nx % size;
    
    int local_Nx = base_Nx + (rank < remainder ? 1 : 0);
    int start_x = rank * base_Nx + std::min(rank, remainder);
    int end_x = start_x + local_Nx;
    
    int local_Nx_halo = local_Nx + 2;
    
    float mpi_start_time = MPI_Wtime();
    
    // Allocate continuous Eigen blocks
    EigenMatrix con = EigenMatrix::Zero(local_Nx_halo, Ny);
    EigenMatrix dfdcon = EigenMatrix::Zero(local_Nx_halo, Ny);
    EigenMatrix lap_con = EigenMatrix::Zero(local_Nx_halo, Ny);
    EigenMatrix dummy_con = EigenMatrix::Zero(local_Nx_halo, Ny);
    EigenMatrix con_new = EigenMatrix::Zero(local_Nx_halo, Ny);
    
    // Initialize with noise
    std::mt19937 rng(rank + 12345);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    
    for (int i = 1; i <= local_Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            con(i, j) = con_0 + noise * (0.5f - dist(rng));
        }
    }
    
    apply_periodic_y(con, local_Nx);
    
    int left  = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int right = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
    
    MPI_Status status;
    const float inv_dxdy = 1.0f / (dx * dy);
    
    //--- Main evolution loop
    for (int tsteps = 1; tsteps <= nsteps; ++tsteps) {
        // Exchange halos for concentration
        exchange_x_halos(con, local_Nx, left, right, MPI_COMM_WORLD, status);
        
        // Compute chemical potential and first Laplacian
        compute_laplacian_and_chemical_potential(
            con, dfdcon, lap_con, dummy_con, 
            local_Nx, dx, dy, A, grad_coef
        );
        
        // Exchange halos for dummy_con
        exchange_x_halos(dummy_con, local_Nx, left, right, MPI_COMM_WORLD, status);
        
        // Compute second Laplacian and update
        update_field(con, dummy_con, con_new, local_Nx, dt, mobility, inv_dxdy);
        
        // Swap instead of copy
        std::swap(con, con_new);
        
        // Apply y-boundary conditions
        apply_periodic_y(con, local_Nx);
        
        if (tsteps % nprint == 0 && rank == 0) {
            std::cout << "Done steps = " << tsteps << std::endl;
        }
    }
    
    float mpi_end_time = MPI_Wtime();
    
    //--- Final Output: Rank 0 prints the top-left 5x5 sub-matrix
    if (rank == 0) {
        std::cout << "---------------------------------" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Total Run Time = " << (mpi_end_time - mpi_start_time) << " seconds." << std::endl;
        std::cout << "---------------------------------" << std::endl;
        std::cout << "Final 5x5 Matrix Output (Top-Left corner of grid):" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        
        // Print first 5 rows (only if we have them!)
        int rows_to_print = std::min(5, local_Nx);
        for (int i = 1; i <= rows_to_print; ++i) {
            for (int j = 0; j < 5; ++j) {
                std::cout << std::setw(10) << con(i, j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "---------------------------------" << std::endl;
        
        // Write subset to file - ONLY if we have enough rows!
        int output_size = std::min(100, local_Nx);
        std::ofstream outfile("ch_small.dat");
        if (outfile) {
            outfile << std::fixed << std::setprecision(6);
            for (int i = 1; i <= output_size; ++i) {
                for (int j = 0; j < std::min(100, Ny); ++j) {
                    //outfile << std::format("{:.6f}{}", con(i, j), (j < 99 ? " " : ""));
					outfile << std::fixed << std::setprecision(6) << con(i, j) << (j < 99 ? " " : "");
                }
                outfile << "\n";
            }
            outfile.close();
            std::cout << "Output written to ch_small.dat (" << output_size << "x" 
                      << std::min(100, Ny) << ")" << std::endl;
        }
    }
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}