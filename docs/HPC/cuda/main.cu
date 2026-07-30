//-------------------------------------------------------------------------------
//
//   CUDA Finite Difference Phase Field Code of Cahn-Hilliard Equation.
//
//   Author:
//               Shahid Maqbool
//
//   Modified:
//                    08 July 2026
//
//   To compile and run:
//                            nvcc -std=c++17 -o main main.cu
//                            ./main
//
//-------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cuda_runtime.h>

using namespace std::chrono;

//--- simulation cell parameters
constexpr int Nx = 128;
constexpr int Ny = 128;
constexpr double dx = 1.0;
constexpr double dy = 1.0;

//--- time integration parameters
constexpr int nsteps = 5000;
constexpr int nprint = 1000;
constexpr double dt = 0.01;

//--- material specific parameters
constexpr double con_0 = 0.4;
constexpr double mobility = 1.0;
constexpr double grad_coef = 0.5;

//--- microstructure parameters
constexpr double noise = 0.02;
constexpr double A = 1.0;

//============================================================================
// CUDA error checking macro
//============================================================================
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error in " << __FILE__ << " line " << __LINE__ \
                      << ": " << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

//============================================================================
// CUDA Kernel: Compute free energy derivative and first Laplacian
//============================================================================
__global__ void compute_first_laplacian_kernel(
    const double* con,
    double* dfdcon,
    double* lap_con,
    int nx, int ny,
    double dx, double dy,
    double A_val, double grad_coef_val)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) {
	    //=================================
        // Periodic boundary conditions
		//=================================
        int ip = (i + 1) % nx;
        int im = (i - 1 + nx) % nx;
        int jp = (j + 1) % ny;
        int jm = (j - 1 + ny) % ny;
        
        double c = con[i * ny + j];
        //========================
        // Free energy derivative
        //========================
        dfdcon[i * ny + j] = A_val * (2.0 * c * (1.0 - c) * (1.0 - c) - 
                                      2.0 * c * c * (1.0 - c));
        //===============================
        //  Laplacian
        //===============================
        lap_con[i * ny + j] = (con[ip * ny + j] + con[im * ny + j] + 
                               con[i * ny + jp] + con[i * ny + jm] - 
                               4.0 * c) / (dx * dy);
    }
}

//=============================================================================
//  CUDA Kernel: Compute second Laplacian and time integration
//=============================================================================
__global__ void compute_second_laplacian_kernel(
    const double* con,
    const double* dfdcon,
    const double* lap_con,
    double* dummy_con,
    double* lap_dummy,
    double* con_new,
    int nx, int ny,
    double dx, double dy,
    double grad_coef_val,
    double dt_val, double mobility_val)
{
    // Global thread index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i < nx && j < ny) {
	    //=================================
        // Periodic boundary conditions
		//=================================
        int ip = (i + 1) % nx;
        int im = (i - 1 + nx) % nx;
        int jp = (j + 1) % ny;
        int jm = (j - 1 + ny) % ny;
        
		//============================
        // Compute dummy concentration
		//============================
        dummy_con[i * ny + j] = dfdcon[i * ny + j] - grad_coef_val * lap_con[i * ny + j];
        
		//======================
        // Laplacian of dummy
		//======================
        lap_dummy[i * ny + j] = (dummy_con[ip * ny + j] + dummy_con[im * ny + j] + 
                                 dummy_con[i * ny + jp] + dummy_con[i * ny + jm] - 
                                 4.0 * dummy_con[i * ny + j]) / (dx * dy);
        //==================================
        // Time integration with clamping
        //==================================
        double new_val = con[i * ny + j] + dt_val * mobility_val * lap_dummy[i * ny + j];
        con_new[i * ny + j] = fmax(0.00001, fmin(0.99999, new_val));
    }
}

//================================================================================
//        Host function to output concentration to file
//================================================================================
void output_concentration_on_file(const double* con, int nx, int ny, const std::string& filename = "ch.dat") {
    std::ofstream outfile(filename);
    
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    outfile << std::fixed << std::setprecision(6);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            outfile << con[i * ny + j] << (j < ny - 1 ? " " : "");
        }
        outfile << "\n";
    }
    
    outfile.close();
    std::cout << "Results written to: " << filename << std::endl;
}

//============================================================================
//                           MAIN PROGRAM 
//============================================================================
int main() {
    //Start timing
    auto start_time = high_resolution_clock::now();
    
    //Allocate host memory (1D flattened arrays for better performance)
    size_t N = Nx * Ny;
    size_t bytes = N * sizeof(double);
    
    double* h_con = new double[N];
    double* h_con_new = new double[N];
    double* h_dfdcon = new double[N];
    double* h_lap_con = new double[N];
    double* h_dummy_con = new double[N];
    double* h_lap_dummy = new double[N];
    
	//=====================================================
    //  Initialize concentration field with noise on host
	//=====================================================
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double r = dist(rng);
            h_con[i * Ny + j] = con_0 + noise * (0.5 - r);
        }
    }
    
	//==================================
    //  Allocate device memory
	//==================================
    double* d_con;
    double* d_con_new;
    double* d_dfdcon;
    double* d_lap_con;
    double* d_dummy_con;
    double* d_lap_dummy;
    
    CUDA_CHECK(cudaMalloc(&d_con, bytes));
    CUDA_CHECK(cudaMalloc(&d_con_new, bytes));
    CUDA_CHECK(cudaMalloc(&d_dfdcon, bytes));
    CUDA_CHECK(cudaMalloc(&d_lap_con, bytes));
    CUDA_CHECK(cudaMalloc(&d_dummy_con, bytes));
    CUDA_CHECK(cudaMalloc(&d_lap_dummy, bytes));
    
	//================================
    //  Copy initial data to device
	//================================
    CUDA_CHECK(cudaMemcpy(d_con, h_con, bytes, cudaMemcpyHostToDevice));
    
	//=======================================
    // Configure kernel launch parameters
	//=======================================
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((Nx + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (Ny + threadsPerBlock.y - 1) / threadsPerBlock.y);
    //====================================
    // Start microstructure evolution
    //====================================
    for (int tsteps = 1; tsteps <= nsteps; tsteps++) {
        //============================================================
        // First kernel: compute free energy derivative and Laplacian
        //============================================================
        compute_first_laplacian_kernel<<<numBlocks, threadsPerBlock>>>(
            d_con, d_dfdcon, d_lap_con,
            Nx, Ny, dx, dy, A, grad_coef);
        CUDA_CHECK(cudaGetLastError());
        
        //=================================================================
        // Second kernel: compute second Laplacian and time integration
        //=================================================================
        compute_second_laplacian_kernel<<<numBlocks, threadsPerBlock>>>(
            d_con, d_dfdcon, d_lap_con,
            d_dummy_con, d_lap_dummy, d_con_new,
            Nx, Ny, dx, dy, grad_coef, dt, mobility);
        CUDA_CHECK(cudaGetLastError());
        
        //============================================================
        // Swap pointers: new becomes current for next iteration
        //============================================================
        std::swap(d_con, d_con_new);
        
        //======================================================
        // Print progress
        //======================================================
        if (tsteps % nprint == 0) {
            std::cout << "Done steps = " << tsteps << std::endl;
        }
    }
    //=================================
    // Copy final result back to host
    //=================================
    CUDA_CHECK(cudaMemcpy(h_con, d_con, bytes, cudaMemcpyDeviceToHost));
    
    //================================
    //  End timing
    //================================
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    
    //===============================
    //  Print results
    //===============================
    std::cout << "---------------------------------" << std::endl;
    std::cout << "  Time       = " << duration.count() / 1000.0 << " seconds." << std::endl;
    
    //==========================
    //  Write results to file
    //==========================
    output_concentration_on_file(h_con, Nx, Ny);
    
    //==================================
    //  Free device memory
    //==================================
    CUDA_CHECK(cudaFree(d_con));
    CUDA_CHECK(cudaFree(d_con_new));
    CUDA_CHECK(cudaFree(d_dfdcon));
    CUDA_CHECK(cudaFree(d_lap_con));
    CUDA_CHECK(cudaFree(d_dummy_con));
    CUDA_CHECK(cudaFree(d_lap_dummy));
    
    //==================================
    //    Free host memory
    //==================================
    delete[] h_con;
    delete[] h_con_new;
    delete[] h_dfdcon;
    delete[] h_lap_con;
    delete[] h_dummy_con;
    delete[] h_lap_dummy;
    
    return EXIT_SUCCESS;
}