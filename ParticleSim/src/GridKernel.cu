#include "GridKernel.cuh"
#include "Config.h"
#include "Particle.h"
#include "CudaUtils.cuh"
#include <cuda_runtime.h>

__device__ int getCellIndex(float x, float y, float z, int Nx, int Ny, int Nz) 
{
    const float invCellSize = 1.0f / Config::smoothingRadius;

    int cx = (int)((x - Config::xmin) * invCellSize);
    int cy = (int)((y - Config::ymin) * invCellSize);
    int cz = (int)((z - Config::zmin) * invCellSize);

    // Clamp to ensure we stay within valid grid bounds
    cx = max(0, min(Nx - 1, cx));
    cy = max(0, min(Ny - 1, cy));
    cz = max(0, min(Nz - 1, cz));

    return cx + cy * Nx + cz * Nx * Ny;
}

__global__ void buildGridKernel(Particle* particles, int numParticles,
                                int* gridHead, int* particleNext,
                                int Nx, int Ny, int Nz) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < numParticles; i += stride)
    {
        // Calculate which cell this particle belongs to
        // We pass Nx, Ny, Nz because grid dimensions depend on the bounds size
        int cellIndex = getCellIndex(particles[i].x, particles[i].y, particles[i].z, Nx, Ny, Nz);
    
        // Insert into linked list (Atomic Exchange)
        int nextParticleID = atomicExch(&gridHead[cellIndex], i);
        particleNext[i] = nextParticleID;
    }
}

GridKernel::GridKernel(int numParticles)
    : numParticles(numParticles)
{
    // Store local copies for getters/host logic, though kernel uses Config directly
    this->cellSize = Config::smoothingRadius;
    this->xmin = Config::xmin;
    this->ymin = Config::ymin;
    this->zmin = Config::zmin;
    this->invCellSize = 1.0f / cellSize;
    
    // Calculate dimensions on host
    this->Nx = (int)std::ceil((Config::xmax - Config::xmin) * invCellSize);
    this->Ny = (int)std::ceil((Config::ymax - Config::ymin) * invCellSize);
    this->Nz = (int)std::ceil((Config::zmax - Config::zmin) * invCellSize);
    
    this->totalCells = Nx * Ny * Nz;

    cudaMalloc(&d_gridHead, totalCells * sizeof(int));
    cudaMalloc(&d_particleNext, numParticles * sizeof(int));

    auto launchConfig = CudaUtils::getOptimalConfig(buildGridKernel, numParticles);
    this->blockSize = launchConfig.blockSize;
    this->gridSize = launchConfig.gridSize;
}

GridKernel::~GridKernel() {
    cudaFree(d_gridHead);
    cudaFree(d_particleNext);
}

void GridKernel::build(Particle* d_particles) {
    cudaMemset(d_gridHead, -1, totalCells * sizeof(int));

    buildGridKernel<<<gridSize, blockSize>>>(
        d_particles, numParticles, 
        d_gridHead, d_particleNext,
        Nx, Ny, Nz
    );
}