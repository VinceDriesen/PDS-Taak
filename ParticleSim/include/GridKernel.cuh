#ifndef GRID_KERNEL_CUH
#define GRID_KERNEL_CUH

#include "CudaUtils.cuh"
#include "Particle.h"
#include <cuda_runtime.h>

class GridKernel {
public:
    GridKernel(int numParticles);
    ~GridKernel();
    GridKernel(const GridKernel&) = delete;
    GridKernel& operator=(const GridKernel&) = delete;

    void build(Particle* d_particles);

    int* getGridHead() const { return d_gridHead; }
    int* getParticleNext() const { return d_particleNext; }
    
    int getNx() const { return Nx; }
    int getNy() const { return Ny; }
    int getNz() const { return Nz; }

private:
    int Nx, Ny, Nz;
    float invCellSize;

    int totalCells;
    int numParticles;

    LaunchConfig launchConfig;

    int blockSize;
    int gridSize;

    int* d_gridHead;      
    int* d_particleNext;  
};

#endif // GRID_KERNEL_CUH