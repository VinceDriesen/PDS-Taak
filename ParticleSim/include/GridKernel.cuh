#ifndef GRID_KERNEL_CUH
#define GRID_KERNEL_CUH

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
    
    float getCellSize() const { return cellSize; }
    int getNx() const { return Nx; }
    int getNy() const { return Ny; }
    int getNz() const { return Nz; }

private:
    int Nx, Ny, Nz;
    float cellSize, invCellSize;
    float xmin, ymin, zmin;

    int totalCells;
    int numParticles;

    int blockSize;
    int gridSize;

    int* d_gridHead;      
    int* d_particleNext;  
};

#endif // GRID_KERNEL_CUH