#ifndef COLLIDER_KERNEL_H
#define COLLIDER_KERNEL_H

#include <cuda_runtime.h>
#include "Particle.h"

class ColliderKernel {
public:
    ColliderKernel(int numParticles); 
    ~ColliderKernel();
    ColliderKernel(const ColliderKernel&) = delete;
    ColliderKernel& operator=(const ColliderKernel&) = delete;

    void update(Particle* d_particles, 
                int* d_gridHead, int* d_particleNext,
                float cellSize, int Nx, int Ny, int Nz, float dt);

private:
    float *d_rho, *d_pressure;
    float *d_fx, *d_fy, *d_fz;
    int m_numParticles;

    int blockSize_density, gridSize_density;
    int blockSize_pressure, gridSize_pressure;
    int blockSize_force, gridSize_force;
    int blockSize_integrate, gridSize_integrate;
};

#endif // COLLIDER_KERNEL_H