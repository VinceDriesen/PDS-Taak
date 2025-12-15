#ifndef COLLIDER_KERNEL_H
#define COLLIDER_KERNEL_H

#include <cuda_runtime.h>
#include "Particle.h"
#include "CudaUtils.cuh"

class ColliderKernel {
public:
    ColliderKernel(int numParticles); 
    ~ColliderKernel();
    ColliderKernel(const ColliderKernel&) = delete;
    ColliderKernel& operator=(const ColliderKernel&) = delete;

    void update(Particle* d_particles, 
                int* d_gridHead, int* d_particleNext, int Nx, int Ny, int Nz, float dt);

private:
    // Arrays to store intermetiate results, allocated on GPU
    float *d_rho, *d_pressure;
    float *d_fx, *d_fy, *d_fz;

    int numParticles;

    LaunchConfig launchConfigDensity;
    LaunchConfig launchConfigPressure;
    LaunchConfig launchConfigFoce;
    LaunchConfig launchConfigIntegrate;
};

#endif // COLLIDER_KERNEL_H