#ifndef SIMULATION_KERNEL_CUH
#define SIMULATION_KERNEL_CUH

#include "Particle.h"
#include "CudaUtils.cuh"
#include <cuda_runtime.h>


struct BoundsCuda {
    float xmin, xmax;
    float ymin, ymax;
    float zmin, zmax;
};

struct ConfigCuda {
    float gravity;
    float maxSpeed;
};

class SimulationKernel
{
public:
    SimulationKernel(Particle* particles, size_t numParticles);
    ~SimulationKernel();
    SimulationKernel(const SimulationKernel&) = delete;
    SimulationKernel& operator=(const SimulationKernel&) = delete;

    Particle* getDevicePtr() const { return d_particles; }

    void copyDeviceToHost();
    void copyHostToDevice();

    void simulationUpdate(float dt);

private:
    Particle* particles;
    size_t numParticles;

    Particle *d_particles;
    
    BoundsCuda bounds;
    ConfigCuda config;

    LaunchConfig launchConfig;
};


#endif // SIMULATION_KERNEL_CUH