#ifndef SIMULATION_KERNEL_CUH
#define SIMULATION_KERNEL_CUH

#include "Particle.h"
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

    Particle* getDevicePtr() const { return _d_particles; }

    void copyDeviceToHost();
    void copyHostToDevice();

    void simulationUpdate(float dt);

private:
    Particle* _particles;
    size_t _numParticles;
    Particle *_d_particles;
    BoundsCuda _bounds;
    ConfigCuda _config;
    int _gridSize;
    int _blockSize;
};


#endif // SIMULATION_KERNEL_CUH