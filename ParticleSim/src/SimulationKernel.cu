#include "SimulationKernel.cuh"
#include "Particle.h"
#include "Config.h"
#include "CudaUtils.cuh"
#include <cuda_runtime.h>

__device__ void updateColour(Particle *p, float maxSpeed)
{
    float speed = sqrtf(p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    
    float t = speed / maxSpeed;
    t = fminf(t, 1.0f);

    float r = 0.0f, g = 0.0f, b = 0.0f;

    if (t < 0.25f) {
        r = 0.0f;
        g = 4.0f * t;
        b = 1.0f;
    } else if (t < 0.5f) {
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - 4.0f * (t - 0.25f);
    } else if (t < 0.75f) {
        r = 4.0f * (t - 0.5f);
        g = 1.0f;
        b = 0.0f;
    } else {
        r = 1.0f;
        g = 1.0f - 4.0f * (t - 0.75f);
        b = 0.0f;
    }

    p->r = r;
    p->g = g;
    p->b = b;
}

__global__ void updateKernel(Particle* particles, int numParticles, float dt, BoundsCuda bounds, ConfigCuda config)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < numParticles; i += stride)
    {
        Particle p = particles[i]; // Load into registers
    
        // Apply Gravity
        p.vy += config.gravity * dt;
    
        // Move Particle
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    
        if (p.x < bounds.xmin) { p.x = bounds.xmin; p.vx *= -0.5f; }
        else if (p.x > bounds.xmax) { p.x = bounds.xmax; p.vx *= -0.5f; }
    
        if (p.y < bounds.ymin) { p.y = bounds.ymin; p.vy *= -0.5f; }
        else if (p.y > bounds.ymax) { p.y = bounds.ymax; p.vy *= -0.5f; }
    
        if (p.z < bounds.zmin) { p.z = bounds.zmin; p.vz *= -0.5f; }
        else if (p.z > bounds.zmax) { p.z = bounds.zmax; p.vz *= -0.5f; }
    
        updateColour(&p, config.maxSpeed);
    
        particles[i] = p; // Write back once
    }

}

SimulationKernel::SimulationKernel(Particle* particles, size_t numParticles)
        : _particles(particles), _numParticles(numParticles) 
{
    cudaMalloc(&_d_particles, numParticles * sizeof(Particle));

    _bounds = BoundsCuda{
        Config::xmin, Config::xmax,
        Config::ymin, Config::ymax,
        Config::zmin, Config::zmax,
    };

    _config = ConfigCuda{
        Config::gravity,
        Config::maxSpeed,
    };

    auto launchConfig = CudaUtils::getOptimalConfig(updateKernel, numParticles, "UpdateSimulationKernel");
    this->_blockSize = launchConfig.blockSize;
    this->_gridSize = launchConfig.gridSize;
};

SimulationKernel::~SimulationKernel()
{
    if (_d_particles)
    {
        cudaFree(_d_particles);
    }
}

void SimulationKernel::copyDeviceToHost()
{
    cudaMemcpy(_particles, _d_particles, _numParticles * sizeof(Particle), cudaMemcpyDeviceToHost);
}
    
void SimulationKernel::copyHostToDevice()
{
    cudaMemcpy(_d_particles, _particles, _numParticles * sizeof(Particle), cudaMemcpyHostToDevice);
}

void SimulationKernel::simulationUpdate(float dt)
{
    copyHostToDevice();
    updateKernel<<<_gridSize, _blockSize>>>(_d_particles, _numParticles, dt, _bounds, _config);
    cudaDeviceSynchronize();
    copyDeviceToHost();
}

// Bekijken van scaling CUDA gebruik
// Scaling multi-GPU's
// single GPU - Streaming Multi-processors (1, 2, etc.)
// single GPU - Multi kernels gebruiken om de helft "Busy" te houden (cudaStreams)
// Verschillende GPU's (kloksnelheden?)