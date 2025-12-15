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
        // Cache Particle
        Particle p = particles[i];
    
        // Apply Gravity
        p.vy += config.gravity * dt;
    
        // Move Particle based on stored velocity
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    
        // Check bounds, and reverse direction if it does collide
        if (p.x < bounds.xmin) { p.x = bounds.xmin; p.vx *= -0.5f; }
        else if (p.x > bounds.xmax) { p.x = bounds.xmax; p.vx *= -0.5f; }
    
        if (p.y < bounds.ymin) { p.y = bounds.ymin; p.vy *= -0.5f; }
        else if (p.y > bounds.ymax) { p.y = bounds.ymax; p.vy *= -0.5f; }
    
        if (p.z < bounds.zmin) { p.z = bounds.zmin; p.vz *= -0.5f; }
        else if (p.z > bounds.zmax) { p.z = bounds.zmax; p.vz *= -0.5f; }
    
        updateColour(&p, config.maxSpeed);
        
        // Writeback Particle
        particles[i] = p;
    }

}

SimulationKernel::SimulationKernel(Particle* particles, size_t numParticles)
        : particles(particles), numParticles(numParticles) 
{
    // Allocate GPU Memory
    cudaMalloc(&_d_particles, numParticles * sizeof(Particle));

    // Copy the Bounds and Config Settings
    bounds = BoundsCuda{
        Config::xmin, Config::xmax,
        Config::ymin, Config::ymax,
        Config::zmin, Config::zmax,
    };

    config = ConfigCuda{
        Config::gravity,
        Config::maxSpeed,
    };

    launchConfig = CudaUtils::getOptimalConfig(updateKernel, numParticles, "UpdateSimulationKernel");
};

SimulationKernel::~SimulationKernel()
{
    if (d_particles)
    {
        cudaFree(d_particles);
    }
}

void SimulationKernel::copyDeviceToHost()
{
    cudaMemcpy(particles, d_particles, numParticles * sizeof(Particle), cudaMemcpyDeviceToHost);
}
    
void SimulationKernel::copyHostToDevice()
{
    cudaMemcpy(d_particles, particles, numParticles * sizeof(Particle), cudaMemcpyHostToDevice);
}

void SimulationKernel::simulationUpdate(float dt)
{
    copyHostToDevice();
    updateKernel<<<launchConfig.gridSize, launchConfig.blockSize>>>(d_particles, numParticles, dt, bounds, config);
    cudaDeviceSynchronize();
    copyDeviceToHost();
}

// Bekijken van scaling CUDA gebruik
// Scaling multi-GPU's
// single GPU - Streaming Multi-processors (1, 2, etc.)
// single GPU - Multi kernels gebruiken om de helft "Busy" te houden (cudaStreams)
// Verschillende GPU's (kloksnelheden?)