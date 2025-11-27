#include "SimulationKernel.cuh"
#include <cstddef>
#include <cstdio>
#include <cuda_runtime.h>
#include "Config.h"
#include "Particle.h"

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
};

SimulationKernel::~SimulationKernel()
{
    if (_d_particles)
    {
        cudaFree(_d_particles);
    }
}

__device__ bool isXOutOfBounds(Particle* p, BoundsCuda &bounds) {
    if (p->x < bounds.xmin) {
        p->x = bounds.xmin;
        return true;
    } else if (p->x > bounds.xmax) {
        p->x = bounds.xmax;
        return true;
    }
    return false;
}

__device__ bool isYOutOfBounds(Particle* p, BoundsCuda &bounds) {
    if (p->y < bounds.ymin) {
        p->y = bounds.ymin;
        return true;
    } else if (p->y > bounds.ymax) {
        p->y = bounds.ymax;
        return true;
    }
    return false;
}

__device__ bool isZOutOfBounds(Particle* p, BoundsCuda &bounds) {
    if (p->z < bounds.zmin) {
        p->z = bounds.zmin;
        return true;
    } else if (p->z > bounds.zmax) {
        p->z = bounds.zmax;
        return true;
    }
    return false;
}

__device__ void checkBounds(Particle* p, BoundsCuda &bounds)
{
    if (isXOutOfBounds(p, bounds)) p->vx *= -0.5f;
    if (isYOutOfBounds(p, bounds)) p->vy *= -0.5f;
    if (isZOutOfBounds(p, bounds)) p->vz *= -0.5f;
}

__device__ void updateColour(Particle *p, float maxSpeed)
{
    float speed = sqrtf(p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    
    float t = speed / maxSpeed;
    if (t > 1.0f) t = 1.0f;

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
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    // if (idx == 0)
    // {
    //     printf("Sanity Check, GPU Werkt!");
    // }

    for (int i = idx; i < numParticles; i += stride) {
        Particle *p = &particles[i];

        // Apply Gravity
        p->vy += config.gravity * dt;

        // Move Particle
        p->x += p->vx * dt;
        p->y += p->vy * dt;
        p->z += p->vz * dt;

        checkBounds(p, bounds);

        updateColour(p, config.maxSpeed);
    }
}

void SimulationKernel::simulationUpdate(float dt)
{
    size_t sizeBytes = _numParticles * sizeof(Particle);

    int minGridSize;
    int blockSize;
    
    // Deze functie berekent de ideale configuratie om de GPU 100% vol te plannen
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, updateKernel, 0, _numParticles);
    int gridSize = (_numParticles + blockSize - 1) / blockSize;

    cudaMemcpy(_d_particles, _particles, sizeBytes, cudaMemcpyHostToDevice);

    updateKernel<<<gridSize, blockSize>>>(_d_particles, _numParticles, dt, _bounds, _config);

    cudaDeviceSynchronize();
    cudaMemcpy(_particles, _d_particles, sizeBytes, cudaMemcpyDeviceToHost);
}

// Bekijken van scaling CUDA gebruik
// Scaling multi-GPU's
// single GPU - Streaming Multi-processors (1, 2, etc.)
// single GPU - Multi kernels gebruiken om de helft "Busy" te houden (cudaStreams)
// Verschillende GPU's (kloksnelheden?)