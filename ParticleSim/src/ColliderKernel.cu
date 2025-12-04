#include "ColliderKernel.cuh"
#include "Config.h"
#include "Particle.h"
#include "CudaUtils.cuh"
#include <cmath>

__device__ inline void checkBoundary(float& pos, float& vel, float minVal, float maxVal, float damping) {
    if (pos < minVal) {
        pos = minVal;
        vel *= damping;
    } else if (pos > maxVal) {
        pos = maxVal;
        vel *= damping;
    }
}

__global__ void computeDensityKernel(
    Particle* particles, int numParticles,
    int* gridHead, int* particleNext,
    float* rho,
    int Nx, int Ny, int Nz,
    float poly6Coeff)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < numParticles; i += stride)
    {
        // Cache position
        float3 pos = make_float3(particles[i].x, particles[i].y, particles[i].z);
        
        // Constants from Config
        const float h2 = Config::smoothingRadius * Config::smoothingRadius;
        const float mass = Config::particleMass;
        const float invCellSize = 1.0f / Config::smoothingRadius;

        float density = 0.0f;

        // Grid coordinates
        int gx = (int)((pos.x - Config::xmin) * invCellSize);
        int gy = (int)((pos.y - Config::ymin) * invCellSize);
        int gz = (int)((pos.z - Config::zmin) * invCellSize);

        // Neighbor search
        for (int k = -1; k <= 1; k++) {
            for (int j = -1; j <= 1; j++) {
                for (int l = -1; l <= 1; l++) {
                    int nx = gx + l; 
                    int ny = gy + j; 
                    int nz = gz + k;

                    if (nx >= 0 && nx < Nx && ny >= 0 && ny < Ny && nz >= 0 && nz < Nz) {
                        int cellIndex = nx + ny * Nx + nz * Nx * Ny;
                        int neighbor = gridHead[cellIndex];

                        while (neighbor != -1) {
                            float3 nPos = make_float3(particles[neighbor].x, particles[neighbor].y, particles[neighbor].z);
                            
                            float dx = pos.x - nPos.x;
                            float dy = pos.y - nPos.y;
                            float dz = pos.z - nPos.z;
                            float r2 = dx*dx + dy*dy + dz*dz;

                            if (r2 < h2) {
                                float term = h2 - r2;
                                density += mass * poly6Coeff * term * term * term;
                            }
                            neighbor = particleNext[neighbor];
                        }
                    }
                }
            }
        }
        
        // Clamp minimal density to avoid division by zero later
        rho[i] = (density < 0.0001f) ? 0.0001f : density; 
    }
}

__global__ void computePressureKernel(float* rho, float* pressure, int numParticles) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < numParticles; i += stride)
    {
        pressure[i] = fmaxf(0.0f, Config::stiffness * (rho[i] - Config::restDensity));
    }
}

__global__ void computeForcesKernel(
    Particle* particles, int numParticles, 
    int* gridHead, int* particleNext,
    float* rho, float* pressure,
    float* fx, float* fy, float* fz,
    int Nx, int Ny, int Nz,
    float spikyCoeff, float viscCoeff)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < numParticles; i += stride)
    {

        // 1. Cache current particle data
        float3 pos = make_float3(particles[i].x, particles[i].y, particles[i].z);
        float3 vel = make_float3(particles[i].vx, particles[i].vy, particles[i].vz);
        float pi = pressure[i];
        
        // 2. Prepare accumulators
        float force_x = 0.0f;
        float force_y = 0.0f;
        float force_z = 0.0f;
    
        // 3. Config Constants
        const float h = Config::smoothingRadius;
        const float h2 = h * h;
        const float mass = Config::particleMass;
        const float viscosity = Config::viscosity;
        const float invCellSize = 1.0f / h;
    
        int gx = (int)((pos.x - Config::xmin) * invCellSize);
        int gy = (int)((pos.y - Config::ymin) * invCellSize);
        int gz = (int)((pos.z - Config::zmin) * invCellSize);
    
        // 4. Neighbor Search
        for (int k = -1; k <= 1; k++) {
            for (int j = -1; j <= 1; j++) {
                for (int l = -1; l <= 1; l++) {
                    int nx = gx + l; 
                    int ny = gy + j; 
                    int nz = gz + k;
    
                    if (nx >= 0 && nx < Nx && ny >= 0 && ny < Ny && nz >= 0 && nz < Nz) {
                        int cellIndex = nx + ny * Nx + nz * Nx * Ny;
                        int neighbor = gridHead[cellIndex];
    
                        while (neighbor != -1) {
                            if (neighbor != i) {
                                float3 nPos = make_float3(particles[neighbor].x, particles[neighbor].y, particles[neighbor].z);
                                
                                float dx = pos.x - nPos.x;
                                float dy = pos.y - nPos.y;
                                float dz = pos.z - nPos.z;
                                float r2 = dx*dx + dy*dy + dz*dz;
    
                                if (r2 < h2 && r2 > 1e-12f) {
                                    float r = sqrtf(r2);
                                    float h_minus_r = h - r;
    
                                    float pj = pressure[neighbor];
                                    float rho_j = rho[neighbor];
                                    float invRhoJ = 1.0f / rho_j;
                                    
                                    // --- Pressure Force ---
                                    float pressTerm = -0.5f * mass * (pi + pj) * invRhoJ;
                                    float gradW = spikyCoeff * h_minus_r * h_minus_r;
                                    float fPress = pressTerm * gradW;
                                    
                                    // Clamp Pressure
                                    if (fPress > Config::maxPressureForce) fPress = Config::maxPressureForce;
                                    if (fPress < -Config::maxPressureForce) fPress = -Config::maxPressureForce;
    
                                    float invR = 1.0f / r;
                                    force_x += fPress * dx * invR;
                                    force_y += fPress * dy * invR;
                                    force_z += fPress * dz * invR;
    
                                    // --- Viscosity Force ---
                                    float3 nVel = make_float3(particles[neighbor].vx, particles[neighbor].vy, particles[neighbor].vz);
                                    float laplacianW = viscCoeff * h_minus_r;
                                    float viscTerm = mass * viscosity * invRhoJ * laplacianW;
    
                                    force_x += viscTerm * (nVel.x - vel.x);
                                    force_y += viscTerm * (nVel.y - vel.y);
                                    force_z += viscTerm * (nVel.z - vel.z);
                                    
                                    // --- Repulsion Force (Collision) ---
                                    if (r < Config::diameter) {
                                        float repForce = Config::repulsionStiffness * (Config::diameter - r);
                                        force_x += repForce * dx * invR;
                                        force_y += repForce * dy * invR;
                                        force_z += repForce * dz * invR;
                                    }
                                }
                            }
                            neighbor = particleNext[neighbor];
                        }
                    }
                }
            }
        }
    
        fx[i] = force_x;
        fy[i] = force_y;
        fz[i] = force_z;
    }
}

__global__ void integrateKernel(Particle* particles, int numParticles,
                                float* fx, float* fy, float* fz,
                                float dt) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < numParticles; i += stride)
    {

        // 1. Read Globals to Registers
        float p_x = particles[i].x;
        float p_y = particles[i].y;
        float p_z = particles[i].z;
        float p_vx = particles[i].vx;
        float p_vy = particles[i].vy;
        float p_vz = particles[i].vz;
    
        float invMass = 1.0f / Config::particleMass;
        const float damping = -0.5f;
    
        // 2. Update Velocity
        p_vx += (fx[i] * invMass) * dt;
        p_vy += (fy[i] * invMass + Config::gravity) * dt;
        p_vz += (fz[i] * invMass) * dt;
    
        // 3. Update Position
        p_x += p_vx * dt;
        p_y += p_vy * dt;
        p_z += p_vz * dt;
    
        // 4. Boundary Checks
        checkBoundary(p_x, p_vx, Config::xmin, Config::xmax, damping);
        checkBoundary(p_y, p_vy, Config::ymin, Config::ymax, damping);
        checkBoundary(p_z, p_vz, Config::zmin, Config::zmax, damping);
    
        // 5. Write Back
        particles[i].x = p_x;
        particles[i].y = p_y;
        particles[i].z = p_z;
        particles[i].vx = p_vx;
        particles[i].vy = p_vy;
        particles[i].vz = p_vz;
    }
}

ColliderKernel::ColliderKernel(int numParticles)
    : m_numParticles(numParticles)
{
    cudaMalloc(&d_rho, m_numParticles * sizeof(float));
    cudaMalloc(&d_pressure, m_numParticles * sizeof(float));
    cudaMalloc(&d_fx, m_numParticles * sizeof(float));
    cudaMalloc(&d_fy, m_numParticles * sizeof(float));
    cudaMalloc(&d_fz, m_numParticles * sizeof(float));

    auto launchConfigDensity = CudaUtils::getOptimalConfig(computeDensityKernel, numParticles);
    this->blockSize_density = launchConfigDensity.blockSize;
    this->gridSize_density = launchConfigDensity.gridSize;

    auto launchConfigPressure = CudaUtils::getOptimalConfig(computePressureKernel, numParticles);
    this->blockSize_pressure = launchConfigPressure.blockSize;
    this->gridSize_pressure = launchConfigPressure.gridSize;

    auto launchConfigFoce = CudaUtils::getOptimalConfig(computeForcesKernel, numParticles);
    this->blockSize_force = launchConfigFoce.blockSize;
    this->gridSize_force = launchConfigFoce.gridSize;

    auto launchConfigIntegrate = CudaUtils::getOptimalConfig(integrateKernel, numParticles);
    this->blockSize_integrate = launchConfigIntegrate.blockSize;
    this->gridSize_integrate = launchConfigIntegrate.gridSize;
}

ColliderKernel::~ColliderKernel() {
    cudaFree(d_rho); 
    cudaFree(d_pressure);
    cudaFree(d_fx); 
    cudaFree(d_fy); 
    cudaFree(d_fz);
}

void ColliderKernel::update(Particle* d_particles, 
                            int* d_gridHead, int* d_particleNext,
                            float cellSize, int Nx, int Ny, int Nz, float dt) 
{
    // Compute SPH Kernel Coefficients
    float h = Config::smoothingRadius;
    float pi = (float)M_PI;
    
    float poly6Coeff = 315.0f / (64.0f * pi * powf(h, 9));
    float spikyCoeff = -45.0f / (pi * powf(h, 6));
    float viscCoeff = 45.0f / (pi * powf(h, 6));

    // 1. Density
    computeDensityKernel<<<gridSize_density, blockSize_density>>>(
        d_particles, m_numParticles, d_gridHead, d_particleNext,
        d_rho, 
        Nx, Ny, Nz, 
        poly6Coeff // Pass coeff as arg, others are in Config
    );

    // 2. Pressure
    computePressureKernel<<<gridSize_pressure, blockSize_pressure>>>(
        d_rho, d_pressure, m_numParticles
    );
    
    // 3. Forces
    computeForcesKernel<<<gridSize_force, blockSize_force>>>(
        d_particles, m_numParticles, d_gridHead, d_particleNext,
        d_rho, d_pressure, d_fx, d_fy, d_fz,
        Nx, Ny, Nz,
        spikyCoeff, viscCoeff
    );

    // 4. Integration
    integrateKernel<<<gridSize_integrate, blockSize_integrate>>>(
        d_particles, m_numParticles, d_fx, d_fy, d_fz,
        dt
    );
    
    cudaDeviceSynchronize();
}