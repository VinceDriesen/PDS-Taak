#ifndef SIMULATION_H
#define SIMULATION_H

#include "Grid.h"
#include "Particle.h"
#include "SimulationKernel.cuh"
#include "ColliderKernel.cuh"
#include "GridKernel.cuh"
#include <memory>
#include <string>
#include <vector>

struct RGB {
  float r, g, b;
};

class Simulation {
public:
  Simulation(int N);
  void runCPU(int frames, float dt, const std::string &outputFolder);

  void runGPU(int frames, float dt, const std::string &outputFolder, bool doIO);

  void testLoops(int frames, float dt, const std::string &outputFolder, bool skipCPU, bool doIO);

private:
  Grid grid;
  void initializeParticles(int N);
  std::vector<Particle> particles;
  void update(float dt);

  std::unique_ptr<SimulationKernel> sm;
  std::unique_ptr<GridKernel> grid_kernel;
  std::unique_ptr<ColliderKernel> collider_kernel;
  const RGB getRGBFromSpeed(float vx, float vy, float vz);

  void updateGpu(float dt);
};

#endif // SIMULATION_H
