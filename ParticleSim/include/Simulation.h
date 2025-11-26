#ifndef SIMULATION_H
#define SIMULATION_H

#include "Grid.h"
#include "Particle.h"
#include "SimulationKernel.cuh"
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

  void runGPU(int frames, float dt, const std::string &outputFolder);

private:
  Grid grid;
  void initializeParticles(int N);
  std::vector<Particle> particles;
  void update(float dt, bool useGpu);

  std::unique_ptr<SimulationKernel> sm;

  const RGB getRGBFromSpeed(float vx, float vy, float vz);
};

#endif // SIMULATION_H
