#ifndef SIMULATION_H
#define SIMULATION_H

#include "Grid.h"
#include "Particle.h"
#include <string>
#include <vector>

struct RGB {
  float r, g, b;
};

class Simulation {
public:
  Simulation(int N);
  void runCPU(int frames, float dt, const std::string &outputFolder);

private:
  Grid grid;
  void initializeParticles(int N);
  std::vector<Particle> particles;
  void update(float dt);

  const RGB getRGBFromSpeed(float vx, float vy, float vz);
};

#endif // SIMULATION_H
