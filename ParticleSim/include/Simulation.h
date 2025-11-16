#ifndef SIMULATION_H
#define SIMULATION_H

#include "Particle.h"
#include <string>
#include <vector>

class Simulation {
public:
  Simulation(int N);
  void runCPU(int frames, float dt, const std::string &outputFolder);

private:
  std::vector<Particle> particles;
  void update(float dt);
};

#endif // SIMULATION_H
