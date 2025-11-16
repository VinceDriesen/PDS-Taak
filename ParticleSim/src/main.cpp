#include "Simulation.h"
#include <string>

int main() {
  const int N = 100;
  const int frames = 300;
  const float dt = 0.01f;
  const std::string outputFolder = "data";

  Simulation sim(N);
  sim.runCPU(frames, dt, outputFolder);

  return 0;
}
