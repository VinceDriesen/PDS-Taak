#include "Simulation.h"
#include <filesystem>
#include <string>

void ensureDirectory(const std::string &path) {
  if (!std::filesystem::exists(path)) {
    std::filesystem::create_directories(path);
  }
}

int main() {
  const int N = 100;
  const int frames = 300;
  const float dt = 0.01f;
  const std::string outputFolder = "data";

  ensureDirectory(outputFolder);

  Simulation sim(N);

  sim.runCPU(frames, dt, outputFolder);

  return 0;
}
