#include "Simulation.h"
#include <filesystem>
#include <string>

void ensureDirectory(const std::string &path) {
  if (!std::filesystem::exists(path)) {
    std::filesystem::create_directories(path);
  }
}

int main() {
  const int N = 1000;
  const int frames = 999;
  const float dt = 0.01f;
  const std::string outputFolder = "data";

  ensureDirectory(outputFolder);

  Simulation sim(N);

  sim.runCPU(frames, dt, outputFolder);

  return 0;
}
