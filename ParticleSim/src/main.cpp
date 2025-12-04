#include "Simulation.h"
#include <filesystem>
#include <iostream>
#include <string>

void ensureDirectory(const std::string &path) {
  if (!std::filesystem::exists(path)) {
    std::filesystem::create_directories(path);
  }
}

int main(int argc, char* argv[]) {
  const int N = 1000000;
  const int frames = 999;
  const float dt = 0.01f;
  const std::string outputFolder = "data";

  ensureDirectory(outputFolder);

  Simulation sim(N);

  bool runTests = false;
  bool skipCPU = false;
  bool doIO = true;
  for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "test") {
        runTests = true;
      }
      if (std::string(argv[i]) == "nocpu") {
        skipCPU = true;
      }
      if (std::string(argv[i]) == "noio") {
        doIO = false;
      }
      if (std::string(argv[i]) == "gpu") {
        sim.runGPU(frames, dt, outputFolder, doIO);
      }
  }

  if (runTests) {
    sim.testLoops(frames, dt, outputFolder, skipCPU, doIO);
  }
  else {
      std::cout << "Use GPU? y or n: ";
    
      char userGPUChar;
      std::cin >> userGPUChar;
      if (userGPUChar == 'y' || userGPUChar == 'Y') {
        sim.runGPU(frames, dt, outputFolder, doIO);
      } else {
        sim.runCPU(frames, dt, outputFolder);
      }
  }

  return 0;
}
