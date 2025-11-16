#include "Simulation.h"
#include "ExporterVTK.h"

Simulation::Simulation(int N) {
  // Initialize particles in a small grid
  for (int i = 0; i < N; i++) {
    Particle p;
    p.x = (i % 10) * 0.1f;
    p.y = 1.0f + (i / 10) * 0.1f;
    p.z = 0.0f;
    p.vx = p.vy = p.vz = 0.0f;
    particles.push_back(p);
  }
}

void Simulation::update(float dt) {
  const float gravity = -9.81f;
  for (auto &p : particles) {
    p.vy += gravity * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;

    if (p.y < 0.0f) {
      p.y = 0.0f;
      p.vy *= -0.6f;
    }
  }
}

void Simulation::runCPU(int frames, float dt, const std::string &outputFolder) {
  for (int frame = 0; frame < frames; ++frame) {
    update(dt);
    ExporterVTK::saveVTKFile(particles, outputFolder, frame);
  }
}
