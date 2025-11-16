#include "Simulation.h"
#include "Config.h"
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

  // Create the bounding box
  ExporterVTK::saveBoxVTK("data/bounding_box.vtk", Config::xmin, Config::xmax,
                          Config::ymin, Config::ymax, Config::zmin,
                          Config::zmax);
}

void Simulation::update(float dt) {
  for (auto &p : particles) {
    p.vy += Config::gravity * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;

    if (p.y < Config::ymin) {
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
