#include "Simulation.h"
#include "Collider.h"
#include "Config.h"
#include "ExporterVTK.h"

Simulation::Simulation(int N) {
  // Initialize particles in a small grid
  initializeParticles(N);

  // Create the bounding box
  ExporterVTK::saveBoxVTK("data/bounding_box.vtk", Config::xmin, Config::xmax,
                          Config::ymin, Config::ymax, Config::zmin,
                          Config::zmax);
}

void Simulation::initializeParticles(int N) {
  particles.clear();
  float spacing = Config::diameter; // gebruik de diameter als spacing

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Particle p;
      p.x = (i % 10) * spacing;
      p.y = 1.0f + (i / 10) * spacing;
      p.z = (j % 10) * spacing;
      p.vx = p.vy = p.vz = 0.0f;
      particles.push_back(p);
    }
  }
}

void Simulation::update(float dt) {
  for (auto &p : particles) {
    p.vy += Config::gravity * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;

    if (Collider::isOutOfBounds(p.x, p.y, p.z)) {
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
