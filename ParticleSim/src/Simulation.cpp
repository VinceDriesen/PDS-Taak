#include "Simulation.h"
#include "Collider.h"
#include "Config.h"
#include "ExporterVTK.h"
#include <iostream>

Simulation::Simulation(int N)
    : grid(Config::diameter, Config::xmin, Config::xmax, Config::ymin,
           Config::ymax, Config::zmin, Config::zmax) {
  // Initialize particles in a small grid
  initializeParticles(N);

  // Create the bounding box
  ExporterVTK::saveBoxVTK("data/bounding_box.vtk", Config::xmin, Config::xmax,
                          Config::ymin, Config::ymax, Config::zmin,
                          Config::zmax);
}

void Simulation::initializeParticles(int N) {
  constexpr int index = 0;
  particles.clear();
  float spacing = Config::diameter; // gebruik de diameter als spacing
  switch (index) {
  case 0:
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
    break;
  case 1:
    std::vector<std::vector<float>> positions = {
        {Config::xmax / 2, Config::ymax * 0.8, Config::zmax / 2},
        {Config::xmax / 2, Config::ymax * 0.5, Config::zmax / 2},
    };
    for (const auto &pos : positions) {
      Particle p;
      p.x = pos[0];
      p.y = pos[1];
      p.z = pos[2];
      p.vx = p.vy = p.vz = 0.0f;
      p.r = rand() % 255;
      p.g = rand() % 255;
      p.b = rand() % 255;
      particles.push_back(p);
    }
    break;
  }
}

void Simulation::update(float dt) {

  // Apply gravity + movement
  for (auto &p : particles) {
    p.vy += Config::gravity * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;

    Collider::isOutOfBounds(p);
  }

  // Build grid
  grid.build(particles);

  // Collider::checkCollions(grid, particles);
  Collider::applyPressure(grid, particles, dt);
  // Handle collisions
}

void Simulation::runCPU(int frames, float dt, const std::string &outputFolder) {
  std::cout << "Running CPU simulation for " << frames << " frames.\n";
  for (int frame = 0; frame < frames; ++frame) {
    update(dt);
    ExporterVTK::saveVTKFile(particles, outputFolder, frame);
  }
}
