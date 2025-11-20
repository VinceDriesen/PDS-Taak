#include "Simulation.h"
#include "Collider.h"
#include "Config.h"
#include "ExporterVTK.h"
#include <cmath>
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
    auto color = getRGBFromSpeed(p.vx, p.vy, p.vz);
    p.r = color.r;
    p.g = color.g;
    p.b = color.b;
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

const RGB Simulation::getRGBFromSpeed(float vx, float vy, float vz) {
  float speed = std::sqrt(vx * vx + vy * vy + vz * vz);
  
  float t = speed / Config::maxSpeed; 
  t = std::max(0.0f, std::min(t, 1.0f));

  RGB c;

  if (t < 0.33f) {
      float local_t = t / 0.33f;
      c.r = 0.0f;
      c.g = local_t;       // 0.0 -> 1.0
      c.b = 1.0f;
  } 
  else if (t < 0.66f) {
      float local_t = (t - 0.33f) / 0.33f;
      c.r = local_t;       // 0.0 -> 1.0
      c.g = 1.0f;
      c.b = 1.0f - local_t; // 1.0 -> 0.0
  } 
  else {
      float local_t = (t - 0.66f) / 0.34f;
      c.r = 1.0f;
      c.g = 1.0f - local_t; // 1.0 -> 0.0
      c.b = 0.0f;
  }

  return c;
}