#include "Simulation.h"
#include "Collider.h"
#include "Config.h"
#include "ExporterVTK.h"
#include "SimulationKernel.cuh"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <chrono>

Simulation::Simulation(int N)
    : grid(Config::diameter, Config::xmin, Config::xmax, Config::ymin,
           Config::ymax, Config::zmin, Config::zmax),
      sm(nullptr), grid_kernel(nullptr), collider_kernel(nullptr) {
  // Initialize particles
  initializeParticles(N);

  // Create the bounding box
  ExporterVTK::saveBoxVTK("data/bounding_box.vtk", Config::xmin, Config::xmax,
                          Config::ymin, Config::ymax, Config::zmin,
                          Config::zmax);
}

void Simulation::initializeParticles(int N) {
  constexpr int index = 0;
  particles.clear();

  // Gebruik iets meer ruimte dan de diameter om explosies bij start te
  // voorkomen
  float spacing = Config::diameter;

  switch (index) {
  case 0: { // Bolvormige start
    // Bepaal het midden van de wereld
    float cx = (Config::xmax - Config::xmin) / 2.0f;
    float cy = (Config::ymax - Config::ymin) / 2.0f;
    float cz = (Config::zmax - Config::zmin) / 2.0f;

    // Straal van de bol
    float sphereRadius =
        spacing * std::pow((3.0f * N) / (4.0f * M_PI), 1.0f / 3.0f);

    // Bounding box rondom de bol berekenen
    int steps = (int)(sphereRadius / spacing);

    for (int i = -steps; i <= steps; i++) {
      for (int j = -steps; j <= steps; j++) {
        for (int k = -steps; k <= steps; k++) {

          // Bereken potentiële positie
          float px = cx + i * spacing;
          float py = cy + j * spacing;
          float pz = cz + k * spacing;

          // Bereken afstand tot het midden (Pythagoras in 3D)
          float distSq = (px - cx) * (px - cx) + (py - cy) * (py - cy) +
                         (pz - cz) * (pz - cz);

          // Als het deeltje binnen de straal valt, voegen we het toe
          if (distSq <= sphereRadius * sphereRadius) {
            Particle p;

            // --- SYMMETRIE BREKEN ---
            // Voeg een minuscuul beetje ruis toe ("Jitter").
            // Dit zorgt ervoor dat buren niet meer PRECIES tegenover elkaar
            // staan, waardoor de pressure force ze opzij kan duwen.
            float jitter = spacing * 0.01f; // 1% variatie
            float rx = ((rand() % 100) / 100.0f - 0.5f) * jitter;
            float ry = ((rand() % 100) / 100.0f - 0.5f) * jitter;
            float rz = ((rand() % 100) / 100.0f - 0.5f) * jitter;

            p.x = px + rx;
            p.y = py + ry;
            p.z = pz + rz;

            p.vx = p.vy = p.vz = 0.0f;
            particles.push_back(p);
          }
        }
      }
    }
    std::cout << "Initialized sphere with " << particles.size()
              << " particles.\n";
    break;
  }
  case 1: { // Twee druppels (Test case)
    std::vector<std::vector<float>> positions = {
        {Config::xmax / 2, Config::ymax * 0.8f, Config::zmax / 2},
        {Config::xmax / 2, Config::ymax * 0.5f, Config::zmax / 2},
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
}

void Simulation::update(float dt) {
  // 1. Grid bouwen (EERST doen, zodat buren kloppen voor drukberekening)
  grid.build(particles);

  // 2. Krachten berekenen (Druk + Viscositeit + Repulsion)
  // Dit past de vx, vy, vz van de particles aan.
  Collider::applyPressure(grid, particles, dt);

  // 3. Bewegen & Zwaartekracht & Grenzen
  for (auto &p : particles) {
    // Zwaartekracht toepassen
    p.vy += Config::gravity * dt;

    // Positie updaten op basis van nieuwe snelheid
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    p.z += p.vz * dt;

    // Controleren of ze de bak uit vliegen
    Collider::isOutOfBounds(p);
    auto color = getRGBFromSpeed(p.vx, p.vy, p.vz);
    p.r = color.r;
    p.g = color.g;
    p.b = color.b;
  }

  // Optioneel: Backup collision check voor als ze door elkaar vliegen
  Collider::checkCollions(grid, particles);
}

void Simulation::updateGpu(float dt)
{
    if(!sm || !grid_kernel || !collider_kernel) 
    {
        throw std::runtime_error("SimulationKernel, GridKernel, or ColliderKernel not defined");
    }

    Particle* d_particles = sm->getDevicePtr();

    grid_kernel->build(d_particles);

    collider_kernel->update(
        d_particles, 
        grid_kernel->getGridHead(),
        grid_kernel->getParticleNext(),
        grid_kernel->getCellSize(),     
        grid_kernel->getNx(), 
        grid_kernel->getNy(), 
        grid_kernel->getNz(), 
        dt
    );

    sm->copyDeviceToHost();
}

void Simulation::runCPU(int frames, float dt, const std::string &outputFolder) {
  std::cout << "Running CPU simulation for " << frames << " frames.\n";
  
  auto start = std::chrono::high_resolution_clock::now();

  for (int frame = 0; frame < frames; ++frame) {
    update(dt);
    ExporterVTK::saveVTKFile(particles, outputFolder, frame);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "CPU Execution completed in: " << diff.count() << " seconds.\n";
}

void Simulation::runGPU(int frames, float dt, const std::string &outputFolder) {
  sm = std::make_unique<SimulationKernel>(particles.data(), particles.size());

  sm->copyHostToDevice();

  grid_kernel = std::make_unique<GridKernel>(
      particles.size()
  );

  collider_kernel = std::make_unique<ColliderKernel>(
      particles.size()
  );
  std::cout << "Running CPU + GPU simulation for " << frames << " frames.\n";
  
  auto start = std::chrono::high_resolution_clock::now();

  for (int frame = 0; frame < frames; ++frame) {
    updateGpu(dt);
    ExporterVTK::saveVTKFile(particles, outputFolder, frame);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "GPU Execution completed in: " << diff.count() << " seconds.\n";
}

const RGB Simulation::getRGBFromSpeed(float vx, float vy, float vz) {
  float speed = std::sqrt(vx * vx + vy * vy + vz * vz);

  float t = speed / Config::maxSpeed;
  t = std::max(0.0f, std::min(t, 1.0f));

  RGB c;

  if (t < 0.33f) {
    float local_t = t / 0.33f;
    c.r = 0.0f;
    c.g = local_t; // 0.0 -> 1.0
    c.b = 1.0f;
  } else if (t < 0.66f) {
    float local_t = (t - 0.33f) / 0.33f;
    c.r = local_t; // 0.0 -> 1.0
    c.g = 1.0f;
    c.b = 1.0f - local_t; // 1.0 -> 0.0
  } else {
    float local_t = (t - 0.66f) / 0.34f;
    c.r = 1.0f;
    c.g = 1.0f - local_t; // 1.0 -> 0.0
    c.b = 0.0f;
  }

  return c;
}