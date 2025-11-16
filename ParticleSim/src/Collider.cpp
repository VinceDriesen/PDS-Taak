#include "Collider.h"
#include "Config.h"
#include <iostream>

void Collider::isOutOfBounds(Particle &particle) {
  if (isXOutOfBounds(particle)) {
    reverseXVelocity(particle);
  }
  if (isYOutOfBounds(particle)) {
    reverseYVelocity(particle);
  }
  if (isZOutOfBounds(particle)) {
    reverseZVelocity(particle);
  }
}

void Collider::checkCollions(const Grid &grid,
                             std::vector<Particle> &particles) {}

bool Collider::isXOutOfBounds(Particle &particle) {
  if (particle.x < Config::xmin) {
    particle.x = Config::xmin;
    return true;
  } else if (particle.x > Config::xmax) {
    particle.x = Config::xmax;
    return true;
  }
  return false;
}

bool Collider::isYOutOfBounds(Particle &particle) {
  if (particle.y < Config::ymin) {
    particle.y = Config::ymin;
    return true;
  } else if (particle.y > Config::ymax) {
    particle.y = Config::ymax;
    return true;
  }
  return false;
}

bool Collider::isZOutOfBounds(Particle &particle) {
  if (particle.z < Config::zmin) {
    particle.z = Config::zmin;
    return true;
  } else if (particle.z > Config::zmax) {
    particle.z = Config::zmax;
    return true;
  }
  return false;
}

void Collider::reverseXVelocity(Particle &particle) {
  particle.vx *= -0.6f; // Demping factor
}

void Collider::reverseYVelocity(Particle &particle) {
  std::cout << "Particle Y velocity reversed\n";
  std::cout << particle.vy << " -> ";
  particle.vy *= -0.6f; // Demping factor
  std::cout << particle.vy << " -> ";
  std::cout << "einde" << std::endl;
}

void Collider::reverseZVelocity(Particle &particle) {
  particle.vz *= -0.6f; // Demping factor
}
