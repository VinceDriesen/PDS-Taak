#ifndef COLLIDER_H
#define COLLIDER_H

#include "Grid.h"
#include "Particle.h"
#include <vector>

class Collider {
public:
  static void isOutOfBounds(Particle &particle);
  static void checkCollions(const Grid &grid, std::vector<Particle> &particles);

  static void applyBoundaryCollision(Particle &p);
  static void applyPressure(const Grid &grid, std::vector<Particle> &particles,
                            float dt);

private:
  static bool isXOutOfBounds(Particle &particle);
  static bool isYOutOfBounds(Particle &particle);
  static bool isZOutOfBounds(Particle &particle);

  static void reverseXVelocity(Particle &particle);
  static void reverseYVelocity(Particle &particle);
  static void reverseZVelocity(Particle &particle);
  static void resolveParticleCollision(Particle &A, Particle &B);
};

#endif // COLLIDER_Hclass Collider {
