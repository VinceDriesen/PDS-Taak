#ifndef COLLIDER_H
#define COLLIDER_H

#include "Grid.h"
#include "Particle.h"
#include <vector>

class Collider {
public:
  static void isOutOfBounds(Particle &particle);
  static void checkCollions(const Grid &grid, std::vector<Particle> &particles);

private:
  static bool isXOutOfBounds(Particle &particle);
  static bool isYOutOfBounds(Particle &particle);
  static bool isZOutOfBounds(Particle &particle);

  static void reverseXVelocity(Particle &particle);
  static void reverseYVelocity(Particle &particle);
  static void reverseZVelocity(Particle &particle);
};

#endif // COLLIDER_H
