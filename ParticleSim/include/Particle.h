#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle {
  Particle()
      : x(0), y(0), z(0), vx(0), vy(0), vz(0), r(1.0f), g(1.0f), b(1.0f) {}
  float x, y, z;
  float vx, vy, vz;
  float r, g, b;
};

#endif // PARTICLE_H
