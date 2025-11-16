#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle {
  Particle() : x(0), y(0), z(0), vx(0), vy(0), vz(0) {}
  float x, y, z;
  float vx, vy, vz;
};

#endif // PARTICLE_H
