#include "Collider.h"
#include "Config.h"
#include <algorithm>
#include <cmath>
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
  particle.vy *= -0.6f; // Demping factor
}

void Collider::reverseZVelocity(Particle &particle) {
  particle.vz *= -0.6f; // Demping factor
}

void Collider::resolveParticleCollision(Particle &A, Particle &B) {
  float dx = B.x - A.x;
  float dy = B.y - A.y;
  float dz = B.z - A.z;

  float dist2 = dx * dx + dy * dy + dz * dz;
  float minDist = Config::diameter;
  float minDist2 = minDist * minDist;

  if (dist2 >= minDist2)
    return;

  float dist = std::sqrt(dist2);
  if (dist == 0.0f)
    dist = 0.0001f; // voorkom div/0

  // Normale vector
  float nx = dx / dist;
  float ny = dy / dist;
  float nz = dz / dist;

  // overlap verdelen
  float overlap = (minDist - dist) * 0.5f;

  A.x -= nx * overlap;
  A.y -= ny * overlap;
  A.z -= nz * overlap;

  B.x += nx * overlap;
  B.y += ny * overlap;
  B.z += nz * overlap;

  // velocity response
  float dvx = B.vx - A.vx;
  float dvy = B.vy - A.vy;
  float dvz = B.vz - A.vz;

  float relVel = dvx * nx + dvy * ny + dvz * nz;
  if (relVel > 0)
    return; // al uit elkaar, geen reactie

  float bounce = 0.6f;
  float j = -(1 + bounce) * relVel * 0.5f;

  A.vx -= j * nx;
  A.vy -= j * ny;
  A.vz -= j * nz;
  B.vx += j * nx;
  B.vy += j * ny;
  B.vz += j * nz;
}

void Collider::checkCollions(const Grid &grid,
                             std::vector<Particle> &particles) {

  const auto &cells = grid.getCells();
  int Nx = grid.getNx();
  int Ny = grid.getNy();
  int Nz = grid.getNz();

  for (int cz = 0; cz < Nz; ++cz) {
    for (int cy = 0; cy < Ny; ++cy) {
      for (int cx = 0; cx < Nx; ++cx) {

        int cellID = cx + cy * Nx + cz * Nx * Ny;

        // Loop over 27 neighbors
        for (int dz = -1; dz <= 1; ++dz) {
          for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {

              int nx = cx + dx;
              int ny = cy + dy;
              int nz = cz + dz;

              if (nx < 0 || nx >= Nx)
                continue;
              if (ny < 0 || ny >= Ny)
                continue;
              if (nz < 0 || nz >= Nz)
                continue;

              int neighborID = nx + ny * Nx + nz * Nx * Ny;

              for (int i : cells[cellID]) {
                for (int j : cells[neighborID]) {
                  if (i < j) { // voorkom dubbele checks
                    resolveParticleCollision(particles[i], particles[j]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

static float W_poly6(float r, float h) {
  if (r >= 0 && r <= h) {
    float hr2 = h * h - r * r;
    return 315.0f / (64.0f * M_PI * pow(h, 9)) * hr2 * hr2 * hr2;
  }
  return 0.0f;
}

static void gradW_spiky(float rx, float ry, float rz, float r, float h,
                        float &gx, float &gy, float &gz) {
  if (r > 0 && r <= h) {
    float coeff = -45.0f / (M_PI * pow(h, 6)) * (h - r) * (h - r) / r;
    gx = coeff * rx;
    gy = coeff * ry;
    gz = coeff * rz;
  } else {
    gx = gy = gz = 0.0f;
  }
}

void Collider::applyPressure(const Grid &grid, std::vector<Particle> &particles,
                             float dt) {
  int N = particles.size();
  float h = Config::smoothingRadius; // SPH smoothing radius
  float mass = Config::particleMass;
  float k = Config::stiffness; // pressure constant
  float rho0 = Config::restDensity;

  std::vector<float> rho(N, 0.0f);
  std::vector<float> pressure(N, 0.0f);
  std::vector<float> fx(N, 0), fy(N, 0), fz(N, 0);

  const auto &cells = grid.getCells();
  int Nx = grid.getNx(), Ny = grid.getNy(), Nz = grid.getNz();

  // 1️⃣ Compute density
  for (int cz = 0; cz < Nz; ++cz)
    for (int cy = 0; cy < Ny; ++cy)
      for (int cx = 0; cx < Nx; ++cx) {
        int cellID = cx + cy * Nx + cz * Nx * Ny;
        for (int dz = -1; dz <= 1; ++dz)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              int nx = cx + dx, ny = cy + dy, nz = cz + dz;
              if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 ||
                  nz >= Nz)
                continue;
              int neighborID = nx + ny * Nx + nz * Nx * Ny;

              for (int i : cells[cellID]) {
                for (int j : cells[neighborID]) {
                  if (i == j)
                    continue;
                  float rx = particles[i].x - particles[j].x;
                  float ry = particles[i].y - particles[j].y;
                  float rz = particles[i].z - particles[j].z;
                  float r = std::sqrt(rx * rx + ry * ry + rz * rz);
                  if (r < h) {
                    rho[i] += mass * W_poly6(r, h);
                  }
                }
              }
            }
      }

  // 2️⃣ Compute pressure
  for (int i = 0; i < N; ++i)
    pressure[i] = k * (rho[i] - rho0);

  // 3️⃣ Compute pressure forces
  for (int cz = 0; cz < Nz; ++cz)
    for (int cy = 0; cy < Ny; ++cy)
      for (int cx = 0; cx < Nx; ++cx) {
        int cellID = cx + cy * Nx + cz * Nx * Ny;
        for (int dz = -1; dz <= 1; ++dz)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              int nx = cx + dx, ny = cy + dy, nz = cz + dz;
              if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 ||
                  nz >= Nz)
                continue;
              int neighborID = nx + ny * Nx + nz * Nx * Ny;

              for (int i : cells[cellID]) {
                for (int j : cells[neighborID]) {
                  if (i >= j)
                    continue; // voorkom dubbele force
                  float rx = particles[i].x - particles[j].x;
                  float ry = particles[i].y - particles[j].y;
                  float rz = particles[i].z - particles[j].z;
                  float r = std::sqrt(rx * rx + ry * ry + rz * rz);
                  if (r < h && r > 0.0f) {
                    float gx, gy, gz;
                    gradW_spiky(rx, ry, rz, r, h, gx, gy, gz);
                    float coeff =
                        -mass * (pressure[i] + pressure[j]) / (2 * rho[j]);
                    fx[i] += coeff * gx;
                    fy[i] += coeff * gy;
                    fz[i] += coeff * gz;
                    fx[j] -= coeff * gx;
                    fy[j] -= coeff * gy;
                    fz[j] -= coeff * gz; // Newton
                  }
                }
              }
            }
      }
  float pressureScale = 0.0001f;
  // 4️⃣ Update velocities
  for (int i = 0; i < N; ++i) {
    fx[i] *= pressureScale;
    fy[i] *= pressureScale;
    fz[i] *= pressureScale;
    particles[i].vx += dt * fx[i] / mass;
    particles[i].vy += dt * fy[i] / mass;
    particles[i].vz += dt * fz[i] / mass;
  }
}
