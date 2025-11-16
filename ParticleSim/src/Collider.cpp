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
  std::cout << "Particle Y velocity reversed\n";
  std::cout << particle.vy << " -> ";
  particle.vy *= -0.6f; // Demping factor
  std::cout << particle.vy << " -> ";
  std::cout << "einde" << std::endl;
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

//
// void Grid::collidePair(Particle &A, Particle &B) {
//   float dx = B.x - A.x;
//   float dy = B.y - A.y;
//   float dz = B.z - A.z;
//
//   float dist2 = dx * dx + dy * dy + dz * dz;
//   float minDist = Config::diameter;
//   float minDist2 = minDist * minDist;
//
//   if (dist2 >= minDist2)
//     return;
//
//   float dist = std::sqrt(dist2);
//   if (dist == 0.0f)
//     dist = 0.0001f;
//
//   float nx = dx / dist;
//   float ny = dy / dist;
//   float nz = dz / dist;
//
//   float overlap = (minDist - dist) * 0.5f;
//
//   A.x -= nx * overlap;
//   A.y -= ny * overlap;
//   A.z -= nz * overlap;
//
//   B.x += nx * overlap;
//   B.y += ny * overlap;
//   B.z += nz * overlap;
//
//   float dvx = B.vx - A.vx;
//   float dvy = B.vy - A.vy;
//   float dvz = B.vz - A.vz;
//
//   float relVel = dvx * nx + dvy * ny + dvz * nz;
//   if (relVel > 0)
//     return;
//
//   float bounce = 0.6f;
//   float j = -(1 + bounce) * relVel * 0.5f;
//
//   A.vx -= j * nx;
//   A.vy -= j * ny;
//   A.vz -= j * nz;
//
//   B.vx += j * nx;
//   B.vy += j * ny;
//   B.vz += j * nz;
// }
//
// void Grid::handleCollisions(std::vector<Particle> &particles) {
//
//   for (int cz = 0; cz < Nz; ++cz) {
//     for (int cy = 0; cy < Ny; ++cy) {
//       for (int cx = 0; cx < Nx; ++cx) {
//
//         int cellID = cx + cy * Nx + cz * Nx * Ny;
//
//         for (int dz = -1; dz <= 1; ++dz) {
//           for (int dy = -1; dy <= 1; ++dy) {
//             for (int dx = -1; dx <= 1; ++dx) {
//
//               int nx = cx + dx;
//               int ny = cy + dy;
//               int nz = cz + dz;
//
//               if (nx < 0 || nx >= Nx)
//                 continue;
//               if (ny < 0 || ny >= Ny)
//                 continue;
//               if (nz < 0 || nz >= Nz)
//                 continue;
//
//               int neighborID = nx + ny * Nx + nz * Nx * Ny;
//
//               for (int i : cells[cellID]) {
//                 for (int j : cells[neighborID]) {
//                   if (i < j)
//                     collidePair(particles[i], particles[j]);
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }
