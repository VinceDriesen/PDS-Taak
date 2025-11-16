#include "Grid.h"
#include "Config.h"
#include <cmath>

Grid::Grid(float cellSize, float xmin, float xmax, float ymin, float ymax,
           float zmin, float zmax)
    : cellSize(cellSize), xmin(xmin), ymin(ymin), zmin(zmin) {
  Nx = std::ceil((xmax - xmin) / cellSize);
  Ny = std::ceil((ymax - ymin) / cellSize);
  Nz = std::ceil((zmax - zmin) / cellSize);

  cells.resize(Nx * Ny * Nz);
}

int Grid::computeCellIndex(float x, float y, float z) const {
  int cx = (int)((x - xmin) / cellSize);
  int cy = (int)((y - ymin) / cellSize);
  int cz = (int)((z - zmin) / cellSize);

  if (cx < 0)
    cx = 0;
  if (cx >= Nx)
    cx = Nx - 1;
  if (cy < 0)
    cy = 0;
  if (cy >= Ny)
    cy = Ny - 1;
  if (cz < 0)
    cz = 0;
  if (cz >= Nz)
    cz = Nz - 1;

  return cx + cy * Nx + cz * Nx * Ny;
}

void Grid::build(const std::vector<Particle> &particles) {
  for (auto &c : cells)
    c.clear();

  for (int i = 0; i < particles.size(); ++i) {
    int cellID =
        computeCellIndex(particles[i].x, particles[i].y, particles[i].z);
    cells[cellID].push_back(i);
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
