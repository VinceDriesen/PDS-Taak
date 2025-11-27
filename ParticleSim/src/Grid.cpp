#include "Grid.h"
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
