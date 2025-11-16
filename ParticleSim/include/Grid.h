#ifndef GRID_H
#define GRID_H

#include "Particle.h"
#include <vector>

class Grid {
public:
  Grid(float cellSize, float xmin, float xmax, float ymin, float ymax,
       float zmin, float zmax);

  void build(const std::vector<Particle> &particles);

private:
  float cellSize;
  int Nx, Ny, Nz;

  float xmin, ymin, zmin;

  std::vector<std::vector<int>> cells; // cell → list of particle indices

  int computeCellIndex(float x, float y, float z) const;
};

#endif // GRID_H
