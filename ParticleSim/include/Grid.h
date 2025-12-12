#ifndef GRID_H
#define GRID_H

#include "Particle.h"
#include <vector>

class Grid {
public:
  Grid(float cellSize, float xmin, float xmax, float ymin, float ymax,
       float zmin, float zmax);

  void build(const std::vector<Particle> &particles);
  std::vector<std::vector<int>> getCells() const { return cells; }
  int getNx() const { return Nx; }
  int getNy() const { return Ny; }
  int getNz() const { return Nz; }

private:
  float cellSize;
  int Nx, Ny, Nz;

  float xmin, ymin, zmin;

  std::vector<std::vector<int>> cells;

  int computeCellIndex(float x, float y, float z) const;
};

#endif // GRID_H
