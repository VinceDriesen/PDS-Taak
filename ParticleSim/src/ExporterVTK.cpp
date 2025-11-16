#include "ExporterVTK.h"
#include <fstream>

void ExporterVTK::saveVTKFile(const std::vector<Particle> &particles,
                              const std::string &folder, const int frame) {
  char filename[256];
  sprintf(filename, "%s/particles_%04d.vtk", folder.c_str(), frame);

  std::ofstream file(filename);
  file << "# vtk DataFile Version 3.0\n";
  file << "Particle simulation\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";
  file << "POINTS " << particles.size() << " float\n";

  for (const auto &p : particles) {
    file << p.x << " " << p.y << " " << p.z << "\n";
  }
}

void ExporterVTK::saveBoxVTK(const std::string &filename, float xmin,
                             float xmax, float ymin, float ymax, float zmin,
                             float zmax) {
  std::ofstream file(filename);
  file << "# vtk DataFile Version 3.0\n";
  file << "Bounding box\n";
  file << "ASCII\n";
  file << "DATASET POLYDATA\n";

  // 8 corners
  file << "POINTS 8 float\n";
  file << xmin << " " << ymin << " " << zmin << "\n"; // 0
  file << xmax << " " << ymin << " " << zmin << "\n"; // 1
  file << xmax << " " << ymax << " " << zmin << "\n"; // 2
  file << xmin << " " << ymax << " " << zmin << "\n"; // 3
  file << xmin << " " << ymin << " " << zmax << "\n"; // 4
  file << xmax << " " << ymin << " " << zmax << "\n"; // 5
  file << xmax << " " << ymax << " " << zmax << "\n"; // 6
  file << xmin << " " << ymax << " " << zmax << "\n"; // 7

  // LINES: 12 edges, elke lijn: npts + indices
  file << "LINES 12 " << 12 * 3 << "\n";
  file << "2 0 1\n";
  file << "2 1 2\n";
  file << "2 2 3\n";
  file << "2 3 0\n";

  file << "2 4 5\n";
  file << "2 5 6\n";
  file << "2 6 7\n";
  file << "2 7 4\n";

  file << "2 0 4\n";
  file << "2 1 5\n";
  file << "2 2 6\n";
  file << "2 3 7\n";

  file << "CELL_DATA 12\n";
  file << "SCALARS lineColor float 1\n";
  file << "LOOKUP_TABLE default\n";

  for (int i = 0; i < 12; i++) {
    file << i << "\n";
  }

  file.close();
}
