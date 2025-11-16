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
