#ifndef EXPORTER_VTK_H
#define EXPORTER_VTK_H

#include "Particle.h"
#include <string>
#include <vector>

class ExporterVTK {
public:
  static void saveVTKFile(const std::vector<Particle> &particles,
                          const std::string &folder, const int frame);
};

#endif // EXPORTER_VTK_H
