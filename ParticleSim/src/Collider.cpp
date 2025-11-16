#include "Collider.h"
#include "Config.h"

bool Collider::isOutOfBounds(float x, float y, float z) {
  return (x < Config::xmin || x > Config::xmax || y < Config::ymin ||
          y > Config::ymax || z < Config::zmin || z > Config::zmax);
}
