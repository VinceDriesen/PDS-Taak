#ifndef CONFIG_H
#define CONFIG_H

struct Config {
  constexpr static float xmin = 0;
  constexpr static float xmax = 10;

  constexpr static float ymin = 0;
  constexpr static float ymax = 10;

  constexpr static float zmin = 0;
  constexpr static float zmax = 10;

  constexpr static float gravity = -9.81;
  constexpr static float diameter = 0.5;
};

#endif // CONFIG_H
