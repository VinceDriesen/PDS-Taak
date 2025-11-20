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
  constexpr static float diameter = 0.1;

  constexpr static float smoothingRadius = diameter * 3.0f;
  constexpr static float particleMass = 10.0f;
  constexpr static float stiffness = 100.0f;
  constexpr static float viscosity = 5.0f;
  constexpr static float repulsionStiffness = 8000.0f;
  constexpr static float maxPressureForce = 800.0f;
  constexpr static float restDensity = 1000.0f;

  constexpr static float maxSpeed = 50.0f;
};

#endif // CONFIG_H
