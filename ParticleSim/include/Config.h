#ifndef CONFIG_H
#define CONFIG_H

struct Config {
  constexpr static float xmin = 0;
  constexpr static float xmax = 30;

  constexpr static float ymin = 0;
  constexpr static float ymax = 30;

  constexpr static float zmin = 0;
  constexpr static float zmax = 30;

  constexpr static float gravity = -9.81;
  constexpr static float diameter = 0.3;

  constexpr static float smoothingRadius = diameter * 3.0f;

  // AANPASSING 1: Massa iets realistischer t.o.v. volume en dichtheid
  constexpr static float particleMass = 14.0f;

  // AANPASSING 2: Stijfheid flink omhoog voor minder samendrukbaarheid
  // Dit zorgt voor meer "bounce" en minder energieverlies in compressie
  constexpr static float stiffness = 1000.0f; // Was 100.0f

  // AANPASSING 3: Viscositeit drastisch omlaag
  // Dit is de belangrijkste factor voor het "tot rust komen"
  constexpr static float viscosity = 0.3f; // Was 2.0f

  constexpr static float repulsionStiffness = 8000.0f;
  constexpr static float maxPressureForce = 800.0f;
  constexpr static float restDensity = 1000.0f;

  // AANPASSING 4: Max snelheid iets omhoog om pieken toe te laten
  constexpr static float maxSpeed = 60.0f;
};

#endif // CONFIG_H
